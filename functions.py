import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import detrend, welch
from scipy import stats
import geopandas as gpd
from shapely.geometry import Point
import geodatasets

ID_VARS = ['RGIId', 'REGION', 'O2Region', 'CenLon', 'CenLat', 'Area', 'WGMS_ID']


def get_climate_oscillation_periods():
    """
    Returns user-defined return periods for climate oscillations.
    These values can be easily modified in this single function.

    Returns:
        tuple: A tuple containing:
            - np.array: Array of ENSO return periods.
            - np.array: Array of PDO return periods.
    """
    enso_return_periods = np.array([2, 2.5, 3.5, 5.5, 3, 3, 5, 3.5, 4, 5, 4.5])
    pdo_return_periods = np.array([21, 16.5])
    return enso_return_periods, pdo_return_periods


def load_glacier_data(file_path):
    """
    Loads glacier mass balance data from a CSV file and transforms it into a long format.

    Args:
        file_path (str): The path to the CSV file.

    Returns:
        tuple: A tuple containing:
            - pandas.DataFrame: The raw glacier data.
            - pandas.DataFrame: The glacier data in long format with 'Year' as numeric.
    """
    try:
        glacier_raw = pd.read_csv(file_path, na_values=['-'])
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}. Please ensure the path is correct.")
        return None, None

    value_vars = [col for col in glacier_raw.columns if col not in ID_VARS]
    glacier_long = glacier_raw.melt(id_vars=ID_VARS,
                                     value_vars=value_vars,
                                     var_name='Year',
                                     value_name='Mass_Balance')
    glacier_long['Year'] = pd.to_numeric(glacier_long['Year'], errors='coerce')
    return glacier_raw, glacier_long

def plot_average_mass_balance(glacier_long_df, line_color='midnightblue'):
    """
    Plots the average annual mass balance across all glaciers.

    Args:
        glacier_long_df (pd.DataFrame): The glacier data in long format.
    """
    # Calculate the average mass balance per year, ignoring NaN values
    average_mass_balance = glacier_long_df.groupby('Year')['Mass_Balance'].mean().reset_index()

    plt.figure(figsize=(10, 6))
    sns.lineplot(data=average_mass_balance, x = 'Year', y = 'Mass_Balance', color = line_color)
    plt.title('Average Annual Mass Balance Across All Glaciers')
    plt.xlabel('Year')
    plt.ylabel('Mass Balance')
    plt.grid(True, linestyle='-', alpha=0.7)
    plt.show()

def prepare_mass_balance_for_psd(glacier_long_df):
    """
    Prepares the average mass balance time series for PSD analysis.

    Args:
        glacier_long_df (pd.DataFrame): The glacier data in long format.

    Returns:
        np.array: Detrended mass balance values.
    """
    avg_mass_balance_series = glacier_long_df.groupby('Year')['Mass_Balance'].mean().sort_index()
    mass_balance_values = avg_mass_balance_series.interpolate(method='linear').dropna().values

    if len(mass_balance_values) < 20:
        print(f"Warning: Too few valid data points ({len(mass_balance_values)} points) for meaningful PSD. At least 20-30 points are recommended.")
        return None

    return detrend(mass_balance_values)

def calculate_psd_welch(detrended_data):
    """
    Calculates Power Spectral Density (PSD) using Welch's method.

    Args:
        detrended_data (np.array): Detrended time series data.

    Returns:
        tuple: A tuple containing:
            - np.array: Frequencies for the PSD.
            - np.array: PSD values.
            - np.array: Periods corresponding to PSD values.
    """
    fs = 1.0  # Sample frequency (1 sample per year)
    nperseg = min(len(detrended_data), 64) # Max nperseg = 64
    if nperseg % 2 != 0: nperseg -= 1 # Ensure nperseg is even

    noverlap_value = nperseg // 2  # 50% overlap for smoothness
    nfft_value = max(128, nperseg) # Padding for visual smoothness

    f, Pxx = welch(detrended_data, fs=fs, nperseg=nperseg, noverlap=noverlap_value,
                    window='hann', detrend='linear', nfft=nfft_value)

    # Filter 0 frequency (DC component) and negative frequencies
    valid_indices = (f > 0) & (f <= 0.5)
    psd_frequencies = f[valid_indices]
    psd_values = Pxx[valid_indices]
    psd_periods = 1 / psd_frequencies
    return psd_frequencies, psd_values, psd_periods

def calculate_kde_scaling(psd_values, psd_periods, min_period, max_period, kde_y_density):
    """
    Helper function to calculate the scaling factor for KDE plots based on PSD values.

    Args:
        psd_values (np.array): Array of Power Spectral Density values.
        psd_periods (np.array): Array of periods corresponding to PSD values.
        min_period (float): Minimum period for the scaling range.
        max_period (float): Maximum period for the scaling range.
        kde_y_density (np.array): Array of KDE density values.

    Returns:
        float: The calculated scaling factor.
    """
    range_psd_values = psd_values[(psd_periods >= min_period) & (psd_periods <= max_period)]
    target_value_psd = 0.05 # Default or base value if no peaks
    if len(range_psd_values) > 0:
        target_value_psd = np.max(range_psd_values)
    
    scaling_factor = target_value_psd / np.max(kde_y_density) if np.max(kde_y_density) > 0 else 1.0
    return scaling_factor * 0.9 # Adjust this multiplier for visual height


def plot_psd_base(ax, psd_periods, psd_values, title_prefix=""):
    """
    Plots the base Power Spectral Density (PSD) line.

    Args:
        ax (matplotlib.axes.Axes): The axes object to plot on.
        psd_periods (np.array): Periods from PSD.
        psd_values (np.array): PSD values.
        title_prefix (str): Prefix for the plot title (e.g., "Average Glacier").
    """
    sns.lineplot(x=psd_periods, y=psd_values, markersize=4, color='midnightblue',
                 label=f'{title_prefix} Mass Balance (PSD - Welch)', ax=ax)

    ax.set_title(f'Spectral Power Density of {title_prefix} Mass Balance')
    ax.set_xlabel('Period (Years)')
    ax.set_ylabel('Power Spectral Density (unit²/Hz)')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_ylim(bottom=0)
    ax.set_xscale('log')

    plot_xmin_periods = 1.0
    plot_xmax_periods = 100
    ax.set_xlim(plot_xmin_periods, plot_xmax_periods)

    # Set x-axis ticks
    desired_periods_on_x_axis = [1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100]
    valid_periods_for_ticks = [p for p in desired_periods_on_x_axis if plot_xmin_periods <= p <= plot_xmax_periods]
    ax.set_xticks(valid_periods_for_ticks)
    ax.set_xticklabels([f'{p:.0f}' for p in valid_periods_for_ticks])


def add_climate_overlays(ax, psd_values, psd_periods):
    """
    Adds KDEs for climate oscillations and reference vertical lines/shaded zones
    to an existing PSD plot.

    Args:
        ax (matplotlib.axes.Axes): The axes object to plot on.
        psd_values (np.array): PSD values (needed for KDE scaling).
        psd_periods (np.array): Periods from PSD (needed for KDE scaling).
    """
    enso_return_periods_user, pdo_return_periods_user = get_climate_oscillation_periods()

    # --- KDE for ENSO return periods ---
    kde_enso = stats.gaussian_kde(enso_return_periods_user)
    kde_enso_x_values_periods = np.linspace(1.0, 7.5, 500)
    kde_enso_y_density = kde_enso(kde_enso_x_values_periods)
    kde_enso_scaling_factor = calculate_kde_scaling(psd_values, psd_periods, 2, 7, kde_enso_y_density)
    kde_enso_y_density_scaled = kde_enso_y_density * kde_enso_scaling_factor

    sns.lineplot(x=kde_enso_x_values_periods, y=kde_enso_y_density_scaled, color='chocolate', linestyle='--', linewidth=2,
                 label='ENSO - KDE', ax=ax)
    ax.fill_between(kde_enso_x_values_periods, kde_enso_y_density_scaled, color='chocolate', alpha=0.1)

    # --- KDE for PDO return periods ---
    kde_pdo = stats.gaussian_kde(pdo_return_periods_user)
    kde_pdo_x_values_periods = np.linspace(15, 25, 500)
    kde_pdo_y_density = kde_pdo(kde_pdo_x_values_periods)
    kde_pdo_scaling_factor = calculate_kde_scaling(psd_values, psd_periods, 20, 30, kde_pdo_y_density)
    kde_pdo_y_density_scaled = kde_pdo_y_density * kde_pdo_scaling_factor

    sns.lineplot(x=kde_pdo_x_values_periods, y=kde_pdo_y_density_scaled, color='purple', linestyle='-.', linewidth=2,
                 label='PDO - KDE', ax=ax)
    ax.fill_between(kde_pdo_x_values_periods, kde_pdo_y_density_scaled, color='purple', alpha=0.1)

    # Add reference ranges
    ax.axvspan(2.0, 7.0, color='chocolate', alpha=0.2, label='ENSO Range (2-7 years)')
    #ax.axvline(x=2.0, color='green', linestyle=':', label='2-year cycle')
    #ax.axvline(x=5.0, color='darkgreen', linestyle=':', label='5-year cycle')

    min_period_pdo_ref = 20
    max_period_pdo_ref = 30
    ax.axvspan(min_period_pdo_ref, max_period_pdo_ref, color='purple', alpha=0.1, label='PDO Reference Range (20-30 years)')
    #ax.axvline(x=20.0, color='red', linestyle='--', label='20-year cycle')

    # Update title with full context if overlays are included
    ax.set_title(f'{ax.get_title()} & Climate Oscillation Return Periods')


def analyze_prominent_peaks(psd_frequencies, psd_values, psd_periods):
    """
    Analyzes and prints information about prominent peaks in specified period ranges.
    Variables are initialized to np.nan to prevent NameError if no peak is found in a range,
    and a try-except block is used for robustness during calculation.

    Args:
        psd_frequencies (np.array): Frequencies from PSD.
        psd_values (np.array): PSD values.
        psd_periods (np.array): Periods from PSD.
    """
    # Initialize variables to NaN (Not a Number) to ensure they are always defined
    peak_psd_pdo = np.nan
    peak_index_psd_pdo = np.nan
    peak_frequency_psd_pdo = np.nan
    peak_period_pdo = np.nan

    peak_psd_enso = np.nan
    peak_index_psd_enso = np.nan
    peak_frequency_psd_enso = np.nan
    peak_period_enso = np.nan


    # PDO range analysis
    min_period_pdo_analysis = 20
    max_period_pdo_analysis = 30
    min_freq_pdo_analysis = 1 / max_period_pdo_analysis
    max_freq_pdo_analysis = 1 / min_period_pdo_analysis

    pdo_indices_psd = np.where((psd_frequencies >= min_freq_pdo_analysis) & (psd_frequencies <= max_freq_pdo_analysis))[0]

    if len(pdo_indices_psd) > 0:
        try:
            peak_psd_pdo = np.max(psd_values[pdo_indices_psd])
            peak_index_psd_pdo = pdo_indices_psd[np.argmax(psd_values[pdo_indices_psd])]
            peak_frequency_psd_pdo = psd_frequencies[peak_index_psd_pdo]
            peak_period_pdo = 1 / peak_frequency_psd_pdo

            print(f"\n--- Results for Average Glacier (PDO Range: {min_period_pdo_analysis}-{max_period_pdo_analysis} years) ---")
            print(f"  Highest peak in PDO range: Period of {peak_period_pdo:.1f} years")
            print(f"  Frequency: {peak_frequency_psd_pdo:.3f} cycles/year")
            print(f"  Power Spectral Density (PSD): {peak_psd_pdo:.4f}")
        except Exception as e:
            print(f"\n--- Error calculating PDO peak in range {min_period_pdo_analysis}-{max_period_pdo_analysis} years: {e} ---")
            print(f"  No prominent peak details can be displayed for PDO range due to an error.")
    else:
        print(f"\nNo prominent peak found in the {min_period_pdo_analysis}-{max_period_pdo_analysis} year range for the average glacier.")

    # ENSO range analysis
    min_period_enso_analysis = 2
    max_period_enso_analysis = 7
    min_freq_enso_analysis = 1 / max_period_enso_analysis
    max_freq_enso_analysis = 1 / min_period_enso_analysis

    enso_indices_psd = np.where((psd_frequencies >= min_freq_enso_analysis) & (psd_frequencies <= max_freq_enso_analysis))[0]

    if len(enso_indices_psd) > 0:
        try:
            peak_psd_enso = np.max(psd_values[enso_indices_psd])
            peak_index_psd_enso = enso_indices_psd[np.argmax(psd_values[enso_indices_psd])]
            peak_frequency_psd_enso = psd_frequencies[peak_index_psd_enso]
            # CORRECTED LINE: Used peak_frequency_psd_enso instead of peak_frequency_enso
            peak_period_enso = 1 / peak_frequency_psd_enso

            print(f"\n--- Results for Average Glacier (ENSO Range: {min_period_enso_analysis}-{max_period_enso_analysis} years) ---")
            print(f"  Highest peak in ENSO range: Period of {peak_period_enso:.1f} years")
            print(f"  Frequency: {peak_frequency_psd_enso:.3f} cycles/year")
            print(f"  Power Spectral Density (PSD): {peak_psd_enso:.4f}")
        except Exception as e:
            print(f"\n--- Error calculating ENSO peak in range {min_period_enso_analysis}-{max_period_enso_analysis} years: {e} ---")
            print(f"  No prominent peak details can be displayed for ENSO range due to an error.")
    else:
        print(f"\nNo prominent peak found in the {min_period_enso_analysis}-{max_period_enso_analysis} year range for the average glacier.")

def plot_spectral_mass_balance(glacier_long_df, include_overlays=True, title_prefix="Average Glacier"):
    """
    Performs Power Spectral Density (PSD) analysis on the glacier mass balance,
    plots the results, optionally adds climate oscillation KDEs and reference zones,
    and analyzes prominent peaks.

    Args:
        glacier_long_df (pd.DataFrame): The glacier data in long format.
        include_overlays (bool): If True, overlays climate oscillation KDEs and reference zones.
        title_prefix (str): A string to prefix the plot title and legend label (e.g., "Average Glacier").
    """
    print(f"\n--- Starting Spectral Power Density (PSD) Analysis for {title_prefix} ---")

    detrended_mass_balance = prepare_mass_balance_for_psd(glacier_long_df)
    if detrended_mass_balance is None:
        return # Exit if data is insufficient

    psd_frequencies, psd_values, psd_periods = calculate_psd_welch(detrended_mass_balance)

    plt.figure(figsize=(10, 6))
    ax = plt.gca()

    plot_psd_base(ax, psd_periods, psd_values, title_prefix=title_prefix)

    if include_overlays:
        add_climate_overlays(ax, psd_values, psd_periods)

    plt.legend(loc='upper right', framealpha=0.9)
    plt.tight_layout()
    plt.show()

    # This function call performs the ENSO and PDO peak analysis and prints to console.
    analyze_prominent_peaks(psd_frequencies, psd_values, psd_periods)
    #print(interpretation_text)


def plot_glacier_locations(glacier_long_df):
    """
    Plots the geographical locations of the glaciers colored by region.

    Args:
        glacier_long_df (pd.DataFrame): The glacier data in long format.
    """
    # Get unique glacier entries for plotting to avoid plotting the same glacier multiple times
    unique_glaciers = glacier_long_df.drop_duplicates(subset=['RGIId'])

    # Create a GeoDataFrame
    geometry = [Point(xy) for xy in zip(unique_glaciers['CenLon'], unique_glaciers['CenLat'])]
    gdf = gpd.GeoDataFrame(unique_glaciers, geometry=geometry, crs="EPSG:4326")

    # Load landmasses
    land_path = geodatasets.get_path('naturalearth.land')
    land = gpd.read_file(land_path)

    # Plotting
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    land.plot(ax=ax, color='lightgray', edgecolor='white')
    
    # Use a scatter plot for glaciers to allow more customization than default gdf.plot
    sns.scatterplot(
        x=gdf.geometry.x,
        y=gdf.geometry.y,
        color='olive',  # Stel hier uw gewenste kleur in
        edgecolor='olive',
        size=gdf['Area'],
        sizes=(20, 200),
        legend=False,
        ax=ax,
        #alpha=0.7
    )

    # Finishing touches
    ax.set_title('Glacier Locations', fontsize=16)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    x_ticks = ax.get_xticks()
    ax.set_xticks(x_ticks) # <--- Voeg deze regel toe
    ax.set_xticklabels([f'{int(tick)}°' for tick in x_ticks])
    y_ticks = ax.get_yticks()
    ax.set_yticks(y_ticks) # <--- Voeg deze regel toe
    ax.set_yticklabels([f'{int(tick)}°' for tick in y_ticks])
    plt.tight_layout()
    plt.show()