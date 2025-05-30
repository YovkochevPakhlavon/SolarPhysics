"""
goes_xray_flux_analysis.py

Created on Mon Aug 26 2024
Author: C. Moinstein

Author: Adapted from Christian Monstein’s code

Adapted and modularized by: Abraham-Alowonle Joseph-judah
Description:
------------
This module loads GOES-16 X-ray flux data from a netCDF file,
extracts specified variables, identifies peak and surrounding minima,
and optionally plots the flux time series with annotated peaks and minima.

Usage Example:
--------------
from goes_xray_flux_analysis import plot_goes_xray_flux

file_path = "path_to_your_file.nc"
plot_goes_xray_flux(file_path, num_vars=2, time_slice=(300, 800),
                   plot_range_dates=("2023-06-12 06:50:00", "2023-06-12 07:20:00"),
                   save_plot=True)

"""

import netCDF4 as nc
import numpy as np
import cftime
import matplotlib.pyplot as plt
from datetime import datetime


def plot_goes_xray_flux(nc_file_path, num_vars=2, time_slice=(300, 800),
                        plot_range_dates=None, make_plot=True, save_plot=False):
    """
    Load GOES-16 X-ray flux data from netCDF and plot with annotated peaks and minima.

    Parameters
    ----------
    nc_file_path : str
        Full path to the GOES netCDF (.nc) file containing X-ray flux data.
    num_vars : int, optional
        Number of energy level variables to plot. Default is 2.
    time_slice : tuple of int, optional
        Tuple (start_idx, end_idx) defining the index range of time data to analyze/plot.
        Default is (300, 800).
    plot_range_dates : tuple of str or None, optional
        Tuple of start and end date/time strings (e.g. "YYYY-MM-DD HH:MM:SS") to set x-axis limits.
        If None, axis limits are not set. Default is None.
    make_plot : bool, optional
        Whether to generate and display the plot. Default is True.
    save_plot : bool, optional
        Whether to save the plot as a PNG file. Default is False.

    Returns
    -------
    None

    Side Effects
    ------------
    - Prints file info and satellite platform.
    - Displays a plot (if make_plot=True).
    - Saves the plot PNG to the same folder as the input file (if save_plot=True).

    """
    # Open netCDF file
    ff = nc.Dataset(nc_file_path)
    
    # Convert time variable to datetime objects
    datetime0 = cftime.num2pydate(ff.variables["time"][:], ff["time"].units)
    
    print(f"Filename: {nc_file_path}")
    print(f"Start time in file [{ff['time'].units}]: {ff.variables['time'][0]}")
    print(f"Start and end times: {datetime0[time_slice[0]]}, {datetime0[time_slice[1]]}")
    
    platform = getattr(ff, "platform", "Unknown")
    print(f"Satellite: {platform}")
    
    # Variable names for X-ray flux channels
    var_name = ["xrsa_flux", "xrsb_flux"][:num_vars]
    
    # Parse plot range dates to datetime objects if provided
    if plot_range_dates is not None:
        t1 = datetime.strptime(plot_range_dates[0], "%Y-%m-%d %H:%M:%S")
        t2 = datetime.strptime(plot_range_dates[1], "%Y-%m-%d %H:%M:%S")
    else:
        t1 = None
        t2 = None
    
    if make_plot:
        chan_color = ["mediumorchid", "green"][:num_vars]
        
        plt.figure(figsize=(10, 6))
        
        for ii in range(num_vars):
            data = ff.variables[var_name[ii]][time_slice[0]:time_slice[1]]
            times = datetime0[time_slice[0]:time_slice[1]]
            
            plt.plot(times, data, linewidth=2, color=chan_color[ii], label=f"{platform} {var_name[ii]}")
            
            # Find max (peak) index and value
            max_idx = np.argmax(data)
            max_val = data[max_idx]
            plt.scatter(times[max_idx], max_val, color='red', s=100, zorder=5)
            plt.text(times[max_idx], max_val, f'{var_name[ii]} Max', fontsize=12, ha='center', color='red', va='bottom')
            plt.axvline(x=times[max_idx], color='red', linestyle='--', linewidth=1.5)
            
            # Find min before peak
            if max_idx > 0:
                min_before_idx = np.argmin(data[:max_idx])
                min_before_val = data[min_before_idx]
                plt.scatter(times[min_before_idx], min_before_val, color='blue', s=100, zorder=5)
                plt.text(times[min_before_idx], min_before_val, f'{var_name[ii]} Min Before Peak',
                         fontsize=12, ha='center', color='blue', va='top')
                plt.axvline(x=times[min_before_idx], color='blue', linestyle='--', linewidth=1.5)
            
            # Find min after peak
            if max_idx < len(data) - 1:
                min_after_idx = np.argmin(data[max_idx:])
                min_after_idx += max_idx  # adjust index relative to slice
                min_after_val = data[min_after_idx]
                plt.scatter(times[min_after_idx], min_after_val, color='purple', s=100, zorder=5)
                plt.text(times[min_after_idx], min_after_val, f'{var_name[ii]} Min After Peak',
                         fontsize=12, ha='center', color='purple', va='top')
                plt.axvline(x=times[min_after_idx], color='purple', linestyle='--', linewidth=1.5)
        
        plt.yscale("log")
        plt.xlabel("Time [UT]", fontsize=14)
        plt.ylabel(f"X-Ray Flux [{ff[var_name[0]].units}]", fontsize=14)
        plt.title(f"GOES-16 X-Ray Flux on {datetime0[0].date()}", fontsize=16)
        
        if t1 is not None and t2 is not None:
            plt.xlim(t1, t2)
        
        plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
        plt.legend(loc="lower left", prop={"size": 12})
        plt.tight_layout()
        
        if save_plot:
            output_filename = nc_file_path + '_peaks_and_minima_corrected.png'
            plt.savefig(output_filename, bbox_inches="tight", dpi=300)
            print(f"Plot saved to: {output_filename}")
        
        plt.show()
    
    print("Done.\n")


# from goes_xray_flux_analysis import plot_goes_xray_flux

# plot_goes_xray_flux(
#     "C:/Users/Joseph-Judah/Desktop/MY_SERIOUS_WORK/COSPAR_2024/WEEK_2/sci_xrsf-l2-avg1m_g18_d20230612_v2-2-0.nc",
#     num_vars=2,
#     time_slice=(300, 800),
#     plot_range_dates=("2023-06-12 06:50:00", "2023-06-12 07:20:00"),
#     make_plot=True,
#     save_plot=True,
# )