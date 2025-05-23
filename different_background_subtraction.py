"""
Module for processing FITS radio data files (Using different background subtraction techniques) and visualizing spectrograms with various enhancements.

Author: Christian Monstein (original script)
Adapted and modularized by: Abraham-Alowonle Joseph-judah
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits
from scipy.ndimage import gaussian_filter

def process_fits_spectrogram(
    file_path,
    output_dir='output',
    fmin=None,
    fmax=None,
    tmin=None,
    tmax=None,
    smoothing_sigma=(0, 1),
    constant_offset=10,
    reference_time=None
):
    """
    Process a FITS (.fit.gz) file to generate spectrograms using various background subtraction methods.

    Parameters
    ----------
    file_path : str
        Path to the FITS (.fit.gz) file.
    output_dir : str, optional
        Directory to save the output images. Default is 'output'.
    fmin : float, optional
        Minimum frequency in MHz for plotting. If None, uses the minimum frequency from the data.
    fmax : float, optional
        Maximum frequency in MHz for plotting. If None, uses the maximum frequency from the data.
    tmin : str, optional
        Start time in 'HH:MM:SS' format for plotting. If None, uses the start time from the data.
    tmax : str, optional
        End time in 'HH:MM:SS' format for plotting. If None, uses the end time from the data.
    smoothing_sigma : tuple of float, optional
        Standard deviations for Gaussian kernel in (frequency, time) directions. Default is (0, 1).
    constant_offset : float, optional
        Constant value to subtract from the data for background subtraction. Default is 10.
    reference_time : float, optional
        Reference time in seconds after the start time to use for individual spectrum subtraction.
        If None, uses the midpoint of the time axis.

    Returns
    -------
    dict
        Dictionary containing the processed data arrays:
        - 'raw': Original data.
        - 'smoothed': Data after Gaussian smoothing.
        - 'constant_subtracted': Data after constant background subtraction.
        - 'average_subtracted': Data after average spectrum subtraction.
        - 'reference_subtracted': Data after individual spectrum subtraction.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Open FITS file
    with fits.open(file_path) as hdu:
        data = hdu[0].data.astype(np.float32)
        freqs = hdu[1].data['Frequency'][0]
        time = hdu[1].data['Time'][0]
        header = hdu[0].header

    # Set frequency and time ranges
    fmin = fmin if fmin is not None else freqs.min()
    fmax = fmax if fmax is not None else freqs.max()
    tmin_sec = 0
    tmax_sec = time[-1]
    if tmin:
        h, m, s = map(int, tmin.split(':'))
        tmin_sec = h * 3600 + m * 60 + s
    if tmax:
        h, m, s = map(int, tmax.split(':'))
        tmax_sec = h * 3600 + m * 60 + s

    # Compute extent for plotting
    extent = (time[0], time[-1], freqs[-1], freqs[0])

    # Raw data plot
    plt.figure(figsize=(12, 8))
    plt.imshow(data, aspect='auto', extent=extent, cmap=cm.inferno)
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [MHz]')
    plt.title('Raw Spectrogram')
    plt.colorbar(label='Intensity')
    raw_path = os.path.join(output_dir, 'raw_spectrogram.png')
    plt.savefig(raw_path)
    plt.close()

    # Gaussian smoothing
    smoothed_data = gaussian_filter(data, sigma=smoothing_sigma)
    plt.figure(figsize=(12, 8))
    plt.imshow(smoothed_data, aspect='auto', extent=extent, cmap=cm.inferno)
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [MHz]')
    plt.title('Smoothed Spectrogram')
    plt.colorbar(label='Intensity')
    smoothed_path = os.path.join(output_dir, 'smoothed_spectrogram.png')
    plt.savefig(smoothed_path)
    plt.close()

    # Constant background subtraction
    constant_subtracted = data - np.mean(data) - constant_offset
    plt.figure(figsize=(12, 8))
    plt.imshow(constant_subtracted, aspect='auto', extent=extent, cmap=cm.inferno)
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [MHz]')
    plt.title(f'Constant Background Subtracted (Offset={constant_offset})')
    plt.colorbar(label='Intensity')
    constant_path = os.path.join(output_dir, 'constant_subtracted_spectrogram.png')
    plt.savefig(constant_path)
    plt.close()

    # Average spectrum subtraction
    avg_spectrum = data.mean(axis=1, keepdims=True)
    average_subtracted = data - avg_spectrum
    plt.figure(figsize=(12, 8))
    plt.imshow(average_subtracted, aspect='auto', extent=extent, cmap=cm.inferno)
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [MHz]')
    plt.title('Average Spectrum Subtracted')
    plt.colorbar(label='Intensity')
    average_path = os.path.join(output_dir, 'average_subtracted_spectrogram.png')
    plt.savefig(average_path)
    plt.close()

    # Individual spectrum subtraction
    if reference_time is None:
        reference_time = time[len(time) // 2]
    time_index = np.argmin(np.abs(time - reference_time))
    reference_spectrum = data[:, time_index]
    reference_subtracted = data - reference_spectrum[:, np.newaxis]
    plt.figure(figsize=(12, 8))
    plt.imshow(reference_subtracted, aspect='auto', extent=extent, cmap=cm.inferno)
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [MHz]')
    plt.title(f'Individual Spectrum Subtracted at {reference_time:.2f} s')
    plt.colorbar(label='Intensity')
    reference_path = os.path.join(output_dir, 'reference_subtracted_spectrogram.png')
    plt.savefig(reference_path)
    plt.close()

    return {
        'raw': data,
        'smoothed': smoothed_data,
        'constant_subtracted': constant_subtracted,
        'average_subtracted': average_subtracted,
        'reference_subtracted': reference_subtracted
    }
