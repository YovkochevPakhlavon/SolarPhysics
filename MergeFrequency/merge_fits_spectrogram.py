"""
merge_fits_spectrogram.py

This script provides tools to merge two spectrogram FITS files by their
frequency axes and visualize the result. It supports optional user-defined
frequency slicing and handles automatic frequency sorting.

Author: Christian Monstein (original)

Adapted and modularized by: Abraham-Alowonle Joseph-judah

Dependencies:
    - numpy
    - matplotlib
    - astropy

Usage:
    - Run as a standalone script or import functions into another module.
    - Optionally pass frequency range tuples for finer control.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits
from typing import Optional, Tuple



def read_fits_data(filepath: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reads a FITS file and extracts the data, frequency, and time axes.

    Parameters:
        filepath (str): Path to the FITS file.

    Returns:
        Tuple of (data, frequency array, time array).
    """
    with fits.open(filepath) as hdu:
        data = hdu[0].data.astype(np.float32)
        freqs = hdu[1].data['Frequency'][0]
        time = hdu[1].data['Time'][0]
    return data, freqs, time


def slice_by_freq_range(
    data: np.ndarray,
    freqs: np.ndarray,
    freq_range: Optional[Tuple[float, float]]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Optionally slices the data and frequency array to a specific frequency range.

    Parameters:
        data (np.ndarray): Spectrogram data array.
        freqs (np.ndarray): Frequency axis.
        freq_range (tuple or None): (min_freq, max_freq) in MHz.

    Returns:
        Sliced (data, freqs) if range is specified; otherwise original.
    """
    if freq_range is None:
        return data, freqs
    fmin, fmax = freq_range
    idx = np.where((freqs >= fmin) & (freqs <= fmax))[0]
    return data[idx, :], freqs[idx]


def merge_fits_by_frequency_auto(
    file1: str,
    file2: str,
    freq_range1: Optional[Tuple[float, float]] = None,
    freq_range2: Optional[Tuple[float, float]] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Merges two FITS files in frequency space. Supports optional slicing and automatic sorting.

    Parameters:
        file1 (str): Path to the first FITS file.
        file2 (str): Path to the second FITS file.
        freq_range1 (tuple or None): Frequency range to slice file1 (min_freq, max_freq).
        freq_range2 (tuple or None): Frequency range to slice file2 (min_freq, max_freq).

    Returns:
        Merged (data, frequency, time) tuple.
    """
    data1, freqs1, time1 = read_fits_data(file1)
    data2, freqs2, time2 = read_fits_data(file2)

    if not np.allclose(time1, time2, atol=1e-3):
        raise ValueError("Time axes of the two FITS files do not align.")

    data1, freqs1 = slice_by_freq_range(data1, freqs1, freq_range1)
    data2, freqs2 = slice_by_freq_range(data2, freqs2, freq_range2)

    data_merged = np.concatenate((data1, data2), axis=0)
    freqs_merged = np.concatenate((freqs1, freqs2))

    # Sort by frequency
    sort_idx = np.argsort(freqs_merged)
    data_merged = data_merged[sort_idx, :]
    freqs_merged = freqs_merged[sort_idx]

    return data_merged, freqs_merged, time1

def plot_spectrogram(
    data: np.ndarray,
    freqs: np.ndarray,
    time: np.ndarray,
    title: str = "Merged Spectrogram",
    vmin: float = -1,
    vmax: float = 40,
    save_path: Optional[str] = None
):
    """
    Plots and optionally saves the spectrogram of merged FITS data.

    Parameters:
        data (np.ndarray): Spectrogram data.
        freqs (np.ndarray): Frequency axis.
        time (np.ndarray): Time axis.
        title (str): Plot title.
        vmin (float): Minimum intensity for color scale.
        vmax (float): Maximum intensity for color scale.
        save_path (str or None): File path to save the figure (e.g., 'output.png').
                                 If None, the plot is shown but not saved.
    """
    data = data - data.mean(axis=1, keepdims=True)
    extent = (time[0], time[-1], freqs[-1], freqs[0])

    plt.figure(figsize=(12, 8))
    plt.imshow(data, aspect='auto', extent=extent, cmap=cm.plasma, vmin=vmin, vmax=vmax)
    plt.colorbar(label='Relative dB')
    plt.xlabel("Time [s]", fontsize=14)
    plt.ylabel("Frequency [MHz]", fontsize=14)
    plt.title(title, fontsize=16)
    plt.tick_params(labelsize=12)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to: {save_path}")
    else:
        plt.show()

# if __name__ == "__main__":
#     file1 = "BIR_20110809_080000_59.fit.gz"
#     file2 = "ALMATY_20110809_080000_59.fit.gz"

#     # Optional frequency ranges (in MHz)
#     freq_range1 = None
#     freq_range2 = None

#     data, freqs, time = merge_fits_by_frequency_auto(file1, file2, freq_range1, freq_range2)

#     # Plot and save
#     plot_spectrogram(
#         data, freqs, time,
#         title="Merged FITS Spectrogram (Auto-Aligned)",
#         save_path="merged_spectrogram.png"
#     )
