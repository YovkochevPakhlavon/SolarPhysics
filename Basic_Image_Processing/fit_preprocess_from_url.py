"""
Download a FITS file from a given URL and perform different background subtraction
methods on the spectrogram data. Saves each processed spectrogram as PNG images.

Background subtraction methods implemented:
1) Subtract a constant value (mean + offset)
2) Subtract average spectrum (mean over time)
3) Subtract an individual spectrum at a specified time marker

Usage:
    Call `process_fit_file(url)` directly from your Jupyter notebook cell.

Author: Adapted from Christian Monstein’s code

Adapted and modularized by: Abraham-Alowonle Joseph-judah
"""


import urllib.request
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

def plot_spectrogram(data, freqs, time, cmap, vmin, vmax, title, filename):
    """
    Plot and save a spectrogram image.

    Parameters
    ----------
    data : 2D np.ndarray
        Spectrogram data to plot.
    freqs : 1D np.ndarray
        Frequency axis values.
    time : 1D np.ndarray
        Time axis values.
    cmap : matplotlib colormap
        Colormap to use for the image.
    vmin : float
        Minimum data value for colormap scaling.
    vmax : float
        Maximum data value for colormap scaling.
    title : str
        Title for the plot.
    filename : str
        File name to save the plot image.
    """
    plt.figure(figsize=(12, 8))
    extent = (time[0], time[-1], freqs[-1], freqs[0])
    plt.imshow(data, aspect='auto', extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.tick_params(labelsize=14)
    plt.xlabel('Time [s]', fontsize=15)
    plt.ylabel('Observed frequency [MHz]', fontsize=15)
    plt.title(title, fontsize=15)
    plt.savefig(filename)
    plt.show()
    plt.close()
    print(f"Saved plot: {filename}")

def process_fit_file(url, offset=10, timemarker_factor=4):
    """
    Download a FITS file from the given URL and apply different background subtraction
    techniques. Saves the results as PNG files and shows the plots inline.

    Parameters
    ----------
    url : str
        URL of the FITS file to download.
    offset : float, optional
        Constant offset for the first background subtraction (default is 10).
    timemarker_factor : int, optional
        Factor to compute the time index for individual spectrum subtraction
        as timemarker = timemarker_factor * 600 seconds (default is 4).
    """
    print(f"Downloading FIT file from: {url}")
    filename = url.rsplit('/', 1)[1]
    urllib.request.urlretrieve(url, filename)

    # Load FITS file
    with fits.open(filename) as hdu:
        data = hdu[0].data.astype(np.float32)
        freqs = hdu[1].data['Frequency'][0]
        time = hdu[1].data['Time'][0]

    # 1) Constant value subtraction (mean + offset)
    data_const_subtracted = data - np.mean(data) - offset
    plot_spectrogram(data_const_subtracted, freqs, time, cm.inferno, -15, 30,
                     f'Constant value ({offset}) subtracted', filename + '_const_subtracted.png')

    # 2) Average spectrum subtraction (mean over time axis)
    avg_spectrum = data.mean(axis=1, keepdims=True)
    data_avg_subtracted = data - avg_spectrum
    plot_spectrogram(data_avg_subtracted, freqs, time, cm.plasma, 3, 40,
                     'Average spectrum subtracted', filename + '_avg_subtracted.png')

    # 3) Individual spectrum subtraction at timemarker
    T = 600  # seconds interval
    timemarker = timemarker_factor * T
    if timemarker >= data.shape[1]:
        timemarker = data.shape[1] // 2  # fallback if out of range
    ref_spectrum = data[:, timemarker]
    data_indiv_subtracted = (data.T - ref_spectrum).T + 10
    plot_spectrogram(data_indiv_subtracted, freqs, time, cm.plasma, 10, 40,
                     f'Individual spectrum subtracted at T={timemarker}', filename + '_indiv_subtracted.png')

# url = "https://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/2021/08/27/EGYPT-Alexandria_20210827_104508_02.fit.gz"
# process_fit_file(url)