"""
Module for processing FITS radio data files and visualizing spectrograms with various enhancements.

Author: Christian Monstein (original script)
Adapted and modularized by: Abraham-Alowonle Joseph-judah
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage import gaussian_filter
import pylab as mplot
from astropy.io import fits

def process_fit_gz_file(filepath, fmin, fmax, tmin, tmax):
    """
    Process a .fit.gz radio data file, generate and save spectrograms with raw, smoothed, and time-formatted axes.

    Parameters:
    -----------
    filepath : str
        Path to the .fit.gz FITS file containing radio spectrogram data.
    fmin : float
        Minimum frequency to display in the plots (MHz).
    fmax : float
        Maximum frequency to display in the plots (MHz).
    tmin : str
        Start time for zoomed-in plots, in 'HH:MM:SS' format (Universal Time).
    tmax : str
        End time for zoomed-in plots, in 'HH:MM:SS' format (Universal Time).

    Returns:
    --------
    None
        Saves three PNG files:
        - '2D_raw.png': Raw spectrogram of the entire frequency and time range.
        - '2D_smooth.png': Gaussian-smoothed spectrogram zoomed in time and frequency.
        - '2D_smooth_with_axes.png': Smoothed plot with axes in HH:MM:SS and frequency in MHz.
    """

    hdu = fits.open(filepath)
    data = hdu[0].data.astype(np.float32) / 255.0 * 2500.0 / 25.4
    rel_dB = data - np.min(data)
    rel_dB = rel_dB - rel_dB.mean(axis=1, keepdims=True)
    rel_dB = rel_dB.clip(-2, 12)

    freqs = hdu[1].data['Frequency'][0]
    time = hdu[1].data['Time'][0]
    extent = (time[0], time[-1], freqs[-1], freqs[0])

    # Raw image plot
    plt.figure(figsize=(15, 8))
    plt.imshow(rel_dB, aspect='auto', extent=extent, cmap=cm.hot)
    plt.tick_params(labelsize=14)
    plt.xlabel('Relative time [s]', fontsize=15)
    plt.ylabel('Frequency [MHz]', fontsize=15)
    plt.title('Raw FIT data visualization', fontsize=15)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('dB above background', fontsize=15)
    plt.savefig('2D_raw.png')

    # Smoothed image plot
    plt.figure(figsize=(15, 8))
    blur = gaussian_filter(rel_dB, sigma=[0, 1], mode='constant')
    plt.imshow(blur, aspect='auto', extent=extent, cmap=cm.hot_r,
               norm=plt.Normalize(vmin=-2, vmax=15))
    plt.xlim((400, 900))
    plt.ylim((fmin, fmax))
    plt.tick_params(labelsize=14)
    plt.xlabel('Relative time [s]', fontsize=15)
    plt.ylabel('Frequency [MHz]', fontsize=15)
    plt.title('Gaussian smoothed FIT data, zoomed view', fontsize=15)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('dB above background', fontsize=15)
    plt.savefig('2D_smooth.png')

    # Time corrected plot with HH:MM:SS labels
    time_obs = hdu[0].header['TIME-OBS']
    dt = hdu[0].header['CDELT1']
    hdu.close()

    hh, mm, ss = map(float, time_obs.split(':'))
    start_time = hh * 3600 + mm * 60 + ss
    ut = time + start_time

    secmin = sum(x * int(t) for x, t in zip([3600, 60, 1], tmin.split(':'))) - ut[0]
    secmax = sum(x * int(t) for x, t in zip([3600, 60, 1], tmax.split(':'))) - ut[0]
    xmin, xmax = (secmin + 1.0) / dt, (secmax + 1.0) / dt

    klow, khigh = next(i for i, f in enumerate(freqs) if f > fmin), next(i for i, f in enumerate(freqs) if f > fmax)

    mplot.figure(figsize=(15, 10))
    ax = mplot.subplot(111)
    mplot.axis([xmin, xmax, klow, khigh])
    mplot.xticks(np.arange(xmin, xmax + 1, (xmax - xmin) / 6))

    vmin, vmax = 0, 12
    mplot.imshow(blur, aspect='auto', cmap=cm.inferno, norm=plt.Normalize(vmin, vmax))
    cb1 = mplot.colorbar(orientation='vertical', shrink=0.99)
    cb1.set_label('dB above background', fontsize=15)
    cb1.ax.tick_params(labelsize=15)

    ahms = ['{:02d}:{:02d}:{:02d}'.format(
        int(ut[int(a)] / 3600),
        int((ut[int(a)] % 3600) / 60),
        int(ut[int(a)] % 60)
    ) for a in mplot.xticks()[0]]

    yticks = mplot.yticks()[0]
    ylabels = ['{:.2f}'.format(freqs[int(y)]) for y in yticks[:-1]]

    ax.set_xticklabels(ahms, fontsize=15)
    ax.set_yticklabels(ylabels, fontsize=15)
    mplot.ylabel('Observed frequency [MHz]', fontsize=12)
    mplot.xlabel('Time [UT]', fontsize=18)
    mplot.title(f"File: {filepath}", fontsize=20)
    mplot.savefig('2D_smooth_with_axes.png', bbox_inches='tight')
