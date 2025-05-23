"""
Callisto CME Analyzer

This script reads one or more FITS files containing radio dynamic spectra
(from Callisto instruments), plots the data, and allows the user to select
a frequency drift (typically of a CME). The script then calculates:

- Electron density
- Solar radial distance (Newkirk model)
- Frequency drift rate (df/dt)
- Velocity of the emission front

Author: Christian Monstein, ETH Zurich


Adapted and modularized by: Abraham-Alowonle Joseph-judah
Date: 2025-05-20

Usage:
    In a script or Jupyter notebook, call:
        run_cme_analysis(
            path='data/',
            file_pattern='BLENSW_20180330*_01.fit.gz',
            harmonic=2,
            newkirk=1.8,
            zoom=[700, 1700, 30, 70],
            vmin=-5,
            vmax=25
        )

Left click to mark points. Right click to finish analysis.
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits

# Constants
RSUN_KM = 695700.0  # Solar radius in kilometers
ELECTRON_CONSTANT = 8.977e-3  # MHz * sqrt(cm^-3)


def run_cme_analysis(path, file_pattern, harmonic=2, newkirk=1.8,
                     zoom=[700, 1700, 30, 70], vmin=-5, vmax=25):
    """
    Main function to process FITS files and perform interactive CME speed analysis.

    Parameters
    ----------
    path : str
        Directory containing FITS files.
    file_pattern : str
        Wildcard pattern for FITS files.
    harmonic : int, optional
        Harmonic mode (1 = fundamental, 2 = 1st harmonic).
    newkirk : float, optional
        Multiplier for Newkirk density model.
    zoom : list of float, optional
        Zoom window [start_time, stop_time, min_freq, max_freq].
    vmin : float, optional
        Minimum value for dynamic spectrum color scaling.
    vmax : float, optional
        Maximum value for dynamic spectrum color scaling.
    """

    paths = glob.glob(path + file_pattern)
    if not paths:
        raise FileNotFoundError("No FITS files found matching pattern.")

    print(f"\nLoaded {len(paths)} FITS file(s).")

    def read_fits_data(f):
        with fits.open(f) as hdu:
            return hdu[0].data.astype(np.float32)

    # Load metadata from first FITS
    with fits.open(paths[0]) as hdu:
        header = hdu[0].header
        title = header['CONTENT']
        T0 = header['TIME-OBS']
        dt = header['CDELT1']
        time_axis = hdu[1].data[0][0]
        frequency = hdu[1].data[0][1]

    h, m, s = map(float, T0.split(":"))
    t_start = h * 3600 + m * 60 + s

    # Load data
    data = np.hstack([read_fits_data(f) for f in paths])
    time_s = np.arange(data.shape[1]) * dt
    time_hours = (t_start + time_s) / 3600.0
    dflat = data - np.mean(data, axis=1, keepdims=True)

    # Plot dynamic spectrum
    extent = (time_s[0], time_s[-1], frequency[-1], frequency[0])
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(dflat, extent=extent, aspect='auto',
                   cmap=cm.CMRmap, norm=plt.Normalize(vmin, vmax))
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(f"Time [s] after {T0} UT", fontsize=12)
    ax.set_ylabel("Frequency [MHz]", fontsize=12)
    ax.axis(zoom)
    ax.tick_params(labelsize=12)

    # Interactive callback setup
    coords = []

    def onclick(event):
        if event.button == 1:
            # Left-click: save coordinate
            coords.append((event.xdata, event.ydata))
            ax.plot(event.xdata, event.ydata, 'wo')
            fig.canvas.draw()
        elif event.button == 3:
            # Right-click: finalize and analyze
            fig.canvas.mpl_disconnect(cid)
            plt.savefig(paths[0] + '.spectrum.png')
            compute_and_plot(coords, frequency, harmonic, newkirk, paths[0])

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    print("Left click to mark points. Right click to analyze.")


def compute_and_plot(coords, frequency, harmonic, newkirk, output_prefix):
    """Compute parameters from selected points and generate plots + results."""

    time, freq = zip(*coords)
    time = np.array(time)
    freq = np.array(freq)

    Ne = (freq / (harmonic * ELECTRON_CONSTANT))**2
    rs = 4.32 / np.log10(Ne / (newkirk * 4.2e4))

    dfdt = np.abs(np.diff(freq) / np.diff(time))
    vr = np.diff(rs) / np.diff(time) * RSUN_KM

    # Summary plots
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))
    fig.suptitle(f'CME Analysis from {output_prefix}', fontsize=14)

    axs[0, 0].plot(time, freq, 'ro-')
    axs[0, 0].set_xlabel("Time [s]"); axs[0, 0].set_ylabel("Freq [MHz]"); axs[0, 0].grid()

    axs[0, 1].plot(time[:-1], dfdt, 'go-')
    axs[0, 1].set_xlabel("Time [s]"); axs[0, 1].set_ylabel("df/dt [MHz/s]"); axs[0, 1].grid()

    axs[0, 2].plot(time, rs, 'bo-')
    axs[0, 2].set_xlabel("Time [s]"); axs[0, 2].set_ylabel("Height [Rsun]"); axs[0, 2].grid()

    axs[1, 0].plot(rs, freq, 'mo-')
    axs[1, 0].set_xlabel("Height [Rsun]"); axs[1, 0].set_ylabel("Freq [MHz]"); axs[1, 0].grid()

    axs[1, 1].plot(rs[:-1], dfdt, 'co-')
    axs[1, 1].set_xlabel("Height [Rsun]"); axs[1, 1].set_ylabel("Drift [MHz/s]"); axs[1, 1].grid()

    axs[1, 2].plot(time[:-1], vr, 'ko-')
    axs[1, 2].set_xlabel("Time [s]"); axs[1, 2].set_ylabel("Speed [km/s]"); axs[1, 2].grid()

    fig.tight_layout()
    plt.savefig(output_prefix + '.results.png')

    # Save results
    with open(output_prefix + '.table.txt', 'w') as f:
        f.write("T[s], F[MHz], Ne[cm^-3], Rs[Rsun]\n")
        for i in range(len(freq)):
            f.write(f"{time[i]:.2f}, {freq[i]:.2f}, {Ne[i]:.1f}, {rs[i]:.2f}\n")

    print("\nSpeed Statistics:")
    print(f"Mean speed   : {np.mean(vr):.1f} km/s")
    print(f"Median speed : {np.median(vr):.1f} km/s")
    v_first_order = (rs[-1] - rs[0]) / (time[-1] - time[0]) * RSUN_KM
    print(f"1st-order est: {v_first_order:.1f} km/s")
