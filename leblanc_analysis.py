"""
LEBLANC CME SPEED ANALYSIS TOOL
--------------------------------

Created on Tue Aug 27 15:34:53 2024


Description:
This script reads Callisto FITS dynamic spectra, allows user interaction
to select CME fronts, and uses the Leblanc density model to compute electron
density, height, and CME speed.

Leblanc Model:
    n(r) = 1.36e6 * r^-2.14 + 1.68e8 * r^-6.13
    Where:
        r: radial distance [Rsun]
        n: electron density [cm^-3]
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from scipy.optimize import fsolve
import os

# Global constants
RSUN_KM = 695700.0  # Solar radius in km
LEBLANC_COEFFS = {'a': 1.36e6, 'b': 1.68e8, 'alpha': 2.14, 'beta': 6.13}
ELECTRON_CONSTANT = 8.977e-3  # MHz * sqrt(cm^-3)

def leblanc_model_equation(r, n_value):
    """Leblanc density model equation to solve for radius."""
    a = LEBLANC_COEFFS['a']
    b = LEBLANC_COEFFS['b']
    alpha = LEBLANC_COEFFS['alpha']
    beta = LEBLANC_COEFFS['beta']
    return a * r**-alpha + b * r**-beta - n_value

def find_r_from_density(n):
    """
    Solves the Leblanc model for the radial distance given electron density.
    
    Parameters
    ----------
    n : float
        Electron density [cm^-3]
        
    Returns
    -------
    float
        Height [Rsun]
    """
    initial_guess = 1.0
    r_solution = fsolve(leblanc_model_equation, initial_guess, args=(n,))
    return r_solution[0]

def run_leblanc_analysis(fits_path, file_name, harmonic=1, zoom=[1200, 1800, 20, 180], vmin=-5, vmax=25):
    """
    Runs the Leblanc-based CME speed analysis on a given FITS file.

    Parameters
    ----------
    fits_path : str
        Directory path to the FITS file.
    file_name : str
        FITS file name.
    harmonic : int, optional
        Excitation mode (1 = fundamental, 2 = harmonic).
    zoom : list of float, optional
        Zoom view [start_time, end_time, min_freq, max_freq].
    vmin : float, optional
        Minimum dB for plot color scale.
    vmax : float, optional
        Maximum dB for plot color scale.
    """
    full_path = os.path.join(fits_path, file_name)
    hdu = fits.open(full_path)
    header = hdu[0].header
    image = hdu[0].data
    T0 = header['TIME-OBS']
    title = header['CONTENT']
    dt = header['CDELT1']
    freq = hdu[1].data['Frequency'][0]
    time_axis = np.arange(image.shape[1]) * dt
    hh, mm, ss = map(float, T0.split(":"))
    t0_sec = hh * 3600 + mm * 60 + ss
    time_sec = t0_sec + time_axis
    extent = (time_axis[0], time_axis[-1], freq[-1], freq[0])
    
    dflat = gaussian_filter(image - np.mean(image, axis=1, keepdims=True), sigma=1)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(dflat, aspect='auto', extent=extent, cmap='jet',
                   norm=plt.Normalize(vmin, vmax))
    ax.set_xlabel(f"Time [s] after {T0} UT")
    ax.set_ylabel("Frequency [MHz]")
    ax.set_title(title)
    ax.axis(zoom)
    fig.colorbar(im, ax=ax, label='dB')

    coords = []

    class MouseMonitor:
        def __init__(self, fig, ax):
            self.fig = fig
            self.ax = ax

        def __call__(self, event):
            if event.button == 3:  # Right-click ends session
                fig.canvas.mpl_disconnect(cid)
                analyze_clicks(coords, harmonic, freq, os.path.splitext(full_path)[0])
            else:
                x, y = event.xdata, event.ydata
                self.ax.plot(x, y, 'wo')
                self.fig.canvas.draw_idle()
                coords.append((x, y))

    mouse = MouseMonitor(fig, ax)
    cid = fig.canvas.mpl_connect('button_press_event', mouse)
    plt.show()

def analyze_clicks(coords, harmonic, freq_axis, out_prefix):
    """
    Process user-selected points and compute plasma parameters and speeds.

    Parameters
    ----------
    coords : list of tuples
        Clicked (time, frequency) points.
    harmonic : int
        Fundamental or harmonic mode.
    freq_axis : np.ndarray
        Frequency array for scale.
    out_prefix : str
        File path prefix for outputs.
    """
    times, freqs = zip(*coords)
    times = np.array(times)
    freqs = np.array(freqs)

    Ne = (freqs / (harmonic * ELECTRON_CONSTANT)) ** 2
    Rs = np.array([find_r_from_density(n) for n in Ne])
    dfdt = np.diff(freqs) / np.diff(times)
    speeds = np.diff(Rs) / np.diff(times) * RSUN_KM

    # Plotting results
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))
    fig.suptitle(f'CME Speed Analysis using Leblanc Model: {out_prefix}', fontsize=14)

    axs[0, 0].plot(times, freqs, 'ro-')
    axs[0, 0].set_title('Selected Frequencies')
    axs[0, 0].set_xlabel('Time [s]'); axs[0, 0].set_ylabel('Freq [MHz]'); axs[0, 0].grid()

    axs[0, 1].plot(times[:-1], dfdt, 'go-')
    axs[0, 1].set_title('Drift Rate (df/dt)')
    axs[0, 1].set_xlabel('Time [s]'); axs[0, 1].set_ylabel('MHz/s'); axs[0, 1].grid()

    axs[0, 2].plot(times, Rs, 'bo-')
    axs[0, 2].set_title('Heights from Leblanc')
    axs[0, 2].set_xlabel('Time [s]'); axs[0, 2].set_ylabel('Height [Rsun]'); axs[0, 2].grid()

    axs[1, 0].plot(Rs, freqs, 'mo-')
    axs[1, 0].set_title('Freq vs Height')
    axs[1, 0].set_xlabel('Height [Rsun]'); axs[1, 0].set_ylabel('Freq [MHz]'); axs[1, 0].grid()

    axs[1, 1].plot(Rs[:-1], dfdt, 'co-')
    axs[1, 1].set_title('Drift vs Height')
    axs[1, 1].set_xlabel('Height [Rsun]'); axs[1, 1].set_ylabel('df/dt [MHz/s]'); axs[1, 1].grid()

    axs[1, 2].plot(times[:-1], speeds, 'ko-')
    axs[1, 2].set_title('Speed')
    axs[1, 2].set_xlabel('Time [s]'); axs[1, 2].set_ylabel('Speed [km/s]'); axs[1, 2].grid()

    fig.tight_layout()
    plt.savefig(out_prefix + '_results.png')

    # Save table
    with open(out_prefix + '_table.txt', 'w') as f:
        f.write('Time[s], Freq[MHz], Ne[cm^-3], Rsun[R]\n')
        for t, f_, n, r in zip(times, freqs, Ne, Rs):
            f.write(f"{t:.2f}, {f_:.2f}, {n:.2e}, {r:.2f}\n")

    # Print summary
    print('\nStatistical CME Velocity:')
    print(f"Mean Speed     = {np.mean(speeds):.1f} km/s")
    print(f"Median Speed   = {np.median(speeds):.1f} km/s")
    v1 = (Rs[-1] - Rs[0]) / (times[-1] - times[0]) * RSUN_KM
    print(f"1st Order Speed = {v1:.1f} km/s")


# from leblanc_analysis import run_leblanc_analysis

# run_leblanc_analysis(
#     fits_path="C:/Users/Analysis/",
#     file_name="LEARMONTH_20230612_064000_01.fit.gz",
#     harmonic=1,
#     zoom=[1200, 1800, 20, 180],
#     vmin=-5,
#     vmax=25
# )
