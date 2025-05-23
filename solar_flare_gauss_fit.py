"""
solar_flare_gauss_fit.py

Created on Wed Feb 07 19:23:25 2018

Author: Adapted from Christian Monstein’s code

Adapted and modularized by: Abraham-Alowonle Joseph-judah

Description:
------------
This module downloads and processes solar flare flux data from a GOES FITS file,
performs Gaussian fitting on a selected time window of the flux, and plots
the results.

Requirements:
-------------
- numpy
- astropy
- matplotlib
- lmfit (install via `pip install lmfit`)

Usage:
------
In Jupyter or Spyder, import the main function and call it with your FITS filename:
>>> from solar_flare_gauss_fit import process_and_fit_flux
>>> process_and_fit_flux('go1520180122.fits')

Or run this script standalone with:
>>> python solar_flare_gauss_fit.py go1520180122.fits

"""

import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from lmfit import Model


def gaussian(x, amp, cen, wid):
    """
    Defines a 1D Gaussian function.

    Parameters
    ----------
    x : array_like
        Independent variable.
    amp : float
        Amplitude of the Gaussian.
    cen : float
        Center (mean) position of the Gaussian.
    wid : float
        Standard deviation (width) of the Gaussian.

    Returns
    -------
    array_like
        Gaussian function evaluated at x.
    """
    return (amp / (np.sqrt(2 * np.pi) * wid)) * np.exp(-(x - cen) ** 2 / (2 * wid ** 2))


def process_and_fit_flux(fits_filename, time_window=(11000, 16000), plot_limits=None):
    """
    Load GOES solar flare FITS file, plot the flux, select a time window,
    fit a Gaussian to the flux in that window, and plot the fit.

    Parameters
    ----------
    fits_filename : str
        Path to the GOES FITS file.
    time_window : tuple, optional
        Start and end time in seconds to select data for Gaussian fitting (default (11000, 16000)).
    plot_limits : dict or None, optional
        Dictionary with keys 'xlim' and 'ylim' for plot axis limits, e.g.,
        {'xlim': (5000, 20000), 'ylim': (0, 1e-6)}. Defaults to None.

    Returns
    -------
    result : lmfit.model.ModelResult
        The fit result object containing fit parameters and statistics.
    """
    # Open FITS file and extract data
    with fits.open(fits_filename) as hdu:
        flux1 = hdu[2].data[0][1][:, 0]
        flux2 = hdu[2].data[0][1][:, 1]

        date = hdu[0].header.get('DATE-OBS', 'Unknown date')
        time_obs = hdu[0].header.get('TIME-OBS', 'Unknown time')
        unit = hdu[2].header.get('TUNIT2', 'unknown unit')
        telescope = hdu[2].header.get('TELESCOP', 'Unknown telescope')

    dT = 2.048  # exposure time in seconds
    taxis = np.arange(flux1.size) * dT

    # Plot original flux
    plt.figure(figsize=(10, 7))
    plt.plot(taxis, flux1, label='Flux1')
    plt.xlabel('Time [s]')
    plt.ylabel(f'Flux [{unit}]')
    if plot_limits and 'xlim' in plot_limits:
        plt.xlim(*plot_limits['xlim'])
    else:
        plt.xlim(5000, 20000)
    if plot_limits and 'ylim' in plot_limits:
        plt.ylim(*plot_limits['ylim'])
    else:
        plt.ylim(0, 1e-6)
    plt.title(f"{telescope} of {date} at {time_obs}")
    plt.grid(True)

    # Select time window for fitting
    mask = (taxis > time_window[0]) & (taxis < time_window[1])
    flux = flux1[mask]
    time = taxis[mask]
    plt.plot(time, flux, '.r', label='Fit window')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Gaussian model fitting
    gmodel = Model(gaussian)
    # Initial parameters guess: amplitude negative because flux dips? Adjust if necessary
    result = gmodel.fit(flux, x=time, amp=-1e-6, cen=time_window[0], wid=2000)

    print(result.fit_report())

    # Plot fit result
    plt.figure(figsize=(10, 7))
    plt.plot(time, flux, 'b', label='Data')
    plt.plot(time, result.best_fit, 'r', label='Gaussian fit')
    plt.xlabel('Time [s]')
    plt.ylabel(f'Flux [{unit}]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return result


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python solar_flare_gauss_fit.py <fits_filename>")
        sys.exit(1)

    fits_file = sys.argv[1]
    process_and_fit_flux(fits_file)

# from solar_flare_gauss_fit import process_and_fit_flux

# result = process_and_fit_flux('go1520180122.fits')
