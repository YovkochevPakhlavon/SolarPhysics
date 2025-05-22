# -*- coding: utf-8 -*-
"""
Background Subtraction by a Specific Spectrum

This script loads a FITS file from the e-Callisto network,
subtracts the radio spectrum at a specific time from the entire dataset,
and plots the result using a selected colormap.

Useful link: http://soleil.i4ds.ch/solarradio/callistoQuicklooks/
Find your file and download it into the same folder or specify the path.

@author: Christian Monstein
"""

# ------------------------------------------------------------------------------

import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

# ------------------------------------------------------------------------------

# Define your local directory and filename (replace with your own file path)
mypath = r'C:\Users\VivoBook\OneDrive\Desktop\python'  # Use raw string (r'') to avoid escape issues
myfile = 'GAURI_20110809_075959_59.fit.gz'
# Combine into full path
full_path = os.path.join(mypath, myfile) 

# Open FITS file
hdu = fits.open(full_path)
data = hdu[0].data.astype(np.float32)          # Spectral data
freqs = hdu[1].data['Frequency'][0]            # Frequency array in MHz
time = hdu[1].data['Time'][0]                  # Time array in seconds
hdu.close()

# ------------------------------------------------------------------------------

# Define the background subtraction time (in seconds from the start)
T = 600  # seconds
timemarker = 4 * T  # 4 samples per second, so 4*T gives the index at time T

# Ensure index is within range
if timemarker >= data.shape[1]:
    raise ValueError(f"Time marker index {timemarker} exceeds data range: {data.shape[1]}")

# Subtract the spectrum at time = T
reference_spectrum = data[:, timemarker]  # A single spectrum at time T
data_subtracted = data - reference_spectrum[:, np.newaxis] + 10  # Subtract, keep values positive

# Set up plot range
extent = (time[0], time[-1], freqs[-1], freqs[0])

# Plot the result
plt.figure(figsize=(12, 8))
plt.imshow(data_subtracted, aspect='auto', extent=extent, cmap=cm.plasma, vmin=10, vmax=40)
plt.plot(T, 400, marker=11, color="black", label=f'Reference @ {T}s')  # Mark reference point
plt.tick_params(labelsize=14)
plt.xlabel(f'Time [s] of FIT file: {myfile}', fontsize=15)
plt.ylabel('Observed Frequency [MHz]', fontsize=15)
plt.title(f'Background Subtracted at T = {T} seconds', fontsize=15)
plt.legend()
plt.colorbar(label="Intensity [relative units]")
#plt.savefig(myfile.replace('.fit.gz', '_subtracted.png'))
plt.show()

# ------------------------------------------------------------------------------
