# -*- coding: utf-8 -*-
"""
Downloading, Reading, and Plotting e-CALLISTO FITS Data

This script demonstrates:
- How to download a FITS file from the e-CALLISTO archive
- How to extract frequency/time metadata and image data
- How to plot and save a radio spectrogram

Author: Christian Monstein

Useful link to find FITS files:
http://soleil.i4ds.ch/solarradio/callistoQuicklooks/
"""

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

import urllib.request                      # For downloading data from a URL
from astropy.io import fits                # For reading FITS files
import matplotlib.pyplot as plt            # For plotting
import numpy as np                         # For data processing
from matplotlib import cm                  # For colormap options

# ------------------------------------------------------------------------------
# Step 1 – Define the URL of the desired FITS file
# ------------------------------------------------------------------------------

# This example uses a file from the MONGOLIA-UB station on 20 August 2024
url = 'http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/2024/08/20/MONGOLIA-UB_20240820_000000_01.fit.gz'

# Extract the filename from the URL
filename = url.rsplit('/', 1)[1]

print("Downloading FITS file from:", url)
urllib.request.urlretrieve(url, filename)

# ------------------------------------------------------------------------------
# Step 2 – Read the FITS data
# ------------------------------------------------------------------------------

# Open the FITS file
fits_file = fits.open(filename)

# Extract the main image data from the primary HDU
# Convert to float32 to ensure consistent plotting
data = fits_file[0].data.astype(np.float32)

# Extract frequency and time metadata from the second HDU
freqs = fits_file[1].data['Frequency'][0]   # Frequency axis [MHz]
time = fits_file[1].data['Time'][0]         # Time axis [s]

# Close the FITS file
fits_file.close()

# ------------------------------------------------------------------------------
# Step 3 – Plot the spectrogram
# ------------------------------------------------------------------------------

plt.figure(figsize=(12, 8))

# Set axes range based on time and frequency arrays
extent = (time[0], time[-1], freqs[-1], freqs[0])  # flip frequency axis for spectrogram look

# Plot the image using imshow
plt.imshow(data, aspect='auto', extent=extent, cmap=cm.jet, vmin=110, vmax=170)

# Optional colormaps you can try: cm.inferno, cm.plasma, cm.cool, cm.magma, etc.

# Add plot labels and formatting
plt.colorbar(label='Intensity')
plt.tick_params(labelsize=14)
plt.xlabel('Time [s]', fontsize=15)
plt.ylabel('Frequency [MHz]', fontsize=15)
plt.title('Radio Spectrogram from MONGOLIA-UB – ' + filename, fontsize=16)

# Save the figure
#plt.savefig(filename + ".png")
plt.show()

# ------------------------------------------------------------------------------
