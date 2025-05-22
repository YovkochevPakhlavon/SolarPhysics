# -*- coding: utf-8 -*-
"""
Background Subtraction by a Constant Value
This script loads a FITS file from a local directory, applies a constant background subtraction,
and visualizes the data using a selected colormap.

Find your file here: http://soleil.i4ds.ch/solarradio/callistoQuicklooks/
Download it to your working folder before running this script.

@author: Christian Monstein
Modified for educational clarity by: Yovkochev Pakhlavon
"""

# ------------------------------------------------------------------------------

import os
from astropy.io import fits                # For reading FITS files
import matplotlib.pyplot as plt            # For plotting
import numpy as np                         # For data processing
from matplotlib import cm                  # For colormap options



# ------------------------------------------------------------------------------

# Define your local directory and filename (replace with your own file path)
mypath = r"C:\Users\VivoBook\OneDrive\Desktop\python"  # Use raw string (r'') to avoid escape issues
myfile = 'GAURI_20110809_075959_59.fit.gz'
# Combine into full path
full_path = os.path.join(mypath, myfile) 

# Open FITS file and extract data
hdu = fits.open(full_path)
data = hdu[0].data.astype(np.float32)
freqs = hdu[1].data['Frequency'][0]  # Frequency axis in MHz
time = hdu[1].data['Time'][0]        # Time axis in seconds
hdu.close()

# ------------------------------------------------------------------------------

# Apply background subtraction (mean + offset)
myoffset = 10  # Background offset to subtract
data_subtracted = data - np.mean(data) - myoffset

# Set the plotting extent (so axes are labeled correctly)
extent = (time[0], time[-1], freqs[-1], freqs[0])

# Create the plot
plt.figure(figsize=(12, 8))
plt.imshow(data_subtracted, aspect='auto', extent=extent, cmap=cm.inferno, vmin=-15, vmax=30)

# Customize plot
plt.tick_params(labelsize=14)
plt.xlabel('Time [s] of FIT file: ' + myfile, fontsize=15)
plt.ylabel('Observed Frequency [MHz]', fontsize=15)
plt.title(f'Background Subtracted: Offset = {myoffset}', fontsize=15)
plt.colorbar(label='Relative Intensity')

# Save the figure
#plt.savefig(myfile + "_subtracted.png")
plt.show()

# ------------------------------------------------------------------------------
