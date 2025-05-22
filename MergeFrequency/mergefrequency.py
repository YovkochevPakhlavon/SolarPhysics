# -*- coding: utf-8 -*-
"""
Merging Two Spectra in Frequency Space

This script merges FITS data from two different CALLISTO observatories
into one combined frequency spectrogram, adjusting for overlapping times
but different frequency bands.

Source: http://soleil.i4ds.ch/solarradio/callistoQuicklooks/
Download and save the files into your working directory.

@author: Christian Monstein
"""

# ------------------------------------------------------------------------------
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

# ------------------------------------------------------------------------------
mypath = r'C:\Users\VivoBook\OneDrive\Desktop\python'  # Use raw string (r'') to avoid escape issues
# === Load first FITS file (BIR - Ireland) ===
myfile1 = 'BIR_20110809_080000_59.fit.gz'
# Combine into full path
full_path1 = os.path.join(mypath, myfile1) 
hdu1 = fits.open(full_path1)
data1 = hdu1[0].data.astype(np.float32)
freqs1 = hdu1[1].data['Frequency'][0]
time1 = hdu1[1].data['Time'][0]
hdu1.close()

# === Load second FITS file (ALMATY - Kazakhstan) ===
myfile2 = 'ALMATY_20110809_080000_59.fit.gz'
# Combine into full path
full_path2 = os.path.join(mypath, myfile2) 
hdu2 = fits.open(full_path2)
data2 = hdu2[0].data.astype(np.float32)
freqs2 = hdu2[1].data['Frequency'][0]
time2 = hdu2[1].data['Time'][0]
hdu2.close()

# ------------------------------------------------------------------------------

# === Merge data in frequency space ===
# Select relevant frequency channels manually to align frequency gaps
data_combined = np.concatenate((data2[17:182, :], data1[10:178, :]), axis=0)
freqs_combined = np.hstack((freqs2[17:182], freqs1[10:178]))
time = time1  # Use time from BIR, assuming time alignment

# === Background subtraction (row-wise mean) ===
data_normalized = data_combined - data_combined.mean(axis=1, keepdims=True)

# ------------------------------------------------------------------------------

# === Plot the combined spectrogram ===
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)

extent = (time[0], time[-1], freqs_combined[-1], freqs_combined[0])

im = ax.imshow(data_normalized, aspect='auto', extent=extent,
               cmap=cm.plasma, vmin=-1, vmax=40)

# Optional: custom y-ticks if combined frequencies are not linearly spaced
# Comment out the next line if you want actual frequency labels
# ax.set_yticklabels([20, 40, 60, 80, 100, 200, 400, 600, 800])

plt.colorbar(im, label="Intensity [relative units]")
plt.tick_params(labelsize=14)
plt.xlabel('Time [s]', fontsize=15)
plt.ylabel('Observed Frequency [MHz]', fontsize=15)
plt.title('Merged Spectrogram: BIR (Ireland) + ALMATY (Kazakhstan)', fontsize=15)

# Save the plot
output_name = myfile1.replace('.fit.gz', '') + '_' + myfile2.replace('.fit.gz', '') + '_merged.png'
#plt.savefig(output_name)
plt.show()

# ------------------------------------------------------------------------------
