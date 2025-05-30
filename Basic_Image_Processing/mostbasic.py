# -*- coding: utf-8 -*-
"""
The most basic Python script to plot a FITS file from e-CALLISTO

This script demonstrates how to:
- Open a FITS file using Astropy
- Extract image data from the FITS file
- Visualize the data using Matplotlib
- Save the resulting plot to an image file

Author: Christian Monstein
"""

# --------------------------------------------------------------------------------
# Import required libraries
# --------------------------------------------------------------------------------
import os
from astropy.io import fits      # For reading FITS files
import matplotlib.pyplot as plt  # For plotting the data

# --------------------------------------------------------------------------------
# Load the FITS file
# --------------------------------------------------------------------------------

# Open the compressed FITS file
 
# Define your local directory and filename (replace with your own file path)
mypath = r'C:\Users\VivoBook\OneDrive\Desktop\python'  # Use raw string (r'') to avoid escape issues
myfile = 'MONGOLIA-UB_20240820_000000_01.fit.gz'
# Combine into full path
full_path = os.path.join(mypath, myfile) 

fits_file   = fits.open(full_path)    # Open the FITS file
data = fits_file[0].data             # Access the image data in the primary HDU
fits_file.close()                    # Close the file to free up resources

# --------------------------------------------------------------------------------
# Plot the data
# --------------------------------------------------------------------------------

plt.imshow(data, aspect='auto')      # Display the data as an image
plt.colorbar(label='Intensity')      # Add a colorbar to show intensity scale
plt.title('e-CALLISTO Radio Spectrogram')  # Add a plot title
plt.xlabel('Time [s]')         
plt.ylabel('Frequency [MHz]')    
#plt.savefig('mostbasic.png')         # Save the plot as a PNG image
plt.show()                           # Show the plot on screen
