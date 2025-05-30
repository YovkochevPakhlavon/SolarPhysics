# Callisto FITS Data Analysis with Python

This repository provides step-by-step Python scripts and Jupyter Notebooks for analyzing solar radio spectrograms in `.fit.gz` format (from instruments like CALLISTO). The aim is to help learners learn how to:

- Load and visualize FITS files
- Perform background subtraction
- Merge data from different instruments
- Improve axis labeling and visual clarity
- Download and handle data from online sources

---

##  Folder Structure

### `Basic_Image_Processing/`
**Description:**  
The simplest working script to load and visualize a `.fit.gz` file.

**Files:**
- `mostbasic.py` — Load and plot a FITS file using Matplotlib.
- `webaccess.py` — Demonstrates how to download CALLISTO data directly from the web and visualize it. Uses `urllib` or similar to fetch and display data.
- `fit_preprocess_from_url.py` -Download a FITS file from a given URL and perform different background subtraction methods on the spectrogram data. Saves each processed spectrogram as PNG images.

---

### `BackgroundSubtraction/`
**Description:**  
Focuses on removing background noise from the data for clearer signal detection.

**Files:**
- `constant.py` — Subtracts the mean signal (a constant background) from the data.
- `individual.py` — Subtracts a user-defined background at a specific time or frequency range.
- `different_background_subtraction.py` - Module for processing FITS radio data files (Using different background subtraction techniques) and visualizing spectrograms with various enhancements.

---

### `MergeFrequency/`
**Description:**  
Combines FITS data from different instruments (e.g., BIR and ALMATY) in the frequency domain.

**Files:**
- `MergeFrequency.py` — Merges and displays spectrograms from different FITS files in one unified plot.
- `merge_fits_spectrogram.py` -This script provides tools to merge two spectrogram FITS files by their frequency axes and visualize the result. It supports optional user-defined frequency slicing and handles automatic frequency sorting.

---

### `Axis/`
**Description:**  
Improves time and frequency axis labeling for better readability.

**Files:**
- `spectrogram_visualization.ipynb` — A Jupyter Notebook that:
  - Converts time axis to human-readable UT format
  - Applies Gaussian filtering for visual smoothing
  - Uses advanced colormaps for clarity

---

### `SolarFlareAnalyze/`
**Description:**  
Analyzes Solar Flares.

**Files:**
- `solar_flare_gauss_fit.py` — This module downloads and processes solar flare flux data from a GOES FITS file, performs Gaussian fitting on a selected time window of the flux, and plots the results.

- `goes_xray_flux_analysis.py` — This module loads GOES-16 X-ray flux data from a netCDF file, extracts specified variables, identifies peak and surrounding minima, and optionally plots the flux time series with annotated peaks and minima.

---

### `CME_Analysis/`
**Description:**  
Callisto CME Analyzer

This scripts read one or more FITS files containing radio dynamic spectra (from Callisto instruments), plot the data, and allow the user to select a frequency drift (typically of a CME). The scripts then calculate:

- Electron density
- Solar radial distance 
- Frequency drift rate (df/dt)
- Velocity of the emission front

**Files:**
- `Newkirk_model.py` — Uses Newkirk model to compute the electron density, frequency drift rate (df/dt), height and CME speed.

- `Leblanc_model.py` — Uses Leblanc model to compute the electron density, frequency drift rate (df/dt), height and CME speed.

---
##  Requirements

You can install the required packages using pip:

```bash
pip install numpy matplotlib astropy scipy
````

Optional (for notebooks):

```bash
pip install jupyter
```

---

## 🚀 How to Run

1. Clone this repository:

```bash
git clone https://github.com/YovkochevPakhlavon/SolarPhysics.git
cd SolarPhysics
```

2. Explore folders and scripts individually:

   * For scripts (`.py`): run them using a Python IDE or the terminal.
   * For Jupyter Notebooks: start Jupyter with:

   ```bash
   jupyter notebook
   ```

3. Download sample FITS data from:
   [CALLISTO Data Quicklook](http://soleil.i4ds.ch/solarradio/callistoQuicklooks/)

---



## Credits

Scripts adapted and extended from the educational material by **Christian Monstein** and other contributors.
Converted and improved for teaching and analysis by Yovkochev Pakhlavon, Joseph Alowonle, Muhammet Ozcan.

---

## 📬 Questions?

If you encounter any issues or have questions, feel free to open an issue or contact the maintainer.

```
