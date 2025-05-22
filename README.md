````markdown
## Callisto FITS Data Analysis with Python

This repository provides step-by-step Python scripts and Jupyter Notebooks for analyzing solar radio spectrograms in `.fit.gz` format (from instruments like CALLISTO). The aim is to help learners learn how to:

- Load and visualize FITS files
- Perform background subtraction
- Merge data from different instruments
- Improve axis labeling and visual clarity
- Download and handle data from online sources

---

##  Folder Structure

### `basic/`
**Description:**  
The simplest working script to load and visualize a `.fit.gz` file.

**Files:**
- `mostbasic.py` — Load and plot a FITS file using Matplotlib.
- `webaccess.py` — Demonstrates how to download CALLISTO data directly from the web and visualize it. Uses `urllib` or similar to fetch and display data.

---

### `BackgroundSubtraction/`
**Description:**  
Focuses on removing background noise from the data for clearer signal detection.

**Files:**
- `constant.py` — Subtracts the mean signal (a constant background) from the data.
- `individual.py` — Subtracts a user-defined background at a specific time or frequency range.

---

### `MergeFrequency/`
**Description:**  
Combines FITS data from different instruments (e.g., BIR and ALMATY) in the frequency domain.

**Files:**
- `MergeFrequency.py` — Merges and displays spectrograms from different FITS files in one unified plot.

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
