
# Project: Nitrogen Spectrum Analysis

This project aims to analyze nitrogen spectra to evaluate the resolving power (RP) of a given nitrogen spectrum. The primary functionality is provided by the `fit_n2` function, which attempts to automatically fit a nitrogen spectrum and calculate the resolving power.

## Files in the Project

### Data Creation and Analysis Scripts

- **create_data.py**: Script for generating synthetic data for testing and calibration.
- **nitrogen_fit.py**: Contains the `fit_n2` function and related fitting routines for analyzing nitrogen spectra.
- **nitrogen_plot.py**: Utilities for plotting the results of the nitrogen spectrum analysis.
- **PGM_calibration.py**: Script for calibrating the Plane Grating Monochromator (PGM) using the nitrogen spectrum.

### Example Scripts

- **example_energy_back_and_forth.py**: Example script to verify the consistency of energy measurements.
- **example_n2_fit.py**: Example script demonstrating the use of the `fit_n2` function in `nitrogen_fit.py`.
- **example_PGMCalibration.py**: Example script for the PGM calibration routines in `PGM_calibration.py`, using data extracted from the `EneryCalibrationPaper.pdf`.
- **example_PGMCalibration_jens_data.py**: Additional example for PGM calibration using specific data sets.

### Configuration and Documentation

- **requirements.txt**: List of Python dependencies required for the project.
- **CalibrationProcedure.md**: Explanation of the calibration procedure as I understood it from the paper.
- **EnergyCalibrationPaper.pdf**: the energy calibration paper

## How to Use the Project

1. **Install Dependencies**

   Before running the scripts, install the necessary dependencies listed in `requirements.txt` using pip:

   ```bash
   pip install -r requirements.txt
   ```

2. **Generate Data**

   Synthetic data obtained by artificially shifting the orignal N2 spectra recorded at PO4 was created using the `create_data.py` script. One can modify the script to create different dataset if necessary, and then execute it with:

   ```bash
   python create_data.py
   ```


## Examples

The project includes several example scripts to demonstrate the usage of the program:

- **example_energy_back_and_forth.py**: Verifies the consistency of energy measurements. At the moment I have some discrepancy, probably due to some approximation.
- **example_n2_fit.py**: Fit the synthetic data obtained shifting PO4 data. Individual plots are saved in `Results/N2_cff*.png` and the fit results are saved in `Results/N2_fit_results.csv`
- **example_PGMCalibration.py**: Demonstrates the PGM calibration routines using the data from the paper.
- **example_PGMCalibration_PO4_data.py**: Additional example for PGM calibration using a the dataset executing the `example_n2_fit.py` script.

Run the examples using Python:

```bash
python example_energy_back_and_forth.py
python example_n2_fit.py
python example_PGMCalibration.py
python example_PGMCalibration_PO4_data.py
```
