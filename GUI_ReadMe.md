# PGM Energy Calibration GUI

This GUI application is designed for calibrating PGM (Plane Grating Monochromator) energy using measured energies and CFF values. The application allows users to load data, set grating values, perform calibration fitting, and save the results.

## Requirements

- Python 3.x
- Tkinter
- Matplotlib
- Pandas
- NumPy
- PGM_calibration module (custom module, ensure it is in the same directory or properly installed)


## How to Use
1. **Load Data File**

    * Click the Browse button *Browse* to select a CSV file containing the data. 
    * Alternatively you can manually input the file path and press *Load File* The selected file path will appear in the entry field.
    * A message will display indicating whether the file was successfully loaded.
2. **Set Grating**
    * Enter the grating value in the next to the button *Set Grating*.
    * Click the *Set Grating* button or press Enter to set the grating value.
    * A message will display indicating the grating value has been set.
3. **Perform Fit**
    * Click the *Fit* button to perform the calibration fitting using the loaded data and set grating value.
    * The fitted parameters (*DTheta, DBeta, E_opt*) will be displayed in the application.
4. **Set Beamline Name**
    * Enter the beamline name in the *Set Beamline Name entry field*.
    * Click the *Set Beamline Name* button or press Enter to set the beamline name. A message will display indicating the beamline name has been set.
5. **Save Results**
    * Click the *Browse* button to select a folder where the results will be saved.
    * Enter the desired file name in the entry field. The current folder path will be displayed in the label below the entry field.
    * Click the Save button to save the results, including the input data, fit results, and a plot as a PDF.
    * The following files will be saved
        * Input File: beamline_name_input_data.csv
        * Fit Results File: beamline_name_fit_results.csv
        * Plot: beamline_name_UE112.pdf

## Example Input File Format
The file must be in csv format, and contain the header as well as all the data. You can see an example in the gui folder, inside the folder *gui_data*
| orders | cff_values | measured_energies |
|--------|-------------|-------------------|
| 1      | 0.2         | 91.451            |
| 1      | 0.5         | 91.593            |
| 1      | 1.5         | 90.681            |
| 1      | 6.3         | 90.941            |
| 2      | 0.2         | 91.383            |
| 2      | 0.5         | 91.498            |
| 2      | 1.5         | 90.776            |
| 2      | 6.3         | 91.004            |
| 3      | 1.5         | 90.834            |
| 3      | 6.3         | 91.036            |

