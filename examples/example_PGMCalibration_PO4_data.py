import numpy as np
import pandas as pd
import os

from PGM_calibration import PGMCalibration

N = 1200  # Line density (grating grooves per mm)

os.makedirs('Results/ExampleCalibrationPO4', exist_ok=True)
# Example data
df = pd.read_csv('Results/ExampleN2Fit/N2_fit_results.csv')

cff_values = np.array(df['cff'])
measured_energies = np.array(df['peak_energy'])
orders = np.full(cff_values.shape[0], 1)

# Instantiate the calibration class
c = PGMCalibration(N)

# Set a reasonable initial guess manually
c.set_initial_guess(DTheta=0, DBeta=0, E=398)  # Example initial guess values

# Perform fitting
dtheta, dbeta, E_opt = c.fit_parameters(measured_energies, cff_values, orders)

c.print_fit_results()
# Plot the fitting results along with initial guess
c.plot_fit(measured_energies, cff_values, orders,
           DTheta=dtheta, DBeta=dbeta, E_opt=E_opt,
           show=False, savepath='Results/ExampleCalibrationPO4/fit_results_PO4_data.png')


