import os
import numpy as np

from PGM_calibration import PGMCalibration

N = 1200  # Line density (grating grooves per mm)
# k = 1  # Order of diffraction

os.makedirs('Results/ExampleCalibration', exist_ok=True)

# Example data
# orders = np.array([           1,      1,     1,      2,      2,      2])
# cff_values = np.array([       1.5,    3,     8,      1.5,    3,      8,])
# measured_energies = np.array([92.819, 91.84, 91.416, 2*92.295, 2*91.641, 2*91.341 ])


orders = np.array([           1,      1,     1])
cff_values = np.array([       1.5,    3,     8])
measured_energies = np.array([92.819, 91.84, 91.416])

# orders = np.array([         2,      2,      2])
# cff_values = np.array([      1.5,    3,      8,])
# measured_energies = np.array([2*92.295, 2*91.641, 2*91.341 ])



# Instantiate the calibration class
c = PGMCalibration(N)

# Set a reasonable initial guess manually
c.set_initial_guess(DTheta=5e-4, DBeta=-1e-4, E=91.1)  # Example initial guess values

# Perform fitting
dtheta, dbeta, E_opt = c.fit_parameters(measured_energies, cff_values, orders)

c.print_fit_results()
# Plot the fitting results along with initial guess
c.plot_fit(measured_energies, cff_values, orders,
           DTheta=dtheta, DBeta=dbeta, E_opt=E_opt,
           show=False, savepath='Results/UE112_2016/fit_results_ue112_2016.png')


