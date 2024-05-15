import numpy as np

from PGM_calibration import PGMCalibration

N = 1200  # Line density (grating grooves per mm)
k = 1  # Order of diffraction

# Example data
cff_values = np.array([1.6, 1.8, 2.0, 2.25, 3.5, 5.0, 8.0, 10.0, 15.0])
measured_energies = np.array([398.3, 398.1, 397.9, 397.7, 397.1, 396.76, 396.4, 396.3, 396.17])

# Instantiate the calibration class
c = PGMCalibration(N, k)

# Set a reasonable initial guess manually
c.set_initial_guess(DTheta=0, DBeta=0, E=398)  # Example initial guess values

# Perform fitting
dtheta, dbeta, E_opt = c.fit_parameters(measured_energies, cff_values)

c.print_fit_results()
# Plot the fitting results along with initial guess
c.plot_fit(measured_energies, cff_values,
           DTheta=dtheta, DBeta=dbeta, E_opt=E_opt,
           show=False, savepath='fit_results.png')


