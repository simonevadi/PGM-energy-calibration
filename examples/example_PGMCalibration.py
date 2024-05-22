import os
import numpy as np
from tabulate import tabulate

from PGM_calibration import PGMCalibration

N = 1200  # Line density (grating grooves per mm)
k = 1  # Order of diffraction

os.makedirs('Results/ExampleCalibration', exist_ok=True)

# Example data
cff_values = np.array([       1.6,   1.8,   2.0,   2.25,  3.5,   5.0,    8.0,   10.0,  15.0])
measured_energies = np.array([398.3, 398.1, 397.9, 397.7, 397.1, 396.76, 396.4, 396.3, 396.17])
orders = np.full(cff_values.shape[0], 1)

# Instantiate the calibration class
c = PGMCalibration(N)

# Set a reasonable initial guess manually
c.set_initial_guess(DTheta=0, DBeta=0, E=401)  # Example initial guess values

# Perform fitting
dtheta, dbeta, E_opt = c.fit_parameters(measured_energies, cff_values, orders, delta_E=3)

c.print_fit_results()
# Plot the fitting results along with initial guess
c.plot_fit(measured_energies, cff_values, orders,
           DTheta=dtheta, DBeta=dbeta, E_opt=E_opt,
           show=False, savepath='Results/ExampleCalibration/fit_results_paper_data.png')


rolf_dtheta = 0.034562
rolf_dbeta = 0.036857
en_rolf = 401.766
# Convert dtheta and dbeta from radians to degrees
dtheta_deg = np.rad2deg(dtheta)
dbeta_deg = np.rad2deg(dbeta)

# Rolf's results
rolf_dtheta = 0.034562
rolf_dbeta = 0.036857
en_rolf = 401.766

# Calculate the percentage difference
diff_dtheta = 100 * (dtheta_deg - rolf_dtheta) / rolf_dtheta
diff_dbeta = 100 * (dbeta_deg - rolf_dbeta) / rolf_dbeta
diff_e_opt = 100 * (E_opt - en_rolf) / en_rolf

# Create the table
headers = ["Parameter", "Fit Results", "Rolf's Results", "Difference"]
data = [
    ["DTheta [deg]", f"{dtheta_deg:.6f}", f"{rolf_dtheta:.6f}", f"{diff_dtheta:.2f}%"],
    ["DBeta  [deg]", f"{dbeta_deg:.6f}", f"{rolf_dbeta:.6f}", f"{diff_dbeta:.2f}%"],
    ["E_opt  [eV]", f"{E_opt:.6f}", f"{en_rolf:.6f}", f"{diff_e_opt:.2f}%"]
]

# Print the table
print(tabulate(data, headers, tablefmt="pretty", colalign=("left", "right", "right", "right")))