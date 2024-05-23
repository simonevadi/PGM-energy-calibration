import os
import numpy as np
from tabulate import tabulate
from PGM_calibration import PGMCalibration

N = 600  # Line density (grating grooves per mm)
# k = 1  # Order of diffraction

os.makedirs('Results/ExampleCalibration', exist_ok=True)

# Example data
orders = np.array([           1,      1,       1,        1,   2,      2,       2,        2,    3,    3   ])
cff_values = np.array([       0.2,    0.5,     1.5,      6.3, 0.2,    0.5,     1.5,      6.3,  1.5,   6.3  ])
measured_energies = np.array([91.451, 91.593, 90.681, 90.941, 91.383, 91.498, 90.776, 91.004,  90.834, 91.036 ])

# c = PGMCalibration(N)

# # Check if a value in cff_values is strictly lower than one
# mask = cff_values < 1
# # If the condition is met, make the corresponding value in orders negative
# orders[mask] *= -1
# print(orders)
# en = c.shifted_energy(cff_values, orders, 0, 0, measured_energies)
# for ind, e in enumerate(en):
#     print(f'original {measured_energies[ind]}, calculated {e}')

# Select only orders 3 and 4
selected_orders = [1,2, 3]
mask = np.isin(orders, selected_orders)

# Filtered data
orders = orders[mask]
cff_values = cff_values[mask]
measured_energies = measured_energies[mask]




# Instantiate the calibration class
c = PGMCalibration(N)

# Set a reasonable initial guess manually
c.set_initial_guess(DTheta=0.0, DBeta=0.0, E=91.84, automatic_guess=True)  # Example initial guess values

# Perform fitting
dtheta, dbeta, E_opt = c.fit_parameters(measured_energies, cff_values, orders)

c.print_fit_results()
# Plot the fitting results along with initial guess
c.plot_fit(measured_energies, cff_values, orders,
           DTheta=dtheta, DBeta=dbeta, E_opt=E_opt,
           show=False, savepath='Results/UE112_2013/fit_results_ue112_2013.png', 
           plot_initial_guess=True)




# Rolf's results 
rolf_dtheta = -np.rad2deg(8.23621e-05)
rolf_dbeta  = np.rad2deg(8.67778e-05)
en_rolf     = 91.1254


# Calculate the percentage difference
diff_dtheta = 100 * (dtheta - rolf_dtheta) / rolf_dtheta
diff_dbeta  = 100 * (dbeta - rolf_dbeta) / rolf_dbeta
diff_e_opt  = 100 * (E_opt - en_rolf) / en_rolf

# Create the table
headers = ["Parameter", "Fit Results", "Rolf's Results", "Difference"]
data = [
    ["DTheta [deg]", f"{dtheta:.6f}", f"{rolf_dtheta:.6f}", f"{diff_dtheta:.2f}%"],
    ["DBeta  [deg]", f"{dbeta:.6f}", f"{rolf_dbeta:.6f}", f"{diff_dbeta:.2f}%"],
    ["E_opt  [eV]", f"{E_opt:.6f}", f"{en_rolf:.6f}", f"{diff_e_opt:.2f}%"]
]

# Print the table
print('\nComparison with Rolfs Results')
print(tabulate(data, headers, tablefmt="pretty", colalign=("left", "right", "right", "right")))



# EMILE Results
e_dtheta = -np.rad2deg(4.89892e-05)
e_dbeta  = np.rad2deg(0.000115096)
en_e     = 91.1447

# Calculate the percentage difference
e_diff_dtheta = 100 * (dtheta - e_dtheta) / e_dtheta
e_diff_dbeta  = 100 * (dbeta - e_dbeta) / e_dbeta
e_diff_e_opt  = 100 * (E_opt - en_e) / en_e

# Create the table
headers = ["Parameter", "Fit Results", "EMILE's Results", "Difference"]
data = [
    ["DTheta [deg]", f"{dtheta:.6f}", f"{e_dtheta:.6f}",
     f"{e_diff_dtheta:.2f}%"],
    ["DBeta  [deg]", f"{dbeta:.6f}", f"{e_dbeta:.6f}", 
     f"{e_diff_dbeta:.2f}%"],
    ["E_opt  [eV]", f"{E_opt:.6f}", f"{en_e:.6f}", 
     f"{e_diff_e_opt:.4f}%"]
]

# Print the table
print('\nComparison with Emiles Results')
print(tabulate(data, headers, tablefmt="pretty", colalign=("left", "right", "right", "right")))

