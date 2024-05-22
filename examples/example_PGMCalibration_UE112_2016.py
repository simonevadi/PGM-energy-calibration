import os
import numpy as np
from tabulate import tabulate
from PGM_calibration import PGMCalibration

N = 600  # Line density (grating grooves per mm)
# k = 1  # Order of diffraction

os.makedirs('Results/ExampleCalibration', exist_ok=True)

# Example data
orders = np.array([           1,      1,     1,      2,      2,      2,      3,      3,      3,      4,      4,      4])
cff_values = np.array([       1.5,    3,     8,      1.5,    3,      8,      1.5,    3,      8,      1.5,    3,      8])
measured_energies = np.array([92.819, 91.84, 91.416, 92.295, 91.641, 91.341, 92.044, 91.545, 91.309, 91.913, 91.489, 91.284 ])

# Select only orders 3 and 4
selected_orders = [1]
mask = np.isin(orders, selected_orders)

# Filtered data
orders = orders[mask]
cff_values = cff_values[mask]
measured_energies = measured_energies[mask]




# Instantiate the calibration class
c = PGMCalibration(N)

# Set a reasonable initial guess manually
c.set_initial_guess(DTheta=0.0, DBeta=0.0, E=91.84)  # Example initial guess values

# Perform fitting
dtheta, dbeta, E_opt = c.fit_parameters(measured_energies, cff_values, orders)

c.print_fit_results()
# Plot the fitting results along with initial guess
c.plot_fit(measured_energies, cff_values, orders,
           DTheta=dtheta, DBeta=dbeta, E_opt=E_opt,
           show=False, savepath='Results/UE112_2016/fit_results_ue112_2016.png', 
           plot_initial_guess=False)



# Convert dtheta and dbeta from radians to degrees
dtheta_deg = np.rad2deg(dtheta)
dbeta_deg = np.rad2deg(dbeta)

# Rolf's results
rolf_dtheta = 0.014690
rolf_dbeta = 0.024430
en_rolf = 91.008522

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

print(np.rad2deg(-0.000519061))