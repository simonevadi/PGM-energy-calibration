import numpy as np
from tabulate import tabulate
from PGM_calibration import PGMCalibration

N = 1200  # Line density (grating grooves per mm)
k = 1  # Order of diffraction

# Instantiate the calibration class
c = PGMCalibration(N, k)

# test the energy calculation back and forth
energy = 400
l      = c.calc_wavelength(energy)
beta   = c.calc_beta(l,2.25)
alpha  = c.calc_alpha(beta,2.25)
theta  = c.calc_theta(alpha, beta)
en     = c.grating_equation(theta, beta)
diff_en = np.round(abs(en-energy)*1000)


# Create data for the table
data = [
    ["Wavelength [m]", l],
    ["Alpha [degrees]", np.round(np.rad2deg(alpha), 2)],
    ["Beta [degrees]", np.round(np.rad2deg(beta), 2)],
    ["Theta [degrees]", np.round(np.rad2deg(theta), 2)],
    ["Initial Energy [eV]", np.round(energy, 3)],
    ["Recalculated Energy [eV]", np.round(en, 3)],
    ["Difference [meV]", diff_en]
]

# Print the table with grid format
print(tabulate(data, headers=["Parameter", "Value"], tablefmt="grid", colalign=("left", "right")))
