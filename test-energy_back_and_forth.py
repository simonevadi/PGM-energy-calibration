import numpy as np

from PGM_calibration import PGMCalibration

N = 1200  # Line density (grating grooves per mm)
k = 1  # Order of diffraction

# Instantiate the calibration class
c = PGMCalibration(N, k)

# test the energy calculation back and forth
energy = 400
l = c.calc_wavelength(energy)
beta = c.calc_beta(l,2.25)
alpha = c.calc_alpha(beta,2.25)
theta = c.calc_theta(alpha, beta)
en = c.grating_equation(l, theta, beta)
print(f'wavelength: {l},\n\
alpha {np.rad2deg(alpha)},\n\
beta {np.rad2deg(beta)},\n\
theta {np.rad2deg(theta)}')
diff_en = np.round(abs(en-energy)*1000)
print(f'initial energy: {np.round(energy,3)}eV, recalculated energy: {np.round(en,3)}eV, difference: {diff_en} meV')
