import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

class PGMCalibration:
    def __init__(self, N, k, h=4.135667696e-15, c=299792458):
        """
        Initialize the PGMCalibration class with given parameters.

        Parameters:
        N (int): Line density (grating grooves per mm)
        k (int): Order of diffraction
        h (float): Planck's constant (eV s)
        c (float): Speed of light (m/s)
        """
        self.N = N
        self.k = k
        self.h = h
        self.c = c

    def grating_equation(self, wavelength, a, b):
        return self.N * self.k * wavelength - 2 * np.cos(a) * np.sin(a + b)

    def calc_b(self, cff, a):
        return -np.arccos(cff * np.cos(a))

    def calc_a(self, l, cff):
        term = np.sqrt(cff**2 + (cff**2 - 1)**2 * (l * self.N * self.k)**2) - 1
        return np.arcsin(l * self.N * self.k / (cff**2 - 1) * term)

    def shifted_wavelength(self, cff, DY, Db, true_energy):
        true_wavelength = self.h * self.c / true_energy
        a = self.calc_a(true_wavelength, cff)
        b = self.calc_b(cff, a)
        Y = (a + b) / 2
        shifted_lambda = 2 * np.cos(Y + DY) * np.sin(Y + DY + b + Db) / (self.N * self.k)
        return shifted_lambda

    def residuals(self, params, measured_energies, cff_values):
        DY, Db, E = params
        residuals = []
        for i, cff in enumerate(cff_values):
            measured_lambda = self.h * self.c / measured_energies[i]
            shifted_lambda = self.shifted_wavelength(cff, DY, Db, E)
            residuals.append(measured_lambda - shifted_lambda)
        return residuals

    def fit_parameters(self, measured_energies, cff_values):
        initial_guess = [0.0, 0.0, measured_energies[0]]
        result = least_squares(self.residuals, initial_guess, args=(measured_energies, cff_values))
        self.DY, self.Db, self.E_opt = result.x
        return self.DY, self.Db, self.E_opt

    def plot_fit(self, measured_energies, cff_values):
        fitted_energies = [self.h * self.c / self.shifted_wavelength(cff, self.DY, self.Db, self.E_opt) for cff in cff_values]
        plt.plot(cff_values, measured_energies, 'ro', label='Measured energies')
        plt.plot(cff_values, fitted_energies, 'b-', label='Fitted energies')
        plt.xlabel('cff')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.show()

# Example usage
if __name__ == "__main__":
    N = 1200  # Line density (grating grooves per mm)
    k = 1  # Order of diffraction

    # Example data
    measured_energies = np.array([400.77, 400.76, 400.89])  # Measured energies at different cff values
    cff_values = np.array([2.1, 2.3, 2.5])  # Corresponding cff values

    # Instantiate the calibration class and perform fitting
    calibration = PGMCalibration(N, k)
    DY, Db, E_opt = calibration.fit_parameters(measured_energies, cff_values)

    # Print fitting results
    print(f"Fitted results:\nDY = {DY:.6f} rad\nDb = {Db:.6f} rad\nE_opt = {E_opt:.2f} eV")

    # Plot the fitting results
    calibration.plot_fit(measured_energies, cff_values)
