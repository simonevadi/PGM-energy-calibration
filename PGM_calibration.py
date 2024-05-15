import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from tabulate import tabulate

class PGMCalibration:
    def __init__(self, N, k):
        """
        Initialize the PGMCalibration class with given parameters.

        Args:
            N (int): Line density (grating grooves per mm)
            k (int): Order of diffraction
        """
        self.N = N * 1000            # lines/m
        self.k = k
        self.h = 4.135667696E-15     # eV*sec
        self.c = 299792458           # m/s 
        self.hc = self.h * self.c
        self.initial_guess = {'DTheta': 0.0, 'DBeta': 0.0, 'E': 0.0}  # Default initial guess

        # Variables that will hold the fit results
        self.DTheta = None
        self.DBeta = None
        self.E_opt = None

    def set_initial_guess(self, DTheta, DBeta, E):
        """
        Set the initial guess for the fitting parameters.
        
        Args:
            DTheta (float): Initial guess for delta theta.
            DBeta (float): Initial guess for delta beta.
            E (float): Initial guess for energy.
        """
        self.initial_guess = {'DTheta': DTheta, 'DBeta': DBeta, 'E': E}

    def grating_equation(self, theta, beta, dtheta=0, dbeta=0):
        """
        Calculate energy using the grating equation.

        Args:
            wavelength (float): Wavelength of the light.
            theta (float): Angle of incidence.
            beta (float): Angle of diffraction.
            dtheta (float, optional): Delta theta. Defaults to 0.
            dbeta (float, optional): Delta beta. Defaults to 0.

        Returns:
            float: Calculated energy.
        """
        theta = np.deg2rad(90) - (theta + dtheta)
        beta = -(np.deg2rad(90) - (beta + dbeta))
        wavelength = (2 * np.cos(theta) * np.sin(theta + beta)) / (self.N * self.k)
        en = self.calc_energy(wavelength)
        return en

    def calc_beta(self, wavelength, cff):
        """
        Calculate beta angle.

        Args:
            wavelength (float): Wavelength of the light.
            cff (float): Fixed focus constant.

        Returns:
            float: Calculated beta angle.
        """
        a = 2 * self.N * self.k * wavelength
        b = cff**2 / (cff**2 - 1)
        beta = np.sqrt(a * b)
        return beta

    def calc_alpha(self, beta, cff):
        """
        Calculate alpha angle.

        Args:
            beta (float): Beta angle.
            cff (float): Fixed focus constant.

        Returns:
            float: Calculated alpha angle.
        """
        a = np.sin(beta) / cff
        alpha = np.arcsin(a)
        return alpha

    def calc_theta(self, alpha, beta):
        """
        Calculate theta angle.

        Args:
            alpha (float): Alpha angle.
            beta (float): Beta angle.

        Returns:
            float: Calculated theta angle.
        """
        theta = (beta + alpha) / 2
        return theta

    def calc_wavelength(self, energy):
        """
        Calculate wavelength from energy.

        Args:
            energy (float): Energy in eV.

        Returns:
            float: Calculated wavelength.
        """
        wavelength = self.hc / energy
        return wavelength

    def calc_energy(self, wavelength):
        """
        Calculate energy from wavelength.

        Args:
            wavelength (float): Wavelength of the light.

        Returns:
            float: Calculated energy.
        """
        energy = self.hc / wavelength
        return energy

    def shifted_energy(self, cff, dtheta, dbeta, true_energy):
        """
        Calculate shifted energy.

        Args:
            cff (float): Fixed focus constant.
            dtheta (float): Delta theta.
            dbeta (float): Delta beta.
            true_energy (float): True energy.

        Returns:
            float: Shifted energy.
        """
        true_wavelength = self.calc_wavelength(true_energy)
        beta = self.calc_beta(true_wavelength, cff)
        alpha = self.calc_alpha(beta, cff)
        theta = self.calc_theta(alpha, beta)
        es = self.grating_equation(theta, beta, dtheta=dtheta, dbeta=dbeta)
        return es

    def residuals(self, params, measured_energies, cff_values):
        """
        Calculate residuals between measured and shifted energies.

        Args:
            params (list): List of parameters [DTheta, DBeta, E].
            measured_energies (list): List of measured energies.
            cff_values (list): List of fixed focus constant values.

        Returns:
            list: List of residuals.
        """
        DTheta, DBeta, E = params
        residuals = []
        for i, cff in enumerate(cff_values):
            measured_energy = measured_energies[i]
            shifted_energy = self.shifted_energy(cff, DTheta, DBeta, E)
            residuals.append(measured_energy - shifted_energy)
        return residuals

    def fit_parameters(self, measured_energies, cff_values):
        """
        Fit parameters using least squares optimization.

        Args:
            measured_energies (list): List of measured energies.
            cff_values (list): List of fixed focus constant values.

        Returns:
            tuple: Optimized parameters (DTheta, DBeta, E_opt).
        """
        initial_guess_list = [self.initial_guess['DTheta'], self.initial_guess['DBeta'], self.initial_guess['E']]
        result = least_squares(self.residuals, initial_guess_list, args=(measured_energies, cff_values))
        self.DTheta, self.DBeta, self.E_opt = result.x
        return self.DTheta, self.DBeta, self.E_opt

    def print_fit_results(self):
        """
        Print the fit results in a formatted table.
        """
        # Convert to urad and round to two decimals
        DTheta_urad = np.round(self.DTheta * 1E6, 2)
        DBeta_urad = np.round(self.DBeta * 1E6, 2)
        E_opt_rounded = np.round(self.E_opt, 2)

        # Prepare data for tabulate
        headers = ["Parameter", "Value"]
        data = [
            ["DTheta [urad]", f"{DTheta_urad:>10}"],
            ["DBeta  [urad]", f"{DBeta_urad:>10}"],
            ["E_opt  [eV]", f"{E_opt_rounded:>10}"]
        ]

        # Print the table with a title
        print("\nFit Results\n")
        print(tabulate(data, headers, tablefmt="pretty", colalign=("left", "right")))

    def plot_fit(self, measured_energies, cff_values, DTheta=0, DBeta=0, E_opt=0, savepath=None, show=True):
        """
        Plot the measured, initial guess, and fitted energies.

        Args:
            measured_energies (list): List of measured energies.
            cff_values (list): List of fixed focus constant values.
            DTheta (float, optional): Delta theta for fit. Defaults to 0.
            DBeta (float, optional): Delta beta for fit. Defaults to 0.
            E_opt (float, optional): Optimized energy. Defaults to 0.
        """
        # Calculate fitted energies
        cff_plot = np.arange(cff_values[0], cff_values[-1], .1)
        fitted_energies = self.shifted_energy(cff_values, DTheta, DBeta, E_opt)        
        # Calculate initial guess energies
        initial_guess_energies = self.shifted_energy(cff_values,
                                                     self.initial_guess['DTheta'],
                                                     self.initial_guess['DBeta'],
                                                     self.initial_guess['E'])
        
        plt.plot(cff_values, measured_energies, 'ro', label='Measured energies')
        plt.plot(cff_values, initial_guess_energies, 'g', label='Initial guess')
        plt.plot(cff_values, fitted_energies, 'b-', label='Fitted energies')
        plt.xlabel('$c_{ff}$')
        plt.ylabel('Energy (eV)')
        plt.xscale('log')
        plt.legend()
        if savepath is not None:
            plt.savefig(savepath)
        if show:
            plt.show()
