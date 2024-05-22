import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from tabulate import tabulate

class PGMCalibration:
    def __init__(self, N):
        """
        Initialize the PGMCalibration class with given parameters.

        Args:
            N (int): Line density (grating grooves per mm)
            k (int): Order of diffraction
        """
        self.N = N * 1000            # lines/m
        # self.k = k
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
            DTheta (float): Initial guess for delta theta in rad
            DBeta (float): Initial guess for delta beta in rad
            E (float): Initial guess for energy in eV
        """
        self.initial_guess = {'DTheta': DTheta, 'DBeta': DBeta, 'E': E}

    def grating_equation(self, theta, beta, k, dtheta=0, dbeta=0):
        """
        Calculate energy using the grating equation. Angles are grazing and in degrees

        Args:
            wavelength (float): Wavelength of the light.
            theta (float): Grazing angle on the mirror in deg
            beta (float): Grazing exit angle grating in deg
            dtheta (float, optional): Delta theta. Defaults to 0.
            dbeta (float, optional): Delta beta. Defaults to 0.

        Returns:
            float: Calculated energy in eV
        """
        theta = np.deg2rad(90) - (theta + dtheta)
        beta = -(np.deg2rad(90) - (beta + dbeta))
        wavelength = (2 * np.cos(theta) * np.sin(theta + beta)) / (self.N * k)
        en = self.calc_energy(wavelength)
        return en

    def calc_beta(self, wavelength, cff, order):
        """
        Calculate beta angle.

        Args:
            wavelength (float): Wavelength of the light.
            cff (float): Fixed focus constant.

        Returns:
            float: Calculated beta angle.
        """
        a = 2 * self.N * order * wavelength
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

    def shifted_energy(self, cff, order, dtheta, dbeta, true_energy):
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
        beta = self.calc_beta(true_wavelength, cff, order)
        alpha = self.calc_alpha(beta, cff)
        theta = self.calc_theta(alpha, beta)
        es = self.grating_equation(theta, beta, order, dtheta=dtheta, dbeta=dbeta)
        return es

    def residuals(self, params, measured_energies, cff_values, orders):
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
            shifted_energy = self.shifted_energy(cff, orders[i], DTheta, DBeta, E)
            residuals.append(measured_energy - shifted_energy)
        return residuals
    
    def check_angles(self, E,shifted_energy, cff,order):
        true_wavelength = self.calc_wavelength(E)
        beta = self.calc_beta(true_wavelength, cff, order)
        alpha = self.calc_alpha(beta, cff)
        theta = self.calc_theta(alpha, beta)
        print(f'E={E}, w={true_wavelength}, cff={cff},  theta={90-np.rad2deg(theta)}, beta={90-np.rad2deg(beta)}, Ecal={shifted_energy}')

    def fit_parameters(self, measured_energies, cff_values, orders, print_fit_report=True, delta_E=3):
        """
        Fit parameters using least squares optimization.

        Args:
            measured_energies (list): List of measured energies.
            cff_values (list): List of fixed focus constant values.
            delta_E (int,float): maximum deviation of the energy from the initial guess parameter

        Returns:
            tuple: Optimized parameters (DTheta, DBeta, E_opt).
        """
        initial_guess_list = [self.initial_guess['DTheta'], self.initial_guess['DBeta'], self.initial_guess['E']]
        
        # Define bounds for the parameters
        lower_bounds = [-np.inf, -np.inf, self.initial_guess['E'] - delta_E]
        upper_bounds = [np.inf, np.inf, self.initial_guess['E'] + delta_E]
        
        result = least_squares(self.residuals, initial_guess_list, args=(measured_energies, cff_values, orders),
                            #bounds=(lower_bounds, upper_bounds),
                            method='lm')
        
        self.DTheta, self.DBeta, self.E_opt = result.x
        if print_fit_report:
            self.print_fit_report(result)
        return self.DTheta, self.DBeta, self.E_opt

    def print_fit_report(self, result):
        """
        Print the goodness of fit report.

        Args:
            result: The result object from the least_squares optimization.
        """
        headers = ["Metric", "Value"]
        data = [
            ["Cost (Sum of squared residuals)", f"{result.cost:.6f}"],
            ["Optimality (First-order optimality measure)", f"{result.optimality:.6f}"],
            ["Number of function evaluations", result.nfev],
            ["Residuals", result.fun]
        ]
        
        # Print the table
        print(tabulate(data, headers, tablefmt="pretty", colalign=("left", "left")))


    def print_fit_results(self):
        """
        Print the fit results in a formatted table.
        """
        # round to two decimals
        DTheta_urad = np.round(self.DTheta * 1E6, 2)
        DBeta_urad = np.round(self.DBeta * 1E6, 2)
        E_opt_rounded = np.round(self.E_opt, 2)

        # Convert to deg and round to two decimals
        DTheta_deg = np.round(np.rad2deg(self.DTheta), 6)
        DBeta_deg = np.round(np.rad2deg(self.DBeta), 6)
        E_opt_rounded = np.round(self.E_opt, 6)

        # Prepare data for tabulate
        headers = ["Parameter", "Value", 'Value']
        data = [
            ["DTheta ", f"{DTheta_urad:>10} [urad]", f"{DTheta_deg:>10} [deg]"],
            ["DBeta", f"{DBeta_urad:>10} [urad]", f"{DBeta_deg:>10} [deg]"],
            ["E_opt  ", f"{E_opt_rounded:>10} [eV]", f"{E_opt_rounded:>10} [eV]"]
        ]

        # Print the table with a title
        print("\nFit Results\n")
        print(tabulate(data, headers, tablefmt="pretty", colalign=("left", "right", 'right')))

    def plot_fit(self, measured_energies, cff_values, 
                 orders, DTheta=0, DBeta=0, E_opt=0, 
                 savepath=None, show=True,
                 plot_initial_guess=True):
        """
        Plot the measured, initial guess, and fitted energies.

        Args:
            measured_energies (list): List of measured energies.
            cff_values (list): List of fixed focus constant values.
            orders (list): List of diffraction orders.
            DTheta (float, optional): Delta theta for fit. Defaults to 0.
            DBeta (float, optional): Delta beta for fit. Defaults to 0.
            E_opt (float, optional): Optimized energy. Defaults to 0.
            savepath (str, optional): Path to save the plot. Defaults to None.
            show (bool, optional): Whether to display the plot. Defaults to True.
            plot_initial_guess (bool, optional): Whether to plot the initial guess. Defaults to True.
        """
        unique_orders = np.unique(orders)
        colors = plt.cm.viridis(np.linspace(0, 1, len(unique_orders)))
        
        plt.figure()
        
        for idx, order in enumerate(unique_orders):
            # Filter data for the current order
            mask = (orders == order)
            cff_order = cff_values[mask]
            measured_order = measured_energies[mask]
            
            # Calculate fitted energies for the current order
            cff_plot = np.arange(cff_order[0], cff_order[-1], .1)
            fitted_energies = self.shifted_energy(cff_plot, 
                                                np.full(cff_plot.shape, order), 
                                                DTheta, 
                                                DBeta, 
                                                E_opt )       

            # Calculate initial guess energies for the current order
            initial_guess_energies = self.shifted_energy(cff_order,
                                                        np.full(cff_order.shape, order),
                                                        self.initial_guess['DTheta'],
                                                        self.initial_guess['DBeta'],
                                                        self.initial_guess['E']) 
            
            plt.scatter(cff_order, measured_order , color=colors[idx], label=f'Measured Order {order}', alpha=0.7)
            if plot_initial_guess:
                plt.plot(cff_order, initial_guess_energies, '--', color=colors[idx], label=f'Initial Guess Order {order}', alpha=0.7)
            plt.plot(cff_plot, fitted_energies, '-', color=colors[idx], label=f'Fitted Order {order}', alpha=0.7)

        plt.xlabel('$c_{ff}$')
        plt.ylabel('Energy (eV)')
        plt.xscale('log')
        
        # Get current tick positions and labels
        x_ticks = plt.gca().get_xticks()
        y_ticks = plt.gca().get_yticks()
        
        # Create custom tick labels
        x_labels = [f'{tick:.2f}' for tick in x_ticks]
        y_labels = [f'{tick:.2f}' for tick in y_ticks]
        
        # Set the custom tick labels
        plt.gca().set_xticks(x_ticks)
        plt.gca().set_xticklabels(x_labels)
        plt.gca().set_yticks(y_ticks)
        plt.gca().set_yticklabels(y_labels)

        plt.xlim((cff_order[0]-.1, cff_order[-1]+1))
        plt.legend()
        plt.tight_layout()
        if savepath is not None:
            plt.savefig(savepath)
        if show:
            plt.show()
        return plt

