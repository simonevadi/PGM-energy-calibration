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

        # indicates if the initial guess should be done by the program
        self.automatic_guess = False

    def set_initial_guess(self, DTheta=0, DBeta=0, E=None, automatic_guess=False):
        """
        Set the initial guess for the fitting parameters.
        
        Args:
            DTheta (float): Initial guess for delta theta in deg
            DBeta (float): Initial guess for delta beta in deg
            E (float): Initial guess for energy in eV
        """
        if automatic_guess:
            self.automatic_guess = True
            return
        elif not automatic_guess and E is not None:
            self.initial_guess = {'DTheta': np.deg2rad(DTheta), 'DBeta': np.deg2rad(DBeta), 'E': E}
        elif not automatic_guess and E is None:
            raise ValueError('If automatic guess is False, E must be given')
    
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

    def calc_beta(self, alpha, cff):
        """
        Calculate beta angle.

        Args:
            wavelength (float): Wavelength of the light.
            cff (float): Fixed focus constant.

        Returns:
            float: Calculated beta angle in rad.
        """
        beta  = abs(-1 * np.arccos(cff * np.cos(alpha)))

        return beta

    def calc_alpha(self, cff, order, E):
        """
        Calculate alpha angle.

        Args:
            beta (float): Beta angle.
            cff (float): Fixed focus constant.

        Returns:
            float: Calculated alpha angle in rad.
        """
        true_wavelength = self.hc*order/E
        temp   = (cff * cff - 1) / (true_wavelength * self.N)
        w      = cff * cff + temp * temp
        x      = true_wavelength * self.N / (cff * cff - 1) * (np.sqrt(w) - 1)
        alpha  = np.arcsin(x)

        return alpha

    def calc_theta(self, alpha, beta):
        """
        Calculate theta angle.

        Args:
            alpha (float): Alpha angle.
            beta (float): Beta angle.

        Returns:
            float: Calculated theta angle in rad.
        """
        theta = (beta + alpha) / 2
        return theta
    

    def calc_wavelength(self, energy):
        """
        Calculate wavelength from energy.

        Args:
            energy (float): Energy in eV.

        Returns:
            float: Calculated wavelength in meter.
        """
        wavelength = self.hc / energy
        return wavelength

    def calc_energy(self, wavelength):
        """
        Calculate energy from wavelength.

        Args:
            wavelength (float): Wavelength of the light.

        Returns:
            float: Calculated energy in eV.
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
        # true_wavelength = self.calc_wavelength(true_energy)
        alpha = self.calc_alpha(cff, order, true_energy)
        beta = self.calc_beta(alpha, cff)
        theta = self.calc_theta(alpha, beta)
        es = self.grating_equation(np.deg2rad(90)-theta, np.deg2rad(90)-beta, order, dtheta=dtheta, dbeta=dbeta)
        # print(f'E={true_energy}, cff={cff}, alpha={np.rad2deg(alpha)}, theta={90-np.rad2deg(theta)}, beta={90-np.rad2deg(beta)}, Ecal={es}')
        return es

    def _check_angles(self, E,shifted_energy, cff,order, DTheta, DBeta):
        true_wavelength = self.calc_wavelength(E)
        alpha = self.calc_alpha(cff, order, shifted_energy)
        beta = self.calc_beta(alpha, cff)
        theta = self.calc_theta(alpha, beta)
        shifted_energy = self.grating_equation(theta, beta, order)#, dtheta=DTheta, dbeta=DBeta)
        print(f'E={E}, w={true_wavelength}, cff={cff}, alpha={np.rad2deg(alpha)},  theta={np.rad2deg(theta)}, beta={np.rad2deg(beta)}, Ecal={shifted_energy}')

    def _check_angles_from_paper(self, E,shifted_energy, cff,order, DTheta, DBeta):
        true_wavelength = self.hc*order/E
        temp   = (cff * cff - 1) / (true_wavelength * self.N)
        w      = cff * cff + temp * temp
        x      = true_wavelength * self.N / (cff * cff - 1) * (np.sqrt(w) - 1)
        alpha  = np.arcsin(x)
        beta  = abs(-1 * np.arccos(cff * np.cos(alpha)))
        theta = 0.5 * (alpha + beta)
        shifted_energy = self.grating_equation(np.deg2rad(90)-theta, np.deg2rad(90)-beta, order)#, dtheta=DTheta, dbeta=DBeta)
        print(f'E={E}, w={true_wavelength}, cff={cff}, alpha={90-np.rad2deg(alpha)}, theta={90-np.rad2deg(theta)}, beta={90-np.rad2deg(beta)}, Ecal={shifted_energy}')
        print('\n')

    def residuals(self, params, measured_energies, cff_values, orders):
        """
        Calculate residuals between measured and shifted energies with higher weights for higher orders.

        Args:
            params (list): List of parameters [DTheta, DBeta, E].
            measured_energies (list): List of measured energies.
            cff_values (list): List of fixed focus constant values.
            orders (list): List of diffraction orders.

        Returns:
            list: List of weighted residuals.
        """
        DTheta, DBeta, E = params
        residuals = []
        for i, cff in enumerate(cff_values):
            measured_energy = measured_energies[i]
            shifted_energy = self.shifted_energy(cff, orders[i], DTheta, DBeta, E)
            residuals.append((measured_energy - shifted_energy))
        return residuals
    
    def fit_parameters(self, measured_energies, cff_values, orders, 
                       print_fit_report=True,
                       return_fit_eval=False):
        """
        Fit parameters using least squares optimization.

        Args:
            measured_energies (list): List of measured energies.
            cff_values (list): List of fixed focus constant values.
            delta_E (int,float): maximum deviation of the energy from the initial guess parameter

        Returns:
            tuple: Optimized parameters (DTheta, DBeta, E_opt).
        """
        orders            = np.array(orders)
        cff_values        = np.array(cff_values)
        measured_energies = np.array(measured_energies)
        
        # Check if a value in cff_values is strictly lower than one
        mask = cff_values < 1
        # If the condition is met, make the corresponding value in orders negative
        orders[mask] *= -1
        print(orders)

        if self.automatic_guess:
            self.set_initial_guess(DTheta=0, DBeta=0, E=np.mean(measured_energies))
       
        initial_guess_list = [self.initial_guess['DTheta'], self.initial_guess['DBeta'], self.initial_guess['E']]
        
        result = least_squares(self.residuals, initial_guess_list, 
                               args=(measured_energies, cff_values, orders),
                               method='lm')
        
        DTheta, DBeta, self.E_opt = result.x
        # self.DTheta, self.DBeta, self.E_opt = result.x

        # convert to degrees
        self.DTheta = np.rad2deg(DTheta)
        self.DBeta = np.rad2deg(DBeta)
        
        if print_fit_report:
            self.print_fit_report(result)
        if return_fit_eval:
            return self.DTheta, self.DBeta, self.E_opt, result.cost, result.optimality, result.nfev
        return self.DTheta, self.DBeta, self.E_opt

    def print_fit_report(self, result):
        """
        Print the goodness of fit report.

        Args:
            result: The result object from the least_squares optimization.
        """
        headers = ["Metric", "Value"]
        residuals_list = [[res] for res in result.fun]
        data = [
            ["Cost (Sum of squared residuals)", f"{result.cost:.6f}"],
            ["Optimality (First-order optimality measure)", f"{result.optimality:.6f}"],
            ["Number of function evaluations", result.nfev],
            ["Residuals", residuals_list]
        ]
        
        # Print the table
        print(tabulate(data, headers, tablefmt="pretty", colalign=("left", "left")))


    def print_fit_results(self):
        """
        Print the fit results in a formatted table.
        """
        # Convert to rad and round to two decimals
        DTheta_urad = np.round(np.deg2rad(self.DTheta) * 1E6, 2)
        DBeta_urad = np.round(np.deg2rad(self.DBeta)* 1E6, 2)

        # Round to two decimals
        DTheta_deg = np.round(self.DTheta, 6)
        DBeta_deg  = np.round(self.DBeta, 6)
        
        E_opt_rounded = np.round(self.E_opt, 6)

        # Prepare data for tabulate
        headers = ["Parameter", "Initial Guess", "Value", 'Value']
        data = [
            ["DTheta ",f"{self.initial_guess['DTheta']:>10} [deg]", f"{DTheta_urad:>10} [urad]", f"{DTheta_deg:>10} [deg]"],
            ["DBeta", f"{self.initial_guess['DBeta']:>10} [deg]",f"{DBeta_urad:>10} [urad]", f"{DBeta_deg:>10} [deg]"],
            ["E_opt  ", f"{np.around(self.initial_guess['E'], 2):>10} [eV]", f"{E_opt_rounded:>10} [eV]", f"{E_opt_rounded:>10} [eV]"]
        ]

        # Print the table with a title
        print("\nFit Results\n")
        print(tabulate(data, headers, tablefmt="pretty", colalign=("left", "right", 'right', 'right')))

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
        # colors = plt.cm.viridis(np.linspace(0, 1, len(unique_orders)))
        colors = ['k', 'royalblue', 'darkgoldenrod', 'g', 'magenta', 'orange', 'red']
        plt.figure()
        
        for idx, order in enumerate(unique_orders):
            # Filter data for the current order
            mask = (orders == order)
            cff_order = cff_values[mask]
            measured_order = measured_energies[mask]
            
            # use more points to plot the fit results
            cff_plot = np.arange(np.min(cff_values), np.max(cff_values), .01)
            mask = (cff_plot < 0.9) | (cff_plot > 1.1)
            cff_plot = cff_plot[mask]

            mask = cff_plot < 1
            orders_plot = np.full(cff_plot.shape, order)
            orders_plot[mask] *= -1

            orders_ig = np.full(cff_order.shape, order)
            mask = cff_order < 1
            orders_ig[mask] *= -1

            fitted_energies = self.shifted_energy(cff_plot, 
                                                orders_plot, 
                                                np.deg2rad(DTheta), 
                                                np.deg2rad(DBeta), 
                                                E_opt )       

            # Calculate initial guess energies for the current order
            initial_guess_energies = self.shifted_energy(cff_order,
                                                        orders_ig,
                                                        np.deg2rad(self.initial_guess['DTheta']),
                                                        np.deg2rad(self.initial_guess['DBeta']),
                                                        self.initial_guess['E']) 
            
            plt.scatter(cff_order, measured_order , c=colors[idx], label=f'Measured Order {order}', alpha=0.7)
            if plot_initial_guess:
                plt.plot(cff_order, initial_guess_energies, '--', color=colors[idx], label=f'Initial Guess Order {order}', alpha=0.7)
            plt.plot(cff_plot, fitted_energies, '-', color=colors[idx], label=f'Fitted Order {order}', alpha=0.7)

        plt.xlabel('$c_{ff}$')
        plt.ylabel('Energy (eV)')
        plt.xscale('log')
        
        xmin = np.min(cff_values)
        xmax = np.max(cff_values)
        plt.xlim((xmin-.01, xmax+.2))
        plt.ylim((np.min(measured_energies)-1, np.max(measured_energies)+1))

        
        
        xticks_positions, xticks_labels = self.generate_x_ticks_pos_and_label(xmin, xmax)
    
        plt.xticks(xticks_positions, labels=xticks_labels)

        plt.legend()
        plt.tight_layout()
        if savepath is not None:
            plt.savefig(savepath)
        if show:
            plt.show()
        return plt
    
    def generate_x_ticks_pos_and_label(self, xmin, xmax):
        xticks_positions = []
        xticks_labels = []
        # Generate x-tick positions and labels
        for ind, i in enumerate(np.arange(xmin, xmax, 0.1)):
            if i < 1 and ind%2==0:
                xticks_positions.append(i)
                xticks_labels.append(f"{i:.1f}")
            else:
                
                i = np.round(i,2)
                if i.is_integer():
                    xticks_positions.append(i)
                    xticks_labels.append(f"{int(i)}")
        return xticks_positions, xticks_labels

