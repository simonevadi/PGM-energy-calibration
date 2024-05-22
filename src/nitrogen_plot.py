import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator, FormatStrFormatter
import numpy as np

class PlotNitrogenFit:

    def plot_results(self, energy, data, out, delta, fix_param, dict_fit, residuals_final, RMS_neg, RMS_pos, show=False, save_img=False):
        """
        Plot the fit results.

        Args:
            energy (np.array): Array of photon energy values.
            data (np.array): Array of ion counts.
            out (ModelResult): Model result from the fit.
            delta (np.array): Fit uncertainty values.
            fix_param (bool): Whether fixed parameters were used in the fit.
            dict_fit (dict): Dictionary of fit parameters.
            residuals_final (np.array): Residuals between data and best fit.
            RMS_neg (float): RMS value for initial guess residuals.
            RMS_pos (float): RMS value for final residuals.
            show (bool, optional): Whether to display the plot. Defaults to False.
            save_img (bool or str, optional): If False, no image is saved. If a string, saves the image with this name. Defaults to False.
        """
        plt.rc("font", size=12, family='serif')
        fig, axes = plt.subplots(3, 1, figsize=(8.0, 16.0))

        # Plot Initial Guess, Data, Fit, and Fit Uncertainty
        self.plot_initial_guess_and_fit(axes[0], energy, data, out, delta)

        # Plot Voigt and Linear Components
        self.plot_voigt_and_linear_components(axes[1], energy, data, out, fix_param, dict_fit)

        # Plot Residuals and Confidence Interval
        self.plot_residuals_and_confidence_interval(axes[2], energy, data, residuals_final, RMS_neg, RMS_pos)

        plt.tight_layout()
        if save_img:
            plt.savefig(f'{save_img}')
        if show:
            plt.show()

    def plot_initial_guess_and_fit(self, ax, energy, data, out, delta):
        """
        Plot initial guess, data, fit, and fit uncertainty.

        Args:
            ax (Axes): Matplotlib Axes object.
            energy (np.array): Array of photon energy values.
            data (np.array): Array of ion counts.
            out (ModelResult): Model result from the fit.
            delta (np.array): Fit uncertainty values.
        """
        ax.plot(energy, out.init_fit, 'orange', label='initial guess')
        ax.scatter(energy, data, label='data', s=10)
        ax.plot(energy, out.best_fit, 'r', label='best fit')
        ax.fill_between(energy, out.best_fit - delta, out.best_fit + delta, color='gray', alpha=0.4, label='fit uncertainty')
        self.configure_plot(ax, 
                            'Photon energy (eV)',
                            'Ion counts (arb. units)',
                            np.min(data) - 0.25, np.max(data) + 0.25)

    def plot_voigt_and_linear_components(self, ax, energy, data, out, fix_param, dict_fit):
        """
        Plot Voigt and Linear components.

        Args:
            ax (Axes): Matplotlib Axes object.
            energy (np.array): Array of photon energy values.
            data (np.array): Array of ion counts.
            out (ModelResult or Parameters): Model result or parameters from the fit.
            fix_param (bool): Whether fixed parameters were used in the fit.
            dict_fit (dict): Dictionary of fit parameters.
        """
        ax.plot(energy, data)

        counter = 0
        # If fixed params are used for the fit
        if not fix_param:
            comps = out.eval_components(x=energy)
            for i in dict_fit.values():
                counter += 1
                # Plot Voigt components
                ax.plot(energy, comps[i[0]], '--', label='Voigt component ' + str(counter))
            # Plot linear component
            ax.plot(energy, comps['lin_'], '--', label='Linear component')
        # If automatic guess for param is used for the fit
        else:
            comps = out.values
            for i in range(1, 8):
                ind = str(i)
                counter += 1
                ax.plot(energy, self.skewed_voigt(energy, 
                                                  amplitude=comps['a' + ind], 
                                                  center=comps['c' + ind], 
                                                  sigma=comps['sigma'], 
                                                  gamma=comps['gamma'], 
                                                  skew=comps['skew']), 
                        '--', 
                        label='Voigt component ' + str(counter))

        self.configure_plot(ax, 
                            'Photon energy (eV)', 
                            'Ion counts (arb. units)')

    def plot_residuals_and_confidence_interval(self, ax, energy, data, residuals_final, RMS_neg, RMS_pos):
        """
        Plot residuals and confidence interval.

        Args:
            ax (Axes): Matplotlib Axes object.
            energy (np.array): Array of photon energy values.
            data (np.array): Array of ion counts.
            residuals_final (np.array): Residuals between data and best fit.
            RMS_neg (float): RMS value for initial guess residuals.
            RMS_pos (float): RMS value for final residuals.
        """
        ax.plot(energy, residuals_final, label='data-best fit')
        ax.fill_between(energy, np.sqrt(data) * -3.0, np.sqrt(data) * 3.0, color='gray', alpha=0.4, label='sqrt(data)*±3')
        self.configure_plot(ax, 'Photon energy (eV)', 'Residual counts (arb. units)')
        self.add_rms_text(ax, RMS_neg, RMS_pos)

    @staticmethod
    def configure_plot(ax, xlabel, ylabel, ymin=None, ymax=None):
        """
        Configure plot aesthetics.

        Args:
            ax (Axes): Matplotlib Axes object.
            xlabel (str): Label for the x-axis.
            ylabel (str): Label for the y-axis.
            ymin (float, optional): Minimum value for the y-axis. Defaults to None.
            ymax (float, optional): Maximum value for the y-axis. Defaults to None.
        """
        ax.legend()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if ymin is not None and ymax is not None:
            ax.set_ylim(ymin, ymax)
        ax.tick_params(axis='both', direction='in', length=6, top=True, right=True)
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(100000))
        ax.tick_params(which='minor', axis='both', direction='in', length=3, top=True, right=True)

    @staticmethod
    def add_rms_text(ax, rms_initial, rms_final):
        """
        Add RMS text to the plot.

        Args:
            ax (Axes): Matplotlib Axes object.
            rms_initial (float): RMS value for initial guess residuals.
            rms_final (float): RMS value for final residuals.
        """
        xmin_ax2, xmax_ax2 = ax.get_xlim()
        ymin_ax2, ymax_ax2 = ax.get_ylim()
        ax.text(((xmax_ax2 - xmin_ax2) * 0.55) + xmin_ax2, ((ymax_ax2 - ymin_ax2) * 0.2) + ymin_ax2, f'RMS initial = {rms_initial:.4f}')
        ax.text(((xmax_ax2 - xmin_ax2) * 0.55) + xmin_ax2, ((ymax_ax2 - ymin_ax2) * 0.15) + ymin_ax2, f'RMS final = {rms_final:.4f}')

    def plot_data_and_initial_guess(self, x, data, initial_guess):
        """
        Plot data and initial guess.

        Args:
            x (np.array): Array of x values.
            data (np.array): Array of y values.
            initial_guess (np.array): Array of initial guess values.
        """
        plt.rc("font", size=12, family='serif')
        fig, axes = plt.subplots(1, 1, figsize=(8.0, 16.0))
        axes.plot(x, initial_guess, 'orange', label='initial guess')
        axes.scatter(x, data, label='data')
        plt.show()
