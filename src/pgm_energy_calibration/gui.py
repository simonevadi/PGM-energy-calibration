import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import pandas as pd
import os
import json
import datetime

from pgm_energy_calibration.PGM_calibration import PGMCalibration

CONFIG_FILE = 'config/config.json'

class PlotApp:
    def __init__(self, root):
        """
        Initialize the PlotApp class.

        Args:
            root (tk.Tk): The root window of the Tkinter application.
        """
        # Store grating value
        self.grating = None
        self.dataframe = None
        self.mono = 'JenOptik'
        self.mono_types = ['JenOptik', 'PM4', 'Andreas']
        self.pgm = PGMCalibration(0) # Initialized with a random value

        # Load last file path from config file
        self.last_file_path = self.load_last_file_path()
        self.last_save_directory = self.load_last_save_directory()

        self.create_root_frame(root, "PGM Energy Calibration")
        plot_frame = self.create_plot_frame(root)
        control_frame = self.create_control_frame(root)

        self.button_font_style = ("Helvetica", 20, 'bold')
        self.text_font_style = ("Helvetica", 20)

        self.create_empty_plot(plot_frame)
        self.create_load_file_frame(control_frame, 1)
        self.create_set_mono_type(control_frame, 2)
        self.create_set_grating_frame(control_frame, 3)
        self.create_fit_frame(control_frame, 4)
        self.create_fit_results_frame(control_frame, 5)
        self.create_set_saving_param_frame(control_frame, 6)
        self.create_save_plot_frame(control_frame, 7)

    def create_root_frame(self, root, title):
        """
        Create the root frame for the application.

        Args:
            root (tk.Tk): The root window.
            title (str): The title of the root window.
        """
        root.title(title)

        # Set initial size of the window
        # root.geometry("1200x800")  # width x height

        # Make the root window scalable
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

    def create_plot_frame(self, root):
        """
        Create the frame for the plot.

        Args:
            root (tk.Tk): The root window.

        Returns:
            ttk.Frame: The created plot frame.
        """
        # Create a frame for the plot
        plot_frame = ttk.Frame(root)
        plot_frame.grid(row=0, column=0, sticky='nsew')
        return plot_frame

    def create_control_frame(self, root):
        """
        Create the frame for the control widgets.

        Args:
            root (tk.Tk): The root window.
        """
        # Create a frame for the controls
        control_frame = ttk.Frame(root)
        control_frame.grid(row=1, column=0, pady=10, sticky='ew')
        control_frame.columnconfigure(0, weight=1)
        return control_frame

    def create_load_file_frame(self, frame, row):
        """
        Create the frame for loading files.

        Args:
            frame (tk.Frame): The parent frame.
            row (int): The row position of the frame.
        """
        # Create the first row frame for Load File button, entry field, and browse button
        load_frame = ttk.Frame(frame)
        load_frame.grid(row=row, column=0, pady=5, sticky='ew')

        # Button: Browse
        browse_button = tk.Button(load_frame, text="Browse", command=self.browse_file, font=self.button_font_style)
        browse_button.pack(side=tk.LEFT, padx=10, pady=5)

        # Button: Load File
        load_button = tk.Button(load_frame, text="Load File", command=self.load_file_button, font=self.button_font_style)
        load_button.pack(side=tk.LEFT, padx=10, pady=5)

        # Create an entry widget for the file path with a default value
        self.entry = tk.Entry(load_frame, width=50, font=self.text_font_style)
        self.entry.insert(0, self.last_file_path or 'gui_data/ue112_2013.csv')
        self.entry.pack(side=tk.LEFT, padx=10, pady=5, fill=tk.X, expand=True)
        self.entry.bind('<Return>', lambda event: load_button.invoke())  # Bind Enter key to Load button

        # Create a label to display messages
        self.load_label = tk.Label(load_frame, text="File not set", font=self.text_font_style)
        self.load_label.pack(side=tk.LEFT, padx=10, pady=5)

    def create_set_grating_frame(self, frame, row):
        """
        Create the frame for setting the grating.

        Args:
            frame (tk.Frame): The parent frame.
            row (int): The row position of the frame.
        """
        # Create the second row frame for Set Grating button and entry field
        grating_frame = ttk.Frame(frame)
        grating_frame.grid(row=row, column=0, pady=5, sticky='ew')

        # Button: Set Grating
        set_grating_button = tk.Button(grating_frame, text="Set Grating", command=self.set_grating, font=self.button_font_style)
        set_grating_button.pack(side=tk.LEFT, padx=10, pady=5)

        # Label: Grating
        grating_label = tk.Label(grating_frame, text="Set Grating [l/mm]:", font=self.text_font_style)
        grating_label.pack(side=tk.LEFT, padx=10, pady=5)

        # Entry: Grating
        self.grating_entry = tk.Entry(grating_frame, width=10, font=self.text_font_style)
        self.grating_entry.insert(0, '600')
        self.grating_entry.pack(side=tk.LEFT, padx=10, pady=5)
        self.grating_entry.bind('<Return>', lambda event: set_grating_button.invoke())  # Bind Enter key to Set Grating button

        # Create a label to display messages
        self.label_grating = tk.Label(grating_frame, text="Grating not Set", font=self.text_font_style)
        self.label_grating.pack(side=tk.LEFT, padx=10, pady=5)

    def create_set_mono_type(self, frame, row):
        """
        Create the frame for setting the monochromator type.

        Args:
            frame (tk.Frame): The parent frame.
            row (int): The row position of the frame.
        """
        # Create the second row frame for Set mono button and entry field
        mono_frame = ttk.Frame(frame)
        mono_frame.grid(row=row, column=0, pady=5, sticky='ew')

        # Label: mono
        mono_label = tk.Label(mono_frame, text="Set Monochromator Type:", font=self.text_font_style)
        mono_label.pack(side=tk.LEFT, padx=10, pady=5)

        # Dropdown menu: mono types
        self.mono_var = tk.StringVar()
        self.mono_var.set(self.mono_types[0])  # Default value
        mono_options = self.mono_types
        self.mono_menu = tk.OptionMenu(mono_frame, self.mono_var, *mono_options, command=self.set_mono_type)
        self.mono_menu.config(font=self.text_font_style)
        self.mono_menu.pack(side=tk.LEFT, padx=10, pady=5)

        # Create a label to display messages
        self.label_mono = tk.Label(mono_frame, text=f"Monochromator type: {self.mono}", font=self.text_font_style)
        self.label_mono.pack(side=tk.LEFT, padx=10, pady=5)

    def create_set_saving_param_frame(self, frame, row):
        """
        Create the frame for setting saving parameters.

        Args:
            frame (tk.Frame): The parent frame.
            row (int): The row position of the frame.
        """
        # Create the second row frame for Set Grating button and entry field
        saving_p_frame = ttk.Frame(frame)
        saving_p_frame.grid(row=row, column=0, pady=5, sticky='ew')

        # Button: Set Beamline Name
        set_saving_p_button = tk.Button(saving_p_frame, text="Set Beamline Name", command=self.set_beamline_name, font=self.button_font_style)
        set_saving_p_button.pack(side=tk.LEFT, padx=10, pady=5)

        # Entry: Beamline Name
        self.saving_p_entry = tk.Entry(saving_p_frame, width=10, font=self.text_font_style)
        self.saving_p_entry.insert(0, 'UE112')
        self.saving_p_entry.pack(side=tk.LEFT, padx=10, pady=5)
        self.saving_p_entry.bind('<Return>', lambda event: set_saving_p_button.invoke())  # Bind Enter key to Set Grating button

        # Create a label to display messages
        self.saving_p_label = tk.Label(saving_p_frame, text="Beamline Name not Set", font=self.text_font_style)
        self.saving_p_label.pack(side=tk.LEFT, padx=10, pady=5)

    def create_fit_frame(self, frame, row):
        """
        Create the frame for the fit button and message label.

        Args:
            frame (tk.Frame): The parent frame.
            row (int): The row position of the frame.
        """
        # Create the third row frame for Fit button and message label
        fit_frame = ttk.Frame(frame)
        fit_frame.grid(row=row, column=0, pady=5, sticky='ew')

        # Create a button to fit the data
        fit_button = tk.Button(fit_frame, text="Fit", command=self.on_fit_button_click, font=self.button_font_style)
        fit_button.pack(side=tk.LEFT, padx=10, pady=5)

        # Create a label to display messages
        self.label_fit = tk.Label(fit_frame, text="", font=self.text_font_style)
        self.label_fit.pack(side=tk.LEFT, padx=10, pady=5)

    def create_fit_results_frame(self, frame, row):
        """
        Create the frame for displaying fit parameters.

        Args:
            frame (tk.Frame): The parent frame.
            row (int): The row position of the frame.
        """
        # Create the fourth row frame for displaying fit parameters
        results_frame = ttk.Frame(frame)
        results_frame.grid(row=row, column=0, pady=5, sticky='ew')

        # Create labels for displaying fit parameters
        self.dtheta_label = tk.Label(results_frame, text="DTheta: ", font=self.text_font_style)
        self.dtheta_label.pack(side=tk.LEFT, padx=10, pady=5)

        self.dbeta_label = tk.Label(results_frame, text="DBeta: ", font=self.text_font_style)
        self.dbeta_label.pack(side=tk.LEFT, padx=10, pady=5)

        self.e_opt_label = tk.Label(results_frame, text="E_opt: ", font=self.text_font_style)
        self.e_opt_label.pack(side=tk.LEFT, padx=10, pady=5)

    def create_empty_plot(self, plot_frame):
        """
        Create an empty plot.

        Args:
            plot_frame (tk.Frame): The frame where the plot will be displayed.
        """
        # Create an empty plot
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlabel('CFF Values')
        self.ax.set_ylabel('Measured Energies')
        self.ax.set_title('CFF Values vs Measured Energies')
        self.ax.set_xscale('log')
        self.ax.legend()
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add the Matplotlib navigation toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        toolbar.update()
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)

        # Make the plot frame expandable
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)

    def create_save_plot_frame(self, frame, row):
        """
        Create the frame for saving the plot.

        Args:
            frame (tk.Frame): The parent frame.
            row (int): The row position of the frame.
        """
        save_frame = ttk.Frame(frame)
        save_frame.grid(row=row, column=0, pady=5, sticky='ew')

        browse_save_button = tk.Button(save_frame, text="Browse", command=self.set_save_folder, font=self.button_font_style)
        browse_save_button.pack(side=tk.LEFT, padx=10, pady=5)

        entry_field = f"{self.last_save_directory}" if self.last_save_directory else "Folder not set"
        self.file_name_entry = tk.Entry(save_frame, width=50, font=self.text_font_style)
        self.file_name_entry.insert(0, entry_field)
        self.file_name_entry.pack(side=tk.LEFT, padx=10, pady=5, fill=tk.X, expand=True)

        save_button = tk.Button(save_frame, text="Save", command=self.save_results, font=self.button_font_style)
        save_button.pack(side=tk.LEFT, padx=10, pady=5)

        last_save_text = "File not Saved"
        self.save_label = tk.Label(save_frame, text=last_save_text, font=self.text_font_style)
        self.save_label.pack(side=tk.LEFT, padx=10, pady=5)

    def set_save_folder(self):
        """
        Set the folder for saving the plot.
        """
        initial_dir = self.load_last_save_directory()
        folder_path = filedialog.askdirectory(initialdir=initial_dir, title="Select folder")
        if folder_path:
            self.save_folder_path = folder_path
            self.save_label.config(text=f"Folder Set")
            self.save_last_save_directory(folder_path)

    def save_results(self):
        """
        Save the fit results and plot to files.
        """
        if self.save_folder_path:
            self.file_name_entry.delete(0, tk.END)
            self.file_name_entry.insert(0, os.path.join(self.save_folder_path))
            # Adding the new columns to the DataFrame
            fit_results_df = pd.DataFrame({
                'dtheta': [self.dtheta],
                'dbeta': [self.dbeta],
                'E_opt': [self.E_opt],
                'cost': [self.cost],
                'opt': [self.opt],
                'nfev': [self.nfev]
            })
            current_time = datetime.datetime.now().strftime("%Y_%m_%d_%H%M")
            data_name = os.path.join(self.save_folder_path, f'{current_time}_{self.beamline_name}')
            self.dataframe.to_csv(data_name+'_input_data.csv', index=False)
            fit_results_df.to_csv(data_name+'_fit_results.csv', index=False)
            self.fig.savefig(data_name+'.pdf')

            self.save_label.config(text=f"Results Saved")

    def set_mono_type(self, selected_type):
        """
        Set the monochromator type from the dropdown menu.

        Args:
            selected_type (str): The selected monochromator type.
        """
        self.mono = selected_type
        self.label_mono.config(text=f"Monochromator type: {self.mono}")

    def set_beamline_name(self):
        """
        Set the beamline name from the entry widget.
        """
        self.beamline_name = self.saving_p_entry.get()
        self.saving_p_label.config(text=f"Name set: {self.beamline_name}")

    def set_grating(self):
        """
        Set the grating value from the entry widget.
        """
        try:
            self.grating = float(self.grating_entry.get())
            self.label_grating.config(text=f"Grating set to: {self.grating} l/mm")
            self.pgm = PGMCalibration(self.grating)
        except ValueError:
            self.label_grating.config(text="Error: Please enter a valid number for grating")

    def load_file_button(self):
        """
        Load the file specified in the entry widget.
        """
        file_path = self.entry.get()
        if os.path.isfile(file_path):
            self.load_label.config(text=f"File loaded")
            self.plot_cff_vs_en(csv_file=file_path)
            self.save_last_file_path(file_path)
        else:
            self.load_label.config(text="Error: File not found")

    def browse_file(self):
        """
        Open a file dialog to browse and select a file to load.
        """
        initial_dir = os.path.dirname(self.last_file_path) if self.last_file_path else '.'
        file_path = filedialog.askopenfilename(initialdir=initial_dir, title="Select file",
                                               filetypes=(("CSV files", "*.csv"), ("all files", "*.*")))
        if file_path:
            self.entry.delete(0, tk.END)
            self.entry.insert(0, file_path)
            self.save_last_file_path(file_path)
            self.load_file_button()

    def save_last_save_directory(self, directory):
        """
        Save the last directory used for saving results to a configuration file.

        Args:
            directory (str): The directory path to save.
        """
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as config_file:
                config = json.load(config_file)
        else:
            config = {}

        config['last_save_directory'] = directory

        with open(CONFIG_FILE, 'w') as config_file:
            json.dump(config, config_file)

    def load_last_save_directory(self):
        """
        Load the last directory used for saving results from a configuration file.

        Returns:
            str: The last directory path, or the current directory if not found.
        """
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as config_file:
                config = json.load(config_file)
                return config.get('last_save_directory', '.')
        return '.'

    def plot_cff_vs_en(self, csv_file='gui_data/ue112_2013.csv'):
        """
        Plot CFF values vs. measured energies from a CSV file.

        Args:
            csv_file (str): The path to the CSV file.
        """
        # Clear previous plot
        self.ax.clear()

        # Read the CSV file
        df = pd.read_csv(csv_file)
        self.dataframe = df  # Store the dataframe for later use

        # Plot the data with different colors for different orders
        orders = df['orders'].unique()
        colors = plt.cm.viridis(np.linspace(0, 1, len(orders)))

        for order, color in zip(orders, colors):
            subset = df[df['orders'] == order]
            self.ax.scatter(subset['cff_values'], subset['measured_energies'], label=f'Order {order}', color=color)

        self.ax.set_xlabel('CFF Values')
        self.ax.set_ylabel('Measured Energies')
        self.ax.set_title('CFF Values vs Measured Energies')
        self.ax.set_xscale('log')
        self.ax.legend()
        self.ax.grid(True)

        # fix limits and labels
        cff_min = np.min(df['cff_values'])
        cff_max = np.max(df['cff_values'])
        en_min = np.min(df['measured_energies'])
        en_max = np.max(df['measured_energies'])

        xticks_positions, xticks_labels = self.pgm.generate_x_ticks_pos_and_label(cff_min, cff_max)

        self.ax.set_xticks(xticks_positions, labels=xticks_labels)

        self.ax.set_xlim((cff_min - .005, cff_max + .2))
        self.ax.set_ylim((en_min - 1, en_max + 1))
        self.canvas.draw()

    def plot_fit(self, measured_energies, cff_values, orders, DTheta=0, DBeta=0, E_opt=0):
        """
        Plot the fit results.

        Args:
            measured_energies (np.array): The measured energies.
            cff_values (np.array): The CFF values.
            orders (np.array): The orders.
            DTheta (float): The delta theta value.
            DBeta (float): The delta beta value.
            E_opt (float): The optimized energy.
        """
        unique_orders = np.unique(orders)
        colors = ['k', 'royalblue', 'darkgoldenrod', 'g', 'magenta', 'orange', 'red']

        for idx, order in enumerate(unique_orders):
            # Filter data for the current order
            mask = (orders == order)
            cff_order = cff_values[mask]
            measured_order = measured_energies[mask]

            # Use more points to plot the fit results
            cff_plot = np.arange(np.min(cff_values), np.max(cff_values), .01)
            mask = (cff_plot < 0.9) | (cff_plot > 1.1)
            cff_plot = cff_plot[mask]

            mask = cff_plot < 1
            orders_plot = np.full(cff_plot.shape, order)
            orders_plot[mask] *= -1

            orders_ig = np.full(cff_order.shape, order)
            mask = cff_order < 1
            orders_ig[mask] *= -1

            fitted_energies = self.pgm.shifted_energy(cff_plot, orders_plot, np.deg2rad(DTheta), np.deg2rad(DBeta), E_opt)

            self.ax.plot(cff_plot, fitted_energies, '-', color=colors[idx], label=f'Fitted Order {order}', alpha=0.7)

        self.ax.legend(loc='upper right')
        self.canvas.draw()

    def on_fit_button_click(self):
        """
        Perform fitting on the loaded data and update the plot and results.
        """
        if self.dataframe is not None:
            print(f'Loaded data\n{self.dataframe}')
            self.load_file_button()
            # Instantiate the calibration class

            # Set a reasonable initial guess manually
            self.pgm.set_initial_guess(automatic_guess=True)  # Example initial guess values
            measured_energies = self.dataframe['measured_energies']
            cff_values = self.dataframe['cff_values']
            orders = self.dataframe['orders']
            # Perform fitting
            self.dtheta, self.dbeta, self.E_opt, \
                self.cost, self.opt, self.nfev = self.pgm.fit_parameters(measured_energies,
                                                                         cff_values, orders,
                                                                         return_fit_eval=True)

            self.pgm.print_fit_results()
            self.plot_fit(measured_energies, cff_values, orders,
                          DTheta=self.dtheta, DBeta=self.dbeta, E_opt=self.E_opt)

            # Update labels with fit parameters
            self.label_fit.config(text=f"Sum Squared Residuals: {np.round(self.cost, 5)}, Optimality: {np.round(self.opt, 5)}, evaluations {self.nfev}")
            if self.mono == 'JenOptik':
                dbeta = -self.dbeta
                dtheta = self.dtheta
                units = 'deg'
            elif self.mono == 'PM4':
                dbeta = self.dbeta
                dtheta = self.dtheta
                units = 'deg'
            elif self.mono == 'Andreas':
                dbeta = abs(self.dbeta) * 3600 / 0.0225
                dtheta = abs(self.dtheta) * 3600 / 0.0225
                units = ''
            self.dtheta_label.config(text=f"DTheta: {np.round(dtheta, 5)} {units}")
            self.dbeta_label.config(text=f"DBeta: {np.round(dbeta, 5)} {units}")
            self.e_opt_label.config(text=f"E_opt: {np.round(self.E_opt, 3)} eV")
        else:
            print("No data loaded.")

    def save_last_file_path(self, path):
        """
        Save the last file path to a configuration file.

        Args:
            path (str): The file path to save.
        """
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as config_file:
                config = json.load(config_file)
        else:
            config = {}

        config['last_file_path'] = path

        with open(CONFIG_FILE, 'w') as config_file:
            json.dump(config, config_file)

    def load_last_file_path(self):
        """
        Load the last file path from a configuration file.

        Returns:
            str: The last file path, or None if not found.
        """
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as config_file:
                config = json.load(config_file)
                return config.get('last_file_path', None)
        return None


def main():
    # Create the main window
    root = tk.Tk()
    app = PlotApp(root)

    # Run the application
    root.mainloop()


if __name__ == '__main__':
    main()
