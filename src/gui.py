import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
import os
import json

from PGM_calibration import PGMCalibration

CONFIG_FILE = 'config/config.json'

class PlotApp:
    def __init__(self, root):
        self.root = root
        self.root.title("PGM Energy calibration")

        # Load last file path from config file
        self.last_file_path = self.load_last_file_path()

        # Make the root window scalable
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Create a frame for the plot
        self.plot_frame = ttk.Frame(root)
        self.plot_frame.grid(row=0, column=0, sticky='nsew')

        # Create a frame for the controls
        self.control_frame = ttk.Frame(root)
        self.control_frame.grid(row=1, column=0, pady=10, sticky='ew')
        self.control_frame.columnconfigure(0, weight=1)

        # Define a font for control frame widgets
        button_font_style = ("Helvetica", 12, 'bold')
        text_font_style = ("Helvetica", 12)

        # Create the first row frame for Load File button, entry field, and browse button
        self.load_frame = ttk.Frame(self.control_frame)
        self.load_frame.grid(row=0, column=0, pady=5, sticky='ew')

        # Button: file path
        load_button = tk.Button(self.load_frame, text="Load File", command=self.on_button_click, font=button_font_style)
        load_button.pack(side=tk.LEFT, padx=10, pady=5)

        # Create an entry widget for the file path with a default value
        self.entry = tk.Entry(self.load_frame, width=50, font=text_font_style)
        self.entry.insert(0, self.last_file_path or 'gui_data/ue112_2013.csv')
        self.entry.pack(side=tk.LEFT, padx=10, pady=5, fill=tk.X, expand=True)

        # Button: Browse
        browse_button = tk.Button(self.load_frame, text="Browse", command=self.browse_file, font=button_font_style)
        browse_button.pack(side=tk.LEFT, padx=10, pady=5)

        # Create the second row frame for Fit button and message label
        self.fit_frame = ttk.Frame(self.control_frame)
        self.fit_frame.grid(row=1, column=0, pady=5, sticky='ew')

        # Create a button to fit the data
        fit_button = tk.Button(self.fit_frame, text="Fit", command=self.on_fit_button_click, font=button_font_style)
        fit_button.pack(side=tk.LEFT, padx=10, pady=5)

        # Create a label to display messages
        self.label = tk.Label(self.fit_frame, text="", font=text_font_style)
        self.label.pack(side=tk.LEFT, padx=10, pady=5)

        # Create the third row frame for displaying fit parameters
        self.results_frame = ttk.Frame(self.control_frame)
        self.results_frame.grid(row=2, column=0, pady=5, sticky='ew')

        # Create labels for displaying fit parameters
        self.dtheta_label = tk.Label(self.results_frame, text="DTheta: ", font=text_font_style)
        self.dtheta_label.pack(side=tk.LEFT, padx=10, pady=5)

        self.dbeta_label = tk.Label(self.results_frame, text="DBeta: ", font=text_font_style)
        self.dbeta_label.pack(side=tk.LEFT, padx=10, pady=5)

        self.e_opt_label = tk.Label(self.results_frame, text="E_opt: ", font=text_font_style)
        self.e_opt_label.pack(side=tk.LEFT, padx=10, pady=5)

        # Create an empty plot
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlabel('cff')
        self.ax.set_ylabel('Energiy [eV]')
        self.ax.set_title('Energy vs cff')
        self.ax.set_xscale('log')
        self.ax.legend()
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.dataframe = None

    def on_button_click(self):
        file_path = self.entry.get()
        if os.path.isfile(file_path):
            self.label.config(text=f"File found: {file_path}")
            self.plot_cff_vs_en(csv_file=file_path)
            self.save_last_file_path(file_path)
        else:
            self.label.config(text="Error: File not found")

    def browse_file(self):
        initial_dir = os.path.dirname(self.last_file_path) if self.last_file_path else '.'
        file_path = filedialog.askopenfilename(initialdir=initial_dir, title="Select file",
                                               filetypes=(("CSV files", "*.csv"), ("all files", "*.*")))
        if file_path:
            self.entry.delete(0, tk.END)
            self.entry.insert(0, file_path)
            self.save_last_file_path(file_path)

    def plot_cff_vs_en(self, csv_file='gui_data/ue112_2013.csv'):
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

        self.canvas.draw()

    def plot_fit(self, pgm_class, measured_energies, cff_values, orders, DTheta=0, DBeta=0, E_opt=0):
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

            fitted_energies = pgm_class.shifted_energy(cff_plot, orders_plot, np.deg2rad(DTheta), np.deg2rad(DBeta), E_opt)

            self.ax.plot(cff_plot, fitted_energies, '-', color=colors[idx], label=f'Fitted Order {order}', alpha=0.7)
        
        self.ax.legend()
        self.canvas.draw()

    def on_fit_button_click(self):
        if self.dataframe is not None:
            print(self.dataframe)
            # Instantiate the calibration class
            c = PGMCalibration(600)

            # Set a reasonable initial guess manually
            c.set_initial_guess(automatic_guess=True)  # Example initial guess values
            measured_energies = self.dataframe['measured_energies']
            cff_values        = self.dataframe['cff_values']
            orders            = self.dataframe['orders']
            # Perform fitting
            dtheta, dbeta, E_opt = c.fit_parameters(measured_energies, cff_values, orders)

            c.print_fit_results()
            self.plot_fit(c, measured_energies, cff_values, orders, DTheta=dtheta, DBeta=dbeta, E_opt=E_opt)

            # Update labels with fit parameters
            self.dtheta_label.config(text=f"DTheta: {np.round(dtheta, 5)} deg")
            self.dbeta_label.config(text=f"DBeta: {np.round(dbeta, 5)} deg")
            self.e_opt_label.config(text=f"E_opt: {np.round(E_opt, 3)} eV")
        else:
            print("No data loaded.")

    def save_last_file_path(self, path):
        with open(CONFIG_FILE, 'w') as config_file:
            json.dump({'last_file_path': path}, config_file)

    def load_last_file_path(self):
        if os.path.exists(CONFIG_FILE):
            with open(CONFIG_FILE, 'r') as config_file:
                config = json.load(config_file)
                return config.get('last_file_path', None)
        return None



# Create the main window
root = tk.Tk()
app = PlotApp(root)

# Run the application
root.mainloop()
