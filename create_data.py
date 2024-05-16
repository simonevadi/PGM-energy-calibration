import pandas as pd
import numpy as np

def apply_energy_shifts(input_file, cff, energy_shift):
    """
    Reads the input file, applies the specified energy shifts, and saves the modified data to new files.

    Args:
        input_file (str): Path to the input file.
        cff (list): List of cff values.
        energy_shift (list): List of energy shift values corresponding to each cff value.

    Returns:
        None
    """
    # Read the file with pandas, skipping the first row
    df = pd.read_csv(input_file, delim_whitespace=True, header=None, names=['energy', 'intensity'], skiprows=1)
    
    # Ensure cff and energy_shift lists are of the same length
    if len(cff) != len(energy_shift):
        raise ValueError("The 'cff' and 'energy_shift' lists must be of the same length.")

    # Loop through cff and energy_shift lists
    for c, shift in zip(cff, energy_shift):
        # Create a copy of the dataframe and apply the energy shift
        df_shifted = df.copy()
        df_shifted['energy'] += shift

        # Save the modified dataframe to a new file
        output_file = f'N2_data/N2_cff{c}.dat'
        df_shifted.to_csv(output_file, sep=' ', index=False, header=False)

        print(f"Saved: {output_file}")

# Example usage
cff          = [1.6, 1.8, 2.0, 2.25, 3.5, 5.0, 8.0, 10.0, 15.0]  # Example cff values
energy_shift = [ 0., -0.2, -0.4, -0.6, -1.2, -1.54, -1.9, -2., -2.13] # Corresponding energy shift values

apply_energy_shifts('N2_data/N2_PO4.dat', cff, energy_shift)
