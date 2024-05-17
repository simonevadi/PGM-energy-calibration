import os
import pandas as pd
from nitrogen_fit import N2_fit

# Initialize N2_fit instance
n = N2_fit(None)

# cff values
cff = [1.6, 1.8, 2.0, 2.25, 3.5, 5.0, 8.0, 10.0, 15.0]

# Initialize an empty DataFrame
results_df = pd.DataFrame(columns=['cff', 'peak_energy', 'RP', 'vp_ratio'])

# Create the Results directory if it doesn't exist
os.makedirs('Results/ExampleN2Fit', exist_ok=True)

# Process each cff value
for c in cff:
    data_path = os.path.join('N2_data', f'N2_cff{c}.dat')
    save_path = os.path.join('Results', 'ExampleN2Fit', f'N2_cff{c}.png')
    
    if not os.path.exists(data_path):
        print(f'File not found: {data_path}')
        continue
    
    print(f'Fitting: {data_path}')
    
    try:
        peak_energy, RP, vp_ratio = n.fit_n2(data_path, print_fit_results=True, show=False, save_img=save_path)
    except Exception as e:
        print(f'Error fitting data for cff={c}: {e}')
        continue

    # Create a new DataFrame for the current results
    new_row = pd.DataFrame({'cff': [c], 'peak_energy': [peak_energy], 'RP': [RP], 'vp_ratio': [vp_ratio]})
    
    # Filter out empty or all-NA rows from new_row
    new_row_filtered = new_row.dropna(how='all')

    # Concatenate the new row to the results DataFrame
    results_df = pd.concat([results_df, new_row_filtered], ignore_index=True)

# Save the DataFrame to a CSV file
results_csv_path = os.path.join('Results', 'ExampleN2Fit', 'N2_fit_results.csv')
results_df.to_csv(results_csv_path, index=False)
