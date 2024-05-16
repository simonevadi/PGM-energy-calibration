import os
import pandas as pd
from nitrogen_fit import N2_fit

# Initialize N2_fit instance
n = N2_fit(None)

# cff values
cff = [1.6, 1.8, 2.0, 2.25, 3.5, 5.0, 8.0, 10.0, 15.0]

# Initialize an empty DataFrame
results_df = pd.DataFrame(columns=['cff', 'peak_energy', 'RP', 'vp_ratio'])

# Process each cff value
for c in cff:
    data_path = os.path.join('N2_data', f'N2_cff{c}.dat')
    save_path = os.path.join('Results', f'N2_cff{c}.png')
    print(f'Fitting: {data_path}')
    peak_energy, RP, vp_ratio = n.fit_n2(data_path, show=False, save_img=save_path)

    # Create a new DataFrame for the current results
    new_row = pd.DataFrame({'cff': [c], 'peak_energy': [peak_energy], 'RP': [RP], 'vp_ratio': [vp_ratio]})

    # Concatenate the new row to the results DataFrame
    results_df = pd.concat([results_df, new_row], ignore_index=True)

# Save the DataFrame to a CSV file
results_csv_path = os.path.join('Results', 'N2_fit_results.csv')
results_df.to_csv(results_csv_path, index=False)
