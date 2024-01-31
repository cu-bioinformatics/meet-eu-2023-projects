import pandas as pd

# Replace 'file1.csv' and 'file2.csv' with your actual file names.
file1 = 'top1000_SPH2_SMILES.csv'  # The first CSV file with columns ID, SMILES.
file2 = 'grouped_data_sph2.csv'  # The second CSV file with columns proteinID, ligandID, mean, median, max, min, var.

# Read the CSV files
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Merge the DataFrames.
# df1's 'ID' column matches df2's 'ligandID' column.
merged_df = pd.merge(df1, df2, left_on='ID', right_on='ligandID')

sorted_df = merged_df.sort_values(by='mean', ascending=True)
sorted_df.drop('ID', axis=1, inplace=True)

# Save the merged DataFrame to a new CSV file
sorted_df.to_csv('results_sph2.csv', index=False)
