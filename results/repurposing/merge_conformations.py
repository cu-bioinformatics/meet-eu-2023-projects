import pandas as pd

# File names of the CSV files
files = ['7RDX/grouped_data.csv', '7RDY/repurposing_results.csv', '7RDZ/grouped_data.csv', '7RE0/grouped_data.csv']

# Read the CSV files and process each one
dataframes = []
for file in files:
    df = pd.read_csv(file, header=0)

    # Rename the first unnamed column to '<ProteinID Prefix> - rank'
    protein_id_prefix = df['proteinID'].iloc[0][:4]  # Extract the first 4 letters of ProteinID
    df.rename(columns={df.columns[0]: protein_id_prefix + ' - rank',
                       'mean': protein_id_prefix + ' - mean'}, inplace=True)

    # Set the index to ligandID for merging
    df.set_index('ligandID', inplace=True)
    dataframes.append(df[[protein_id_prefix + ' - rank', protein_id_prefix + ' - mean']])

# Merge all dataframes based on LigandID using 'outer' join
merged_df = pd.concat(dataframes, axis=1, join='outer')

# Save the merged dataframe to a new CSV file
merged_df.to_csv('merged_dataset.csv')

print("Merging completed. Output saved to 'conformations_merged_dataset.csv'")

