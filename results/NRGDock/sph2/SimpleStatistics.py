import pandas as pd

# Step 1: Read the CSV file
# Replace 'your_data.csv' with the path to your CSV file
df = pd.read_csv('merged_results_sph2.csv')

# Step 2: Group the data by 'ligandID'
grouped = df.groupby('ligandID')

# Step 3: Calculate the statistics
stats = grouped['affinity'].agg(['mean', 'median', 'max', 'min', 'var']).reset_index()

# Adding 'proteinID' to the new DataFrame
# Assuming all entries within a group have the same 'proteinID'
stats['proteinID'] = grouped['proteinID'].agg('first').values

# Reordering the columns to match your specified format
stats = stats[['proteinID', 'ligandID', 'mean', 'median', 'max', 'min', 'var']]

# Sort the DataFrame by the mean affinity in descending order
sorted_stats = stats.sort_values(by='mean', ascending=False)

# Step 4: Save the new data to a CSV file
# Replace 'output_data.csv' with your desired output file name
sorted_stats.to_csv('grouped_data_sph2.csv', index=False)

