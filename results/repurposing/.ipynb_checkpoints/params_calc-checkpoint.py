import pandas as pd

# Load the dataset
df = pd.read_csv('conformations_merged_dataset.csv')

# Calculate Variance of Affinity Scores and Rankings
mean_columns = [col for col in df.columns if 'mean' in col]
rank_columns = [col for col in df.columns if 'rank' in col]

df['Variance of Affinity Scores'] = df[mean_columns].var(axis=1)
df['Variance of Rankings'] = df[rank_columns].var(axis=1)

# Calculate Mean Ranking and Mean Affinity Score
df['Mean Ranking'] = df[rank_columns].mean(axis=1)
df['Mean Affinity Score'] = df[mean_columns].mean(axis=1)

# Find IDs with Max/Min Affinity and Max/Min Rank
df['Max Affinity ID'] = df[mean_columns].idxmax(axis=1).str.replace(' - mean', '')
df['Min Affinity ID'] = df[mean_columns].idxmin(axis=1).str.replace(' - mean', '')
df['Max Rank ID'] = df[rank_columns].idxmax(axis=1).str.replace(' - rank', '')
df['Min Rank ID'] = df[rank_columns].idxmin(axis=1).str.replace(' - rank', '')

# Save the updated dataframe
df.to_csv('updated_dataset.csv', index=False)

print("Dataset updated with new columns. Saved as 'updated_dataset.csv'")
