import dask.dataframe as dd
import sys

def merge_files(file1_path, file2_path, output_file_path):
    # Load both files with Dask DataFrame
    df1 = dd.read_csv(file1_path)
    df2 = dd.read_csv(file2_path)

    # Sorting the second file by 'idnumber' (if it's not already sorted)
    df2 = df2.sort_values('idnumber')

    # Merge the DataFrames on the ID columns
    merged_df = df1.merge(df2, left_on='ID', right_on='idnumber', how='left')

    # Fill missing values with a placeholder if needed
    merged_df = merged_df.fillna('NA')

    # Save the merged DataFrame to a new CSV file
    merged_df.to_csv(output_file_path, single_file=True, index=False)

    print(f"Merged file saved as {output_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_csv.py <file1.csv> <file2.csv> <output.csv>")
        sys.exit(1)

    file1_path, file2_path, output_file_path = sys.argv[1], sys.argv[2], sys.argv[3]
    merge_files(file1_path, file2_path, output_file_path)
