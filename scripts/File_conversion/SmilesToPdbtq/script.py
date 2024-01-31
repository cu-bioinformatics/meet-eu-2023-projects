import argparse
import pandas as pd

def process_csv(input_file, output_file):
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Remove the 'ID' column containing the chemical codes
    df.drop('ID', axis=1, inplace=True)

    # Rename the first column (previously unnamed) to 'ID'
    df.rename(columns={df.columns[0]: 'ID'}, inplace=True)

    # Save the modified DataFrame to the output CSV file
    df.to_csv(output_file, index=False)

    print(f"CSV file has been modified and saved as '{output_file}'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a CSV file.")
    parser.add_argument('input_file', type=str, help='Path to the input CSV file')
    parser.add_argument('output_file', type=str, help='Path to the output CSV file')

    args = parser.parse_args()

    process_csv(args.input_file, args.output_file)
