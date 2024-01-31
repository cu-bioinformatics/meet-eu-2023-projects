import os
import csv
import re 

# Directory containing the subdirectories
directory_path = 'user_predictions_small'

# Regular expression to extract the index number
index_pattern = re.compile(r'index(\d+)')

# Prepare a CSV file to write the data
with open('extracted_data.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Index', 'PDB Code', 'SMILES Chemical Structure'])

    # Iterate over each item in the directory
    for item in os.listdir(directory_path):
        if os.path.isdir(os.path.join(directory_path, item)):
            # Use regular expression to find the index number
            match = index_pattern.search(item)
            if match:
                index = match.group(1)  # Extracting just the numeric part of the index

                # Extracting the pdb code and SMILES chemical structure
                parts = item.split('_')
                pdb_code = parts[1]
                smiles = item.split('____')[1]

                # Write the extracted data to the CSV file
                writer.writerow([index, pdb_code, smiles])

print("Data extraction complete. Check the extracted_data.csv file.")

