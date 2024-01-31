import csv
import sys

# Function to convert .cxsmiles to .csv
def convert_cxsmiles_to_csv(cxsmiles_file, csv_file):
    with open(cxsmiles_file, 'r') as infile, open(csv_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t') # Assuming tab-separated values in .cxsmiles
        writer = csv.writer(outfile)
        
        # Write the header to the CSV file
        writer.writerow(['smiles', 'idnumber', 'Type'])
        
        # Iterate over the lines in the .cxsmiles file and write to the .csv file
        for row in reader:
            if len(row) == 3: # Each row is expected to have 3 columns: smiles, idnumber, and Type
                writer.writerow(row)

# Main function to handle command line arguments
def main():
    if len(sys.argv) != 3:
        print("Usage: python convert_cxsmiles_to_csv.py [input.cxsmiles] [output.csv]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    convert_cxsmiles_to_csv(input_file, output_file)

if __name__ == "__main__":
    main()

