#!/bin/bash

# Input file (CSV format with columns: ID, SMILES)
input_file="top1000_SPH2_SMILES.csv"

# Output directory
output_dir="output_SPH2_sdf"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Count the number of files processed
count=0

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file $input_file not found."
    exit 1
fi

# Read CSV file and convert SMILES to 3D SDF
while IFS=, read -r id smiles; do
    # Get the base name of the ID without spaces
    id_no_spaces=$(echo "$id" | tr -d '[:space:]')

    # Define output file path
    output_file="$output_dir/$id_no_spaces.sdf"

    # Perform the conversion using OpenBabel with 3D generation
    if obabel -:"$smiles" --gen3d -osdf -O "$output_file"; then
        echo "Conversion of $id to 3D SDF successful."
        ((count++))
    else
        echo "Error converting $id to 3D SDF. Check input data and try again."
    fi
done < "$input_file"

if [ $count -gt 0 ]; then
    echo "Batch conversion complete. $count 3D SDF files converted."
else
    echo "No 3D SDF files were converted."
fi
