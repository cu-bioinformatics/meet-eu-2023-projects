#!/bin/bash


# Input and output directories
input_dir="struc_files"
output_dir="output_pdbqt"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Convert .sdf files to .pdbqt
for input_file in "$input_dir"/*.sdf; do
    # Get the base name of the file without extension
    base_name=$(basename -- "$input_file")
    file_name="${base_name%.*}"

    # Define output file path
    output_file="$output_dir/$file_name.pdbqt"

    # Perform the conversion using OpenBabel
    touch "$output_file"
    obabel -isdf "$input_file" -opdbqt "$output_file"

    # Check if the conversion was successful
    if [ $? -eq 0 ]; then
        echo "Conversion of $base_name successful."
    else
        echo "Error converting $base_name. Check input file and try again."
    fi
done

echo "Batch conversion complete."
