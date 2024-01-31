#!/bin/bash

# Input and output directories
input_dir="output_SPH2_sdf"
output_dir="output_SPH2_pdbqt"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Count the number of files processed
count=0

# Convert .sdf files to .pdbqt
for input_file in "$input_dir"/*.sdf; do
    # Check if the file exists
    if [ ! -f "$input_file" ]; then
        echo "No .sdf files found in $input_dir."
        break
    fi

    # Get the base name of the file without extension
    base_name=$(basename -- "$input_file")
    file_name="${base_name%.*}"

    # Define output file path
    output_file="$output_dir/$file_name.pdbqt"

    # Perform the conversion using OpenBabel
    if obabel -isdf "$input_file" -opdbqt -O "$output_file"; then
        echo "Conversion of $base_name successful."
        ((count++))
    else
        echo "Error converting $base_name. Check input file and try again."
    fi
done

if [ $count -gt 0 ]; then
    echo "Batch conversion complete. $count files converted."
else
    echo "No files were converted."
fi
