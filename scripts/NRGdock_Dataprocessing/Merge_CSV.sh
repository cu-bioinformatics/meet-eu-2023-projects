#!/bin/bash

# Define file names
file1="top5000_SPH1.csv"
file2="Diverse_REAL.csv"
output_file="$(basename "$file1" .csv)_SMILES.csv"

# Check if files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
    echo "Error: File(s) not found."
    exit 1
fi

# Get the number of CPU cores for parallel processing
num_cores=$(sysctl -n hw.ncpu)

# Store headers in variables
header1=$(head -n 1 $file1)
header2=$(head -n 1 $file2)

# Sort the second file by the idnumber column (2nd column) with parallel processing, excluding header
sorted_second_file="sorted_$(basename "$file2")"
tail -n +2 $file2 | sort -t, -k2,2 --parallel=$num_cores > "$sorted_second_file"

# Check if sorted file is empty
if [ ! -s "$sorted_second_file" ]; then
    echo "Error: Sorted file is empty. Sorting failed."
    exit 1
fi

# Join the files excluding headers from the first file, then add the combined header
output_without_header="temp_output.csv"
tail -n +2 $file1 | join -t, -1 1 -2 2 -a 1 -o '1.1,1.2,2.1,2.3' - "$sorted_second_file" > "$output_without_header"

# Check if output file is empty
if [ ! -s "$output_without_header" ]; then
    echo "Error: Output file is empty. Joining failed."
    exit 1
fi

# Combine header and output
echo "$header1,smiles,Type" > "$output_file"
cat "$output_without_header" >> "$output_file"

# Clean up
rm "$sorted_second_file" "$output_without_header"

echo "Merge complete. Output saved to $output_file"

