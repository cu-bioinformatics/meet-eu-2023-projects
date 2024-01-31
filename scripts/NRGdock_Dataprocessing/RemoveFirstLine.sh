#!/bin/bash

# Check if a file name is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
fi

# Assign the file name to a variable
FILE=$1

# Check if the file exists
if [ ! -f "$FILE" ]; then
    echo "Error: File not found!"
    exit 1
fi

# Use tail to skip the first line and save to a temporary file
tail -n +2 "$FILE" > "$FILE.tmp"

# Move the temporary file to original file name
mv "$FILE.tmp" "$FILE"

echo "First line removed from $FILE"

