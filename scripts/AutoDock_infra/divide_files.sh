#!/bin/bash

if [ $# -ne 2 ]; then # argument check
    echo "Usage: $0 [source directory] [number of groups]"
    exit 1
fi

# Directory containing files
source_directory=$1
# Number of groups
num_groups=$2

# Create an array of all files in the source directory
cd "$source_directory"
files=(*)

# Calculate number of files per group
num_files=${#files[@]}
files_per_group=$(( (num_files + num_groups - 1) / num_groups ))

# Create groups and distribute files
for ((i=1; i<=num_groups; i++)); do
    group_dir="${i}"
    mkdir -p "$group_dir"
    # Calculate array slice indices
    start=$(( (i - 1) * files_per_group ))
    end=$(( start + files_per_group ))
    # Copy files to group directory
    echo "$group_dir"
    for file in "${files[@]:$start:$files_per_group}"; do
        # mv "$source_directory/$file" "$group_dir"
        cp "$file" "$group_dir"
    done
done
