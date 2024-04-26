#!/bin/bash

# This script recursively merges the contents of two archives into a single archive.
# If there are conflicting files, the script renames the conflicting files in the second archive.
# All files are extracted to temporary directories, merged, and then compressed into a new archive.
# Directory structure within the archives is preserved during the merge.

# The script requires the 7z utility to be installed on the system.
# Usage: ./zipper_merge.bash <first_archive> <second_archive> <merged_archive>


if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <first_archive> <second_archive> <merged_archive>"
    exit 1
fi

first_archive=$1
second_archive=$2
merged_archive=$3

# Function to find the maximum replicate number in files within a directory
find_max_replicate_number() {
    local directory=$1
    local max_R=0
    for file in "$directory"/*.csv; do
        # Extract the replicate number from the filename
        rep_no=$(basename "$file" | grep -oE 'R_[0-9]+' | cut -d'_' -f2)
        if [ "$rep_no" -gt "$max_R" ]; then
            max_R=$rep_no
        fi
    done
    echo "$max_R"
}

# Function to handle conflicts by renaming conflicting files in temp_dir2
handle_conflicts() {
    local temp_dir1=$1
    local temp_dir2=$2

    # Find all CSV files in temp_dir2 recursively
    while IFS= read -r -d '' conflicting_file; do
        # Extract the corresponding file path in temp_dir1
        corresponding_file="${conflicting_file/$temp_dir2/$temp_dir1}"
        parent_directory=$(dirname "$corresponding_file")

        # Check if the file also exists in temp_dir1
        if [ -f "$corresponding_file" ]; then
            max_R=$(find_max_replicate_number "$parent_directory")
            new_name="${conflicting_file%_R_*}_R_$((max_R + 1)).csv"
            mv "$conflicting_file" "$new_name"
        fi
    done < <(find "$temp_dir2" -type f -name '*.csv' -print0)
}

# Extract contents of the first archive
mkdir temp_dir1
7z x "$first_archive" -o./temp_dir1

# Extract contents of the second archive
mkdir temp_dir2
7z x "$second_archive" -o./temp_dir2

# Merge contents of temporary directories into a new directory
mkdir merged_dir
rsync -a ./temp_dir1/ merged_dir/

# Handle conflicts in temp_dir2 recursively
handle_conflicts "temp_dir1" "temp_dir2"

# Merge contents of temp_dir2 into merged_dir
rsync -a --ignore-existing ./temp_dir2/ merged_dir/

# Create a new archive from the merged directory
7z a "$merged_archive" merged_dir

# Clean up temporary directories
rm -r temp_dir1 temp_dir2 merged_dir