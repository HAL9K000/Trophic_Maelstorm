#!/bin/bash

# This script is used to reorganise directories using Python scripts reorganise_dir.py and reorganise_prelims_dir.py.
# Each of these scripts can accept necessary arguements as compiler flags using argparse, and these 
# input parameters are read from the input files ( paths to which are provided as $1 and $2 respectively).
# The script will run the Python scripts for each line in the input files, and compress the output directories
# It will also optionally run additional commands provided in an optional file (if provided as $5).
# The script will run for a specified number of iterations (specified by $4) with a time gap between each iteration (specified by $3).
# The script will log any errors to error_log.txt.
# Finally, this script will compress all found outdirs and rename any existing outdirs with a 
# timestamp to avoid overwrites of previously existing data.

# Main script
if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 <Path/to/dir/inputs/reorganise_dir.txt> <Path/to/dir/inputs/reorganise_prelims_dir.txt> <time_gap_in_hours> <stop_iter> [Optional: Path/to/optionalcommands.txt]"
    exit 1
fi

in_par_frame_file="$1"
in_par_prelims_file="$2"
time_gap_hours="$3"
stop_iter="$4"
commands_file="$5"         # Optional file with commands to execute (will be empty if not provided)
current_iter=0

# Convert time gap from hours to seconds
time_gap_seconds=$(( time_gap_hours * 3600 ))

# Function to log errors
log_error() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') ERROR: $1" >> error_log.txt
}

# Function to compress the output directory using 7z
compress_dir() {
    local outdir="$1"
    local timestamp=$(date +'%d%b%Y_%H%M')
    local archive_name="${outdir}_${timestamp}.7z"

    # Compress using 7z
    7z a "$archive_name" "$outdir" &>/dev/null
    if [[ $? -ne 0 ]]; then
        log_error "Failed to compress $outdir"
    else
        echo "Compressed $outdir to $archive_name"
    fi
}

# Function to execute additional commands (optional fourth argument)
execute_additional_commands() {
    local commands_file="$1"
    
    # Check if the file exists
    if [[ ! -f "$commands_file" ]]; then
        log_error "Commands file not found: $commands_file"
        return 1
    fi

    # Execute each command line by line
    while IFS= read -r command || [[ -n "$command" ]]; do
        echo "Executing command: $command"

        # Check if the command requires additional input (e.g., password/key)
        if [[ "$command" =~ ^(rsync|ssh) ]]; then
            echo "Command might require additional input..."

            # Try to read the next line for additional input (password or key)
            read -r additional_input || {
                log_error "No additional input found for command: $command"
                continue
            }

            echo "Providing additional input: $additional_input"

            # Use `echo` to simulate entering the input into the command
            echo "$additional_input" | eval "$command"
        else
            eval "$command"
        fi

        if [[ $? -ne 0 ]]; then
            log_error "Failed to execute command: $command"
        fi
    done < "$commands_file"
}

# Loop for the number of iterations specified by stop_iter
while [[ $current_iter -lt $stop_iter ]]; do
    # Read the params.txt file line by line
    declare -A renamed_outdirs  # Associative array to track renamed outdirs in this iteration
    declare -A found_outdirs    # Associative array to track all found outdirs in this iteration

    # Execute additional commands at the start of the outer loop (if the fourth argument is provided)
    if [[ -n "$commands_file" ]]; then
        echo "Executing additional commands from $commands_file..."
        execute_additional_commands "$commands_file"
    fi

    while IFS= read -r line || [[ -n "$line" ]]; do
        # Extract input parameters from the line (assuming well-formed lines)
        indir=$(echo "$line" | grep -oP '(?<=--indir )[^ ]+')
        outdir=$(echo "$line" | grep -oP '(?<=--outdir )[^ ]+')
        
        if [[ -z "$indir" || -z "$outdir" ]]; then
            log_error "Invalid line: $line"
            continue
        fi

        #Remove trailing whitespaces and new lines from outdir
        # Remove trailing whitespaces and backslash (/) from outdir and indir
        outdir=$(echo "$outdir" | sed 's/[[:space:]]*\/$//' | sed 's/[[:space:]]*$//')
        indir=$(echo "$outdir" | sed 's/[[:space:]]*\/$//' | sed 's/[[:space:]]*$//')

        # Check if outdir exists and has non-zero size
        if [[ -d "$outdir" && $(du -sb "$outdir" | awk '{print $1}') -gt 0 ]]; then
            if [[ -z "${renamed_outdirs[$outdir]}" ]]; then
                #Only rename outdir once in this iteration (i.e, iff not already renamed in this iteration and it exists)
                # If so, rename the directory with the timestamp
                timestamp=$(date +'%H%M')
                mv "$outdir" "${outdir}_${timestamp}"
                if [[ $? -ne 0 ]]; then
                    log_error "Failed to rename $outdir"
                    continue
                fi
                echo "Renamed existing $outdir to ${outdir}_${timestamp}"
                log_error "Renamed existing $outdir to ${outdir}_${timestamp}"
                # Adding outdir to renamed_outdirs dictionary.
                renamed_outdirs["${outdir}"]="${outdir}_${timestamp}"
            fi
                
        fi

        # Run the Python script
        python3.11 reorganise_dir.py $line
        if [[ $? -ne 0 ]]; then
            log_error "Python script reorganise_dir.py failed for line: $line"
            continue
        fi

        # Add outdir to found_outdirs dictionary only if it is not already present
        if [[ -z "${found_outdirs[$outdir]}" ]]; then
            found_outdirs["$outdir"]=1
        fi

    done < "$in_par_frame_file"

    while IFS= read -r line || [[ -n "$line" ]]; do
        # Extract input parameters from the line (assuming well-formed lines)
        indir=$(echo "$line" | grep -oP '(?<=--indir )[^ ]+')
        outdir=$(echo "$line" | grep -oP '(?<=--outdir )[^ ]+')
        
        if [[ -z "$indir" || -z "$outdir" ]]; then
            log_error "Invalid line: $line"
            continue
        fi

        # Remove trailing whitespaces and backslash (/) from outdir and indir
        outdir=$(echo "$outdir" | sed 's/[[:space:]]*\/$//' | sed 's/[[:space:]]*$//')
        indir=$(echo "$outdir" | sed 's/[[:space:]]*\/$//' | sed 's/[[:space:]]*$//')

        # Check if outdir exists and has non-zero size
        if [[ -d "$outdir" && $(du -sb "$outdir" | awk '{print $1}') -gt 0 ]]; then
            if [[ -z "${renamed_outdirs[$outdir]}" ]]; then
                #Only rename outdir once in this iteration (i.e, iff not already renamed in this iteration and it exists)
                # If so, rename the directory with the timestamp
                timestamp=$(date +'%H%M')
                mv "$outdir" "${outdir}_${timestamp}"
                if [[ $? -ne 0 ]]; then
                    log_error "Failed to rename $outdir"
                    continue
                fi
                echo "Renamed existing $outdir to ${outdir}_${timestamp}"
                log_error "Renamed existing $outdir to ${outdir}_${timestamp}"
                # Adding outdir to renamed_outdirs dictionary.
                renamed_outdirs["${outdir}"]="${outdir}_${timestamp}"
            fi
                
        fi

        # Run the Python script
        python3.11 reorganise_prelims_dir.py $line
        if [[ $? -ne 0 ]]; then
            log_error "Python script reorganise_prelims_dir.py failed for line: $line"
            continue
        fi

        # Add outdir to found_outdirs dictionary only if it is not already present
        if [[ -z "${found_outdirs[$outdir]}" ]]; then
            found_outdirs["$outdir"]=1
        fi

    done < "$in_par_prelims_file"

    # Compress all found outdirs
    for outdir in "${!found_outdirs[@]}"; do
        compress_dir "$outdir"
    done

    # Increment iteration and wait for the next time gap
    current_iter=$((current_iter + 1))
    if [[ $current_iter -lt $stop_iter ]]; then
        # First echo all entries in renamed_outdirs
        for key in "${!renamed_outdirs[@]}"; do
            echo "Renamed ${key} to ${renamed_outdirs[$key]}"
        done
        # Then reset the renamed_outdirs and found_outdirs dictionaries for the next iteration
        unset renamed_outdirs
        unset found_outdirs
        echo "Waiting for $time_gap_hours hours before the next iteration...."
        sleep "$time_gap_seconds"
    fi
done

echo "Script completed after $stop_iter iterations."