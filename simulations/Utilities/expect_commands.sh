#!/bin/bash

# Function to log errors
log_error() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') ERROR: $1" >> error_log.txt
}

# Function to log updates
log_update() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') UPDATE: $1" >> error_log.txt
}

# Function to execute commands with automatic input
execute_commands() 
{
    local default_file="Automata_Inputs/addn_commands.txt"
    # Check if input file is provided as a command-line argument
    local input_file="${1:-$default_file}"

    # Notify the user which file is being used
    if [ "$input_file" = "$default_file" ]; then
        echo "No input file provided. Using default file: $default_file"
        # Check if the default file exists
        if [[ ! -f "$default_file" ]]; then
            echo "Default additonal commands file not found: $commands_file"
            log_error "Default Additonal Commands file not found: $commands_file"
        fi
    else
        echo "Using provided input file: $input_file"
        # Check if the input file exists
        if [[ ! -f "$input_file" ]]; then
            echo "Additonal Commands Input file not found: $input_file"
            log_error "Additonal Commands Input file not found: $input_file"
        fi
    fi

    # Temporary expect script
    expect_script="execute_commands.exp"

    # Clean up previous expect script if exists
    [ -f "$expect_script" ] && rm "$expect_script"

    # Start writing the expect script
    echo '#!/usr/bin/expect -f' > "$expect_script"
    echo 'set timeout -1' >> "$expect_script"
    # Initialize variables
    command=""
    input=""
    # Read the commands file line by line
    while IFS= read -r line || [ -n "$line" ]; do
        # First remove leading and trailing whitespaces from the line
        line=$(echo "$line"  | tr -d '\r' |  xargs)

        # Debugging: Display the read line
        echo "Processing line: $line"

        if [[ "$line" == "#CMD" ]]; then
            echo "Found command marker (#CMD)"
            read -r command
            command=$(echo "$command" | tr -d '\r' | xargs)
            # Expand `~` or `$HOME` in the command
            command=$(eval echo "$command")
            echo "Running command: $command"
        elif [[ "$line" == "#INPUT" ]]; then
            read -r input
            echo "Found input marker (#INPUT) With Input: $input"
            input=$(echo "$input" | tr -d '\r' |  xargs)
            

        # Write expect logic for each command and its input
        cat >> "$expect_script" <<EOL
# Run the command
spawn $command
expect {
    "*password*" {
    send "$input\r"
    interact
    }
eof
}   
EOL
        fi
        echo "Done"
    done < "$input_file"

    # Make the expect script executable
    chmod +x "$expect_script"

    # Run the expect script
    ./"$expect_script"

    # Clean up
    rm "$expect_script"

    log_update "Additional Commands executed successfully."
}

#execute_commands "$1"