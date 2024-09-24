#!/bin/bash

# This script is used in conjunction with the Linux utility screen.
# It is used to run a series of jobs in parallel, each in a separate screen session.
# It accepts two arguments: the input file, and SpB (the number of biotic species).
# The input file is a tab delineated txt file with the following columns:
#dt	t_max	L	R	a_st	a_end	div	dP	PREFIXES
# where the seventh arguement on each line (corresponding to div) gives the number of threads.
# Each line in the input file corresponds to a separate job, which will be run in a separate screen session.
# An array of screen session names is provided, and the script will check if each screen session is running before starting a new one.
# If the screen session is not running, the script will start a new screen session and run the job.
# If the screen session is running, the script will skip to the next screen session name in the array.

screen_names=("casablanca" "tangiers" "tunis" "fes" "maghreb" "marrakesh" "rabat" "algiers" "medea" "gibraltar" "dakar")

# Check if the number of arguments is at least 1
if [ $# -lt 2 ]; then
    echo "Usage: $0 <Path/to/init_file.txt> <SpB>"
    exit 1
fi

# Check if the screen utility is installed.
if ! command -v screen &> /dev/null; then
    echo "Error: screen utility is not installed. Please install screen and try again."
    exit 1
fi

#Check if the input file exists
if [ ! -f $1 ]; then
    echo "Error: Input file not found. Are you sure the path is correct?. Exiting."
    exit 1
fi

# Store the current directory
curr_dir=$(pwd)

# Function that checks if starts a new screen session and runs the current job after first checking if a screen with the same name is already running.
# If a screen session with the same name is running, it will increment the screen_idx_mod variable to use the next screen name in the array.
# It will also issue a warning message to the user.
# The function takes the screen session name, the current screen_index, i as inputs, read the parameters from the input file, and starts the screen session.
# It will then execute the compiled program with the input parameters in that screen and finally return the updated screen_index.
start_screen(){
    local file_path=$1
    local spb=$2
    local screen_index=$3
    local i=$4
    # If screen_index exceeds the length of the screen_names array, issue an error message and exit
    if [ $((screen_index)) -ge ${#screen_names[@]} ]; then
        echo "Error: Insufficient screen names provided. Exiting."
        exit 1
    fi
    read -r p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 <<< $(sed -n "$((i+2))p" $1)
    #Remove trailing whitespaces and new lines from p9
    if [ $3 -eq 1 ]; then
        p9=\$(echo \$p9 | tr -d '[:space:]')
    else
        p10=\$(echo \$p10 | tr -d '[:space:]')
    fi
    if ! screen -list | grep -q "${screen_names[$screen_index]}"; then
        # Start a new screen session and run the job
        screen -dmS ${screen_names[$screen_index]} bash -c "cd .. ; g++-12 -DSPB=${spb} rietkerk_bjork_basic.cpp order_${spb}stoc_test_rietkerk.cpp -fopenmp -o ${screen_names[$screen_index]}_${p10}.out -std=c++23; ./${screen_names[$screen_index]}_${p10}.out $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 &> stderr_${screen_names[$screen_index]}.txt; cd $curr_dir"
        screen_index=$((screen_index+1))
    else
        #echo "Screen session ${screen_names[$screen_index]} is already running. Using the next screen name."
        screen_index=$((screen_index+1))
        # Recursively call the function with the updated screen_idx
        screen_index=$(start_screen $file_path $spb $screen_index $i)
    fi
    echo "$screen_index"
}

# Get the number of screen sessions from the number of lines in the input file
num_screens=$(wc -l < $1)
screen_idx=0

# Read the input file and start a screen session for each job
for ((i=0; i<num_screens; i++)); do
    screen_idx=$(start_screen $1 $2 $screen_idx $i)
    #echo "Screen index: $screen_idx"
    LINE=$(sed -n "$((i+2))p" $1)
    echo "Starting screen session ${screen_names[$((screen_idx-1))]} for job $i with parameters: $LINE"
done