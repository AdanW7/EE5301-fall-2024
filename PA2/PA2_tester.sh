#!/bin/bash

mkdir -p results
PRG_NAME=./placer

run_test() {
    local circuit_name=$1
    local output_file="./results/${circuit_name}_out.txt"
    echo "==========================================="
    echo "Testing ${circuit_name}" 
    echo ""
    
    time $PRG_NAME "$circuit_name" | tee "$output_file"
    # Run the make run command with the file parameter
    time make run file="$circuit_name"
}

# Check if an argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <circuit_name>"
    exit 1
fi

# Pass the provided argument to run_test
run_test "$1"