#!/bin/bash

directory=$(pwd)

# Activate the conda environment for DD if necessary

if [[ $(conda env list | grep "*" | awk '{print$1}') != 'diffdock' ]]; then
    echo "Activating conda environment"
	conda activate diffdock
fi

# create necessary directories

# Check if /results directory exists
if [ -d "output" ]; then
    # If it exists, erase its contents
    rm -r output/*
    mkdir output/figures
else
    # If it doesn't exist, create it
    mkdir output
    mkdir output/figures
fi

# Iterate over each zip file found and unzip it
for zip_file in $(find "$directory" -type f -name "*.zip")
do
    # Unzip the file
    unzip "$zip_file" -d "${zip_file%.*}"
    echo "Unzipped: $zip_file"
done