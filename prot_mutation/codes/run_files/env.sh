#!/bin/bash

directory=$(pwd)

# Activate the conda environment for DD if necessary

if [[ $(conda env list | grep "*" | awk '{print$1}') != 'diffdock' ]]; then
    echo "Activating conda environment"
	conda activate diffdock
fi
echo ''
echo '----------------------------------------------------------------'
# create necessary directories

# Check if /output directory exists
if [ -d "proteins" ]; then
    # If it exists, erase its contents
    echo -e  "The \e]8;;file:$directory/proteins\a proteins \e]8;;\a directory contains:"
    echo "$(ls -hs proteins/protein_structures)"
    echo "$(ls -hs proteins/paired_structures)"  
    # rm -r output/*
    # mkdir output/figures
else
    # If it doesn't exist, create it
    mkdir proteins
    mkdir proteins/protein_structures
    mkdir proteins/paired_structures
fi

echo '----------------------------------------------------------------'
echo ''
echo '----------------------------------------------------------------'

if [ -d "inputs" ]; then
    echo -e  "The \e]8;;file:$directory/inputs\a inputs \e]8;;\a directory are:"
    echo "$(ls -hs inputs)"
else
    # If it doesn't exist, create it
    mkdir inputs
fi
echo '----------------------------------------------------------------'
echo ''
echo '----------------------------------------------------------------'

if [ -d "results" ]; then
    echo -e  "The \e]8;;file:$directory/results\a results \e]8;;\a directory are:"
    echo "$(ls -hs results)"
else
    # If it doesn't exist, create it
    mkdir results
fi
echo '----------------------------------------------------------------'
