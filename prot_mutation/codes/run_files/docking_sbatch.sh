#!/bin/bash
source input.dat # Source the input variables

clust=$(which sbatch) # check if cluster

if [ -z "$clust" ]; then

    if [ $AF == 1 ]; then
        ./codes/run_files/docking.sh AF
    elif [ $DF == 1 ]; then
        ./codes/run_files/docking.sh DF
    elif [ $OF == 1 ]; then
        ./codes/run_files/docking.sh OF
    fi
else
    if [ $AF == 1 ]; then
        sbatch -J "$run_name"_AF codes/run_files/docking.sh AF
    elif [ $DF == 1 ]; then
        sbatch -J "$run_name"_DF codes/run_files/docking.sh DF
    elif [ $OF == 1 ]; then
        sbatch -J "$run_name"_DF codes/run_files/docking.sh OF
    fi
fi