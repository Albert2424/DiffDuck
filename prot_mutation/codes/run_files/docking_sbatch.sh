#!/bin/bash
source input.dat # Source the input variables

clust=$(which sbatch) # check if cluster

if [ -z "$clust" ]; then
    ./codes/run_files/docking.sh
else
    sbatch -J $run_name codes/run_files/docking.sh
fi