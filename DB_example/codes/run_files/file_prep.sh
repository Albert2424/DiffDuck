#!/bin/bash

# source codes/run_files/env.sh # Activate env if necessary
source input.dat # Source the input variables

python codes/input_merge.py --chain_A $data_A --chain_B $data_B
echo "--> Merge input file created"

python codes/file_preparation_multi.py --merge_input "output/merge_input.csv" --replicas $n_replica
echo "--> protein_ligand.csv created (Input for DiffDock)"
echo "--> File preparation done"