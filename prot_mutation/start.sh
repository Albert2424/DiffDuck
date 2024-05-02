#!/bin/bash
source codes/run_files/env.sh # Activate env if necessary
source input.dat #source the input variables   

if [ -z $SMILES ]; then 
    echo "Use input.dat to introduce the SMILES of the ligand to start."
fi

if [ -z $run_name ]; then
    random_int=$((RANDOM % 9000 + 1000))
    run_name="job_$random_int"
    sed -i "s/^run_name=.*/run_name=\"$run_name\"/" input.dat
    echo "Assigning random name to job since none was provided in input.dat."
fi
echo ""
echo ""

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "                $run_name"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo ""
alias peligan=make
make $@