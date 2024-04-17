#!/bin/bash

source input.dat # Source the input variables
export PYTHONPATH=$DD_path #path to where DiffDock is installed
dir=$(pwd)
dd_dir=$DD_path

#SBATCH -J $run_name
#SBATCH -e error_%j.err
#SBATCH -o output_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1
##SBATCH --constraint=gpu8

cd $dd_dir

if [ $DD == true ]; then

    echo "Starting Diffdock..."
    date

    python -W ignore inference.py \
    --protein_ligand_csv $dir/output/protein_ligand.csv \
    --out_dir $dir/output/results_dd/ \
    --inference_steps 20 \
    --samples_per_complex $samples \
    --batch_size 10 \
    --actual_steps 18 \
    --no_final_step_noise \
    --model_dir workdir/paper_score_model \
    --confidence_model_dir workdir/paper_confidence_model
fi

if [ $DD_L == true ]; then

    echo "Starting Diffdock-L..."
    date

    python -W ignore inference.py \
    --config default_inference_args.yaml \
    --protein_ligand_csv $dir/output/protein_ligand.csv \
    --out_dir $dir/output/results_dd/ 
fi

cd $dir