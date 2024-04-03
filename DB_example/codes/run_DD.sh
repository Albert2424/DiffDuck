#!/bin/bash
#SBATCH -J DD_tot
#SBATCH -e error_%j.err
#SBATCH -o output_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1
##SBATCH --constraint=gpu8

./all_run.sh --protein_ligand input_DD.csv --replicas 1
