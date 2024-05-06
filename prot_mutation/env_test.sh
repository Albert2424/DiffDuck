#!/bin/bash

# Define the name of the Conda environment
env_name="test_env"

# Define the list of packages to install
packages=("numpy=1.23.5" "matplotlib=3.8.3" "pandas=2.2.1" "mdtraj=1.9.9")

# Check if the Conda environment exists
if ! conda env list | grep -q "\<$env_name\>"; then
    echo "Creating Conda environment: $env_name"
    # Create the Conda environment
    conda create --name "$env_name" python=3.9.18
    conda config --add channels conda-forge
    
    # Activate the environment
    conda activate "$env_name"
    # Install packages into the environment
    for pkg in "${packages[@]}"; do
        conda install -y "$pkg"
    done
    # conda install -y $packages

else
    echo "Conda environment $env_name already exists."
    conda list
fi
