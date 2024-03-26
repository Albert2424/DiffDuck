# DiffDuck

This repository contains all the necessary files to run DiffDuck simulations along with the programs to reproduce the results in (link to DiffDuck paper?).

## Table of contents

## Requirements

## Contains

## Install

## Execution

## Results

## Reproduce figures

First of all clone the directory containing all the necessary files by doing:

```bash
git clone https://github.com/Albert2424/DiffDuck/tree/master/DB_example
```
Then move to the downloaded directory and use `source run.sh` to generate the output directory which includes both the figures and the output files. The code is split into three main parts. The first part consists of cleaning the selected database (in the example is `Database_example.zip` which is automatically unzipped when you run the code). After running `DB_cleanning.py` from the `/codes` directory, the following files are generated inside the `/output` directory:


1. `AF_input.csv`       --> In case you don't have the PDB structure of the studied proteins, it may be useful to run a folding program (such as AlphaFold) to obtain them. This file is a valid input for AlphaFold, DeepFold...

2. `data_prot.csv`      --> Contains the selected data from the database (i.e. all data that has a valid $K_i$ value). In this data, you may encounter sequences used in the DiffDock training so you may consider excluding them.

Here the code needs a file containing the result of a blast in the _PDBBind_ database to exclude said proteins. In case the user doesn't want to exclude said proteins, just comment the `drop_if_in_training('pdbbind_match_example.csv',df)` line of the `DB_cleanning.py` code. In the other case, a file with a format like _pdbbind_match_example.csv_ where the only column it contains are the indices of the rows to drop will be needed. Whichever is decision taken, the following files will be generated:

3. `data_A.csv`         --> Paired with data_B.csv will provide the different pairs of competing proteins for every ligand.

4. `data_B.csv`

5. `clean_data.csv`     --> This file will only appear if proteins from the training database have been excluded. Is the same file as _data_prot.csv_ but excluding the proteins.

Now everything is almost ready to run the DiffDock simulations. First of all an input for DiffDock is generated using ... \[JUAN ESCRIU AQUESTA PART :)\]

... The following file is generated after the DiffDock simulations:

6. `combined.csv` --> A file containing the ID, confidence and chain to which every sample has been attached.

This file must be included in the `/DD` directory that will appear in the directory where you are executing the `run.sh` file. If not then the program will ask you to add the file inside the corresponding directory. It is allowed to have more than one file in the `/DD` directory and all of them will be evaluated.

The program will now call `rank_formula.py` from the `/codes` directory to evaluate all the samples and establish whether the result of the docking simulation is protein A or B. In addition, a plot showing the reliability of every guess with $log(k_i^A/k_i^B)$ is provided. 

Finally, `plots.py` from the `/codes` directory will generate the rest of the plots that are placed inside `output/figures`. 
