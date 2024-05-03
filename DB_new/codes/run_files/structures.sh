#!/bin/bash

source input.dat # Source the input variables

input="/inputs/folding_${run_name}_input.csv"
cd prot_structures/
mkdir -p ${run_name}
cd ${run_name}

#########################################
#              ALPHAFOLD                #
#########################################
mkdir -p AF
colabfold_batch $input AF

for file in AF/*; do
    if [[ -f "$file" ]]; then 
        filename=$(basename "$file")
        new_filename=$(echo "$filename" | cut -d '_' -f 1)
        mv "$file" "$new_filename.pdb.w" 
    fi
done
 
find /AF -type f ! -name '*.w' -exec rm -f {} \;

for file in /AF/*.w; do
    mv "$file" "${file%.w}.pdb"
done

#########################################
#               DEEPFOLD                #
#########################################
mkdir -p DF
colabfold_batch --model-type deepfold_v1 $input DF

for file in DF/*; do
    if [[ -f "$file" ]]; then 
        filename=$(basename "$file")
        new_filename=$(echo "$filename" | cut -d '_' -f 1)
        mv "$file" "$new_filename.pdb.w" 
    fi
done
 
find /DF -type f ! -name '*.w' -exec rm -f {} \;

for file in /DF/*.w; do
    mv "$file" "${file%.w}.pdb"
done

#########################################
#              OMEGAFOLD                #
#########################################
mkdir -p OF

while IFS=, read -r id sequence
do
    if [[ $id != "id" ]]
    then
        echo ">${id}"
        echo "${sequence}"
    fi
done < $input > input.fasta

omegafold input.fasta OF
rm input.fasta