#!/bin/bash

source input.dat # Source the input variables
set -e # Exit if any command fails (check_fold_input.py)

dir=$(pwd)

input="$dir/inputs/folding_${run_name}_input.csv"
python codes/check_fold_input.py $input
echo "continue ... "
cd protein_structures/
mkdir -p ${run_name}
cd ${run_name}

#########################################
#              ALPHAFOLD                #
#########################################
mkdir -p AF
mkdir -p AF/AF_res

colabfold_batch $input "$dir"/protein_structures/"$run_name"/AF

cd AF && $(cp `ls | grep rank_001 | grep .pdb | xargs` AF_res)
cd AF_res
for file in $(ls); do
    if [[ -f "$file" ]]; then 
        filename=$(basename "$file")
        new_filename=$(echo "$filename" | cut -d '_' -f 1)
        mv "$file" "$new_filename.pdb" 
    fi
done
cd ..
rm -f *. 
rm -r $run_name*
mv AF_res/* .
rm -r AF_res
cd ..

#########################################
#               DEEPFOLD                #
#########################################
mkdir -p DF
mkdir -p DF/DF_res
colabfold_batch --model-type deepfold_v1 $input "$dir"/protein_structures/"$run_name"/DF

cd DF && $(cp `ls | grep rank_001 | grep .pdb | xargs` DF_res)
cd DF_res
for file in $(ls); do
    if [[ -f "$file" ]]; then 
        filename=$(basename "$file")
        new_filename=$(echo "$filename" | cut -d '_' -f 1)
        mv "$file" "$new_filename.pdb" 
    fi
done
cd ..
rm -f *. 
rm -r $run_name*
mv DF_res/* .
rm -r DF_res
cd ..

#########################################
#              OMEGAFOLD                #
#########################################
mkdir -p OF

echo "" >> $input
while IFS=, read -r num id sequence
do
    if [[ $id != "id" ]]
    then
        echo ">${id}"
        echo "${sequence}"
    fi
done < $input > input.fasta

omegafold input.fasta "$dir"/protein_structures/"$run_name"/OF
# rm input.fasta
