#!/bin/bash

protein_ligand=$1
csv_file=$2
replicas=$3
n_replicas=$4

# ./all_run.sh --protein_ligand input.csv --replicas nº_replicas


python /home/ramon/juan/important_files/file_preparation_multi.py $protein_ligand $csv_file $replicas $n_replicas
echo "Starting Diffdock..."

exec &> analysis.log
date

python -W ignore /home/ramon/progs/DiffDock/inference.py \
--protein_ligand_csv protein_ligand.csv \
--out_dir results/ \
--inference_steps 20 \
--samples_per_complex 5 \
--batch_size 10 \
--actual_steps 18 \
--no_final_step_noise \
--model_dir /home/ramon/progs/DiffDock/workdir/paper_score_model \
--confidence_model_dir /home/ramon/progs/DiffDock/workdir/paper_confidence_model 

analysis_counter=0
# complex_names=$(cat complex_names.txt)
result_dir=results/
# for complex in $complex_names; do
for folder in $result_dir*; do
    if [ -d "$folder" ]; then
        echo "Analysing "$folder""
        file_name=$(basename "$folder")
        /home/ramon/juan/important_files/sdf_to_pdb.sh "$folder"
        python /home/ramon/juan/important_files/ranking_afinity_def.py "${file_name}.pdb" "$folder" $analysis_counter
        analysis_counter=$((analysis_counter+1))
    fi



echo "Done"
mv results/"$file_name" .
mv "${file_name}.pdb" "$file_name"
echo "Results are in $"$file_name""
echo "###################################################################"
done

exec &> /dev/tty

mv protein_ligand.csv *log "$file_name"
mkdir clean_results
awk 'FNR==1 && NR!=1{next;}{print}' $(find . -name 'result.csv') > combined.csv
mv combined.csv clean_results

echo "Analysis finished"