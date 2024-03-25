#!/bin/bash

source env.sh #activate env if necessary

result_list=$(ls DD)

python codes/DB_cleaning.py

for result in $result_list; do
    python codes/rank_formula.py --results_file DD/"$result" --output_file pred_"$result"
done

pred_list=$(ls output | grep pred_*)
for pred in $pred_list; do
    python codes/plots.py --data_a 'output/data_A.csv' --data_b output/data_B.csv --threshold 0.1 --predictions output/"$pred"
done