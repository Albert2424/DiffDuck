#!/bin/bash

source env.sh #activate env if necessary

result_list=$(ls DD)

python codes/DB_cleaning.py

for result in $result_list; do
    python codes/rank_formula.py --results_file DD/"$result" --output_file pred_"$result"
done

python codes/plots.py
