#!/bin/bash

source input.dat # Source the input variables

result_list=$(ls DD | sed 's/\..*//')

echo ''

if [ -z "$result_list" ]; then
    echo "/DD is empty."
    echo " --> Use output/Af_input.csv to run the predicted structures for the proteins"
    echo " --> Use the predicted proteins to generate the DD predictions."
    echo " --> If you already have a results file (i.e results_DD_AF.csv) place it in the /DD directory."
    echo " --> When there is at least one file in the /DD directory the code will generate the ranking and the plots."

fi

# Rank the predictions and plot ki ratio vs rel

if [ $model_name == "None" ]; then
    for result in $result_list; do
        python codes/rank_formula.py --results_file DD/"$result".csv --output_file pred_"$result" --data_a "$data_A" --data_b "$data_B" --model_name "$model_name"
    done
else
    for result in $result_list; do
        python codes/rank_formula.py --output_file pred_"$result".csv --data_a "$data_A" --data_b "$data_B" --model_name "$model_name" --n_samples "$n_samples" --validation_data training_data/"$result"_train/validation_data.csv --labels "$labels"
    done
fi