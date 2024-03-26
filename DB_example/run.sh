#!/bin/bash

source env.sh # Activate env if necessary

python codes/DB_cleaning.py # Clean database

result_list=$(ls DD)

echo ''

if [ -z "$result_list" ]; then
    echo "/DD is empty."
    echo " --> Use output/Af_input.csv to run the predicted structures for the proteins"
    echo " --> Use the predicted proteins to generate the DD predictions."
    echo " --> If you already have a results file (i.e results_DD_AF.csv) place it in the /DD directory."
    echo " --> When there is at least one file in the /DD directory the code will generate the ranking and the plots."

fi

# Rank the predictions and plot ki ratio vs rel
for result in $result_list; do
    python codes/rank_formula.py --results_file DD/"$result" --output_file pred_"$result"
done

echo ''

pred_list=$(ls output | grep pred_*)

# plot the rest of the plots
INDEX=0
for pred in $pred_list; do
    if [ $INDEX == 0 ]; then
        python codes/plots.py --data_a 'output/data_A.csv' --data_b output/data_B.csv --threshold 0.1 --predictions output/"$pred" --counts True
    else
        python codes/plots.py --data_a 'output/data_A.csv' --data_b output/data_B.csv --threshold 0.1 --predictions output/"$pred"
    fi
    let INDEX=${INDEX}+1
done