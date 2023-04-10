#!/bin/bash

echo "-----Processing with Scaden-----"
scaden process "$1_training_data.h5ad" "$1_prediction_data.txt"

echo "-----Training with Scaden-----"
scaden train processed.h5ad --steps 5000 --model_dir "$1_model"

echo "-----Prediction with Scaden-----"
scaden predict --model_dir "$1_model" "$1_prediction_data.txt"

echo "-----Evaluating Scaden predictions-----"
python3 evaluate_predictions.py $1

mkdir "$1_scaden"
mv scaden_predictions.txt "$1_scaden"
mv "$1_prediction_data.txt" "$1_scaden"
mv "$1_prediction_ground_truths.csv" "$1_scaden"
mv "$1_training_data.h5ad" "$1_scaden"
mv "$1_gene_cell_matrix.h5ad" "$1_scaden"
mv "$1_model" "$1_scaden"
mv cell_gene_matrices "$1_scaden"
mv cell_type_files "$1_scaden"
mv genome_comparison_files "$1_scaden"
mv processed.h5ad "$1_scaden"
mv "$1_cell_type_names.csv" "$1_scaden"

echo "-----Completed successfully-----"
