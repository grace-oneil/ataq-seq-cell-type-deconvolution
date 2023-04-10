#!/bin/bash

echo "-----Processing with Scaden-----"
scaden process "$1_training_data.h5ad" "$1_prediction_data.txt"

echo "-----Training with Scaden-----"
scaden train processed.h5ad --steps 5000 --model_dir "$1_model"

echo "-----Prediction with Scaden-----"
scaden predict --model_dir "$1_model" "$1_prediction_data.txt"

echo "-----Evaluating Scaden predictions-----"
python3 evaluate_predictions.py $1

echo "-----Completed successfully-----"
