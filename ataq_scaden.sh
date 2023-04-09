#!/bin/bash

TRAINING_FILE="$1_training_data.h5ad"
PREDICTION_FILE="$1_prediction_data.txt"

echo "-----Processing with Scaden-----"
scaden process $TRAINING_FILE $PREDICTION_FILE

echo "-----Training with Scaden-----"
scaden train processed.h5ad --steps 5000 --model_dir model

echo "-----Prediction with Scaden-----"
scaden predict --model_dir model $PREDICTION_FILE

echo "-----Evaluating Scaden predictions-----"
python3 evaluate_predictions.py $1

echo "-----Completed successfully-----"
