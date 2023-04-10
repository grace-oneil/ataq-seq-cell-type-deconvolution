import os
import logging
import sys
import gc
import numpy as np
import pandas as pd
import anndata as ad
import collections
from rich.progress import Progress, BarColumn
import scanpy as sc
import sklearn.metrics
from collections import Counter

sample_name = ''.join(sys.argv[1])

scaden_preds = pd.read_csv("scaden_predictions.txt", sep="\t", header=0, index_col=0)
scaden_preds.reset_index(drop=True, inplace=True)

actual_pred_filename = str(sample_name) + "_prediction_ground_truths.csv"
actual_preds = pd.read_csv(actual_pred_filename, sep=",", header=0, index_col=0)

y_preds = np.array(scaden_preds)
y_true = np.array(actual_preds)

diff = y_preds - y_true

mse = sklearn.metrics.mean_squared_error(y_true, y_preds)
mae = sklearn.metrics.mean_absolute_error(y_true, y_preds)

print("MSE = " + str(mse))
print("MAE = " + str(mae))
