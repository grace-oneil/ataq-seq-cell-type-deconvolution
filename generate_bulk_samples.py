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

matrix_filename = ''.join(sys.argv[1])
sample_name = ''.join(sys.argv[2])
number_of_samples = int(''.join(sys.argv[3]))
num_pred_samples = int(''.join(sys.argv[4]))
number_cells_per_sample = int(''.join(sys.argv[5]))

matrix = sc.read_h5ad(matrix_filename)

cell_types = set(matrix.var_names.values.tolist())

occurrences_per_cell_types = Counter(matrix.var_names.values)
cell_type_counts = list(occurrences_per_cell_types.values())
total_cells = len(matrix.var_names)

matrix.var_names_make_unique(join='_')

def generate_random_num_per_cell_type(cell_types, total_num_cells_per_sample):
    random_cell_nums = np.random.rand(1, len(cell_types))
    sum_cells = np.sum(random_cell_nums)
    fraction_cells = random_cell_nums / sum_cells
    num_cell_types = fraction_cells * total_num_cells_per_sample
    num_cell_types = np.around(num_cell_types).astype(int)
    if (np.sum(num_cell_types)) < total_num_cells_per_sample:
        diff = total_num_cells_per_sample - np.sum(num_cell_types)
        for i in range(diff):
            n = np.random.randint(0, len(cell_types))
            num_cell_types[0,n] += 1
    elif (np.sum(num_cell_types)) > total_num_cells_per_sample:
        diff = np.sum(num_cell_types) - total_num_cells_per_sample
        i = 0
        while i < diff:
            n = np.random.randint(0, len(cell_types))
            if num_cell_types[0,n] > 0: # don't want to generate negative numbers
                num_cell_types[0,n] -= 1
                i += 1
    return num_cell_types

def generate_samples(rng, keys, cell_types, total_num_per_sample):
    sample_gene_cell_mat = np.zeros([matrix.X.shape[0], 0])
    num_cell_types = generate_random_num_per_cell_type(cell_types, total_num_per_sample)
    for cell_type_i, cell_type in enumerate(keys):
        num_cells = occurrences_per_cell_types[cell_type]
        variable_names = [cell_type + "_" + str(j) if j in range(1, num_cells) else cell_type for j in range(num_cells)]
        variable_names = np.array(variable_names)
        choices = rng.choice(variable_names, num_cell_types[0, cell_type_i]).tolist()
        options = sc.get.obs_df(matrix, keys=choices).to_numpy()
        if cell_type_i == 0:
            sample_gene_cell_mat = options
        else:
            sample_gene_cell_mat = np.concatenate((sample_gene_cell_mat, options), axis=1)
    sample_ground_truth = num_cell_types / np.sum(num_cell_types)
    return sample_gene_cell_mat, sample_ground_truth

def generate_bulk_samples(n_samples, cell_types, n_per_sample):
    bulk_samples = np.zeros([matrix.obs.shape[0], 1], dtype=int)
    ground_truths = dict()
    sample_names = []

    rng = np.random.default_rng()
    keys = list(occurrences_per_cell_types.keys())
    for i in range(n_samples):
        sample, sample_ground_truth = generate_samples(rng, keys, cell_types, n_per_sample)
        bulk_sample_values = sample.sum(axis=1, dtype='int')
        new_sample_name = "sample_"+str(i+1)
        if i == 0:
            bulk_samples[:,0] = bulk_sample_values
            ground_truths[new_sample_name] = sample_ground_truth
        else:
            bulk_sample_values = np.reshape(bulk_sample_values, (matrix.obs.shape[0], 1))
            bulk_samples = np.concatenate((bulk_samples, bulk_sample_values), axis=1)
            ground_truths[new_sample_name] = sample_ground_truth
        sample_names.append(new_sample_name)
    adata = sc.AnnData(X=bulk_samples, obs=matrix.obs, var=sample_names)
    return adata, ground_truths

print("Generating simulated bulk samples for training")
adata, ground_truth_dict = generate_bulk_samples(number_of_samples, cell_types, number_cells_per_sample)
adata.var.columns = ["sample_name"]

list_version_ground_truth_dict = dict()
for key, val in ground_truth_dict.items():
    [value] = list(val)
    list_version_ground_truth_dict[key] = value
    
varm_data = pd.DataFrame.from_dict(list_version_ground_truth_dict)

test = varm_data.explode(list(varm_data.columns.values), ignore_index=True)
test.index = sorted(cell_types)
test.index.name = "cell_types"
labels = test.T
labels.reset_index(drop=True, inplace=True)
new_obs = labels

new_uns = dict()
new_uns['cell_types'] = new_obs.columns.values
new_uns['unknown'] = np.array(['unknown'], dtype='O')

adata.obs.drop(["n_cells", "n_counts"], axis=1, inplace=True)
new_var = adata.obs

new_data = sc.AnnData(X = adata.X.T, var = new_var, obs = new_obs, uns = new_uns)

print("Generating simulated bulk samples for prediction")
prediction_data, pred_ground_truth_dict = generate_bulk_samples(num_pred_samples, cell_types, number_cells_per_sample)
prediction_data.var.columns = ["sample_name"]

list_version_ground_truth_dict_preds = dict()
for key, val in pred_ground_truth_dict.items():
    [value] = list(val)
    list_version_ground_truth_dict_preds[key] = value
    
varm_data2 = pd.DataFrame.from_dict(list_version_ground_truth_dict_preds)

pred_truths = varm_data2.explode(list(varm_data2.columns.values), ignore_index=True)
pred_truths.index = sorted(cell_types)
pred_truths.index.name = "cell_types"
pred_labels = pred_truths.T
pred_labels.reset_index(drop=True, inplace=True)
pred_new_obs = pred_labels

pred_truth_filename = str(sample_name) + "_prediction_ground_truths.csv"

pred_txt_file_data = pd.DataFrame(prediction_data.X)
pred_txt_file_data.set_index(new_data.var_names, inplace=True)
pred_txt_file_data.columns = ["sample_"+str(i) for i in range(pred_txt_file_data.shape[1])]
pred_txt_filename = str(sample_name) + "_prediction_data.txt"

print("Writing files")
output_filename = str(sample_name) + "_training_data.h5ad"
new_data.write_h5ad(output_filename, compression='gzip')
pred_new_obs.to_csv(pred_truth_filename, sep=',')
pred_txt_file_data.to_csv(pred_txt_filename, sep="\t", header=True, index=True)
print("Training data is in file " + output_filename)
print("Prediction data for model is in file " + pred_txt_filename)
print("Prediction data ground truth cell proportions are in file " + pred_truth_filename)
print("Completed successfully")
