import os
import logging
import sys
import gc
import numpy as np
import pandas as pd
from anndata import read_h5ad
import collections
from rich.progress import Progress, BarColumn
import scanpy as sc

sample_name = ''.join(sys.argv[1])
cell_type = ''.join(sys.argv[2])
genome_filename = ''.join(sys.argv[3])
chromosomes_of_interest = ','.join(sys.argv[4:])
chromosomes_of_interest = chromosomes_of_interest.split(",")
if chromosomes_of_interest == "":
	chromosomes_of_interest = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
#print("Sample Name: " + sample_name)
#print("Cell Type: " + cell_type)
#print("Chromosomes of Interest: " + str(chromosomes_of_interest))

print("Reading files")
genome = pd.read_csv("grch38.bed", sep="\t", header=None)
g_header = ["g_chrom", "g_start", "g_end", "gene_name"]
genome.columns = g_header
chr_genome = genome[genome["g_chrom"].isin(chromosomes_of_interest)]

cell_type = cell_type.replace(" ", "")
read_filename = "./cell_type_files/"+str(sample_name)+ "_cell_type_" + str(cell_type) + ".csv"

test = pd.read_csv(read_filename, sep=",")

gene_accessibility_sample = chr_genome['gene_name']
gene_accessibility_sample = gene_accessibility_sample.to_frame()

test.sort_values(by=['name', 'gene_name'], inplace=True)

test_range = test['name'].unique()
test_range = test_range.tolist()

num_chr_genes = chr_genome.shape[0]
num_cells = len(test_range)

cell_gene_mat_sample = np.zeros([num_cells, num_chr_genes], dtype=int)
cell_gene_mat_sample.shape

cell_gene_mat_gene_labels = gene_accessibility_sample.sort_values(by=["gene_name"])
cell_gene_mat_gene_labels.reset_index(drop=True, inplace=True)

for i,value in enumerate(test_range):
    cur = test[test['name'] == value]
    for row in cur.itertuples():
        gene = row.gene_name
        idx = cell_gene_mat_gene_labels[cell_gene_mat_gene_labels['gene_name']==gene].index.to_list()[0]
        cell_gene_mat_sample[i][idx] = 1
print("Writing matrix to file")
os.system("mkdir cell_gene_matrices")
output_filename = "./cell_gene_matrices/" + str(sample_name) + "_" + str(cell_type)+"_cell_gene_mat.csv"
pd.DataFrame(cell_gene_mat_sample).to_csv(output_filename, sep=",", header=None, index=None)
print("Completed successfully")
