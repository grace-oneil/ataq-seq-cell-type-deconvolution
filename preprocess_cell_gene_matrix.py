import os
import logging
import sys
import gc
import numpy as np
import pandas as pd
import anndata as ad
# from anndata import read_h5ad
import collections
from rich.progress import Progress, BarColumn
import scanpy as sc
from collections import Counter

cell_type_filename = ''.join(sys.argv[1])
genome_filename = ''.join(sys.argv[2])
sample_name = ''.join(sys.argv[3])
chromosomes_of_interest = ','.join(sys.argv[4:])
chromosomes_of_interest = chromosomes_of_interest.split(",")
if chromosomes_of_interest == "":
        chromosomes_of_interest = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

cell_types = pd.read_csv(cell_type_filename, sep="\t", header=None)
cell_types.columns = ["cell_type"]
cell_types
obs_names = pd.DataFrame()


genome = pd.read_csv(genome_filename, sep="\t", header=None)
g_header = ["g_chrom", "g_start", "g_end", "gene_name", "0", "1", "2", "3", "4", "5", "6", "7"]
genome.columns = g_header
chr_genome = genome[genome["g_chrom"].isin(chromosomes_of_interest)]
chr_genome.drop(["g_chrom", "g_start", "g_end", "0", "1", "2", "3", "4", "5", "6", "7"], axis=1, inplace=True)
chr_genome = chr_genome.sort_values(by=["gene_name"])
chr_genome.reset_index(drop=True, inplace=True)

print("Compiling matrix with all cell types")
directory = "cell_gene_matrices"
i=0
concat_ad = None
prev_filename = ""
cell_types.sort_values(by="cell_type", ignore_index=True, inplace=True)
sorted_filelist = sorted(os.listdir(directory))
for filename in sorted_filelist:
        f = os.path.join(directory, filename)
        f_pd = pd.read_csv(f, sep=',', header=None, on_bad_lines='skip')
        for j in range(f_pd.shape[0]):
                cur = pd.DataFrame([cell_types.iloc[i,0]])
                obs_names = pd.concat([obs_names, cur], axis=0, ignore_index=True)
        if i==0:
                prev_filename = f
        elif i==1:
                prev_f_pd = pd.read_csv(prev_filename, sep=',', header=None, on_bad_lines='skip')
                f_sc = sc.AnnData(f_pd, dtype=int)
                prev_f_sc = sc.AnnData(prev_f_pd, dtype=int)
                concat_ad = ad.concat([prev_f_sc, f_sc], join="inner")
        else:
                f_sc = sc.AnnData(f_pd, dtype=int)
                concat_ad = ad.concat([concat_ad, f_sc], join="inner")
        i+=1
obs_names.columns = ["cell_type"]
concat_ad.obs_names = obs_names["cell_type"].to_list()
concat_ad.var_names = chr_genome["gene_name"].to_list()
adata = concat_ad

print("Preprocessing cell gene matrix")
sc.pp.filter_genes(adata, max_cells=5) # remove genes that are available in more than 5 cells
sc.pp.filter_genes(adata, min_counts=1) # remove genes that are not available in at least one cell

sc.pp.filter_cells(adata, max_genes=500) # remove cells that have more than 500 genes available

sc.pp.normalize_total(adata)

adata = adata.T

print("Writing matrix to file")
output_filename = str(sample_name) + "_gene_cell_matrix.h5ad"
adata.write(output_filename, compression='gzip') # write output
print("Completed successfully")
