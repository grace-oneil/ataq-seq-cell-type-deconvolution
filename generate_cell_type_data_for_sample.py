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

ataq_filename = ''.join(sys.argv[1])
end_file = ataq_filename.split('/')
end = end_file[-1]
intersect_filename = "genome_comparison_files/"+end+".bed"
metadata_filename = ''.join(sys.argv[2])
sample_name = ''.join(sys.argv[3])
cell_type = ''.join(sys.argv[4])
chromosomes_of_interest = ''.join(sys.argv[5:])
chromosomes_of_interest = chromosomes_of_interest.split("_")
if chromosomes_of_interest == [""]:
	chromosomes_of_interest = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
print("ATAQ_seq Data Filename: " + ataq_filename)
print("Genome Comp Filename: " + intersect_filename)
print("Sample Name: " + sample_name)
print("Cell Type: " + cell_type)
print("Chromosomes of Interest: " + str(chromosomes_of_interest))

print("Reading files")
ataq_data = pd.read_csv(ataq_filename, sep='\t', on_bad_lines='skip', header=None)
ataq_data.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
ataq_data = ataq_data.drop(['score', 'strand'], axis=1)

chr_data = ataq_data[ataq_data["chrom"].isin(chromosomes_of_interest)]

genome = pd.read_csv("grch38.bed", sep="\t", header=None)
g_header = ["g_chrom", "g_start", "g_end", "gene_name"]
genome.columns = g_header
chr_genome = genome[genome["g_chrom"].isin(chromosomes_of_interest)]

metadata = pd.read_csv(metadata_filename, sep=",", on_bad_lines="skip")

df_intersect = pd.read_csv(intersect_filename, sep='\t', comment='t', header=None, on_bad_lines='skip')

header_intersect = ['ataq_chrom', 'ataq_start', 'ataq_end', 'g_chrom', 'g_start','g_end']
df_intersect.columns = header_intersect
chr_intersect = df_intersect[df_intersect["ataq_chrom"].isin(chromosomes_of_interest)]

chr_accessibility_w_names = pd.merge(chr_intersect, chr_genome,  how='left', left_on=['g_chrom','g_start','g_end'], right_on = ['g_chrom','g_start','g_end'])

chr_accessibility_w_names.drop(["g_chrom", "g_start", "g_end"], axis=1, inplace=True)

chr_accessibility = pd.merge(chr_accessibility_w_names, chr_data,  how='left', left_on=['ataq_chrom','ataq_start','ataq_end'], right_on = ['chrom','start','end'])

chr_accessibility = pd.merge(chr_accessibility, metadata,  how='left', left_on=['name'], right_on = ['sample_name'])

chr_accessibility.drop(["ataq_chrom", "ataq_start", "ataq_end", "sample_name", "1", "2", "3", "4", "age"], axis=1, inplace=True)

chr_accessibility_sample = chr_accessibility[chr_accessibility['sample'] == sample_name]

chr_accessibility_sample.drop(["chrom", "start", "end"], axis=1, inplace=True)

cell_type_range = chr_accessibility_sample['cell_type'].unique()
cell_type_range = cell_type_range.tolist()

print("Writing data for individual cell types within sample")
os.system("mkdir cell_type_files")
for i,value in enumerate(cell_type_range):
	ctn = value.replace(" ", "")
	file = "./cell_type_files/"+str(sample_name)+"_cell_type_"+str(ctn)+".csv"
	chr_accessibility_sample[chr_accessibility_sample['cell_type'] == value].to_csv(file,index = False, na_rep = 'N/A')

print("Completed successfully")
