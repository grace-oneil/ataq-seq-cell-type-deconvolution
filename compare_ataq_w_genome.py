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
genome_filename = ''.join(sys.argv[2])
chromosomes_of_interest = ','.join(sys.argv[3:])
chromosomes_of_interest = chromosomes_of_interest.split(",")
if chromosomes_of_interest == "":
        chromosomes_of_interest = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]


print("---------Reading in files---------")
df = pd.read_csv(ataq_filename, sep='\t', comment='t', header=None)
df_genome = pd.read_csv(genome_filename, sep='\t', comment='t', header=None, on_bad_lines='skip')

header = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
df.columns = header[:len(df.columns)]
df = df.drop(['score', 'strand', 'name'], axis=1)
df = df[df["chrom"].isin(chromosomes_of_interest)]
df.end = df.end.astype(int)
df.start = df.start.astype(int)

header_genome = ['chrom', 'start', 'end', 'name', '?1','strand', 'cdsStart', 'cdsEnd', '?2', 'exonCount', 'exonStarts', 'exonEnds', 'proteinID', 'alignID']
df_genome.columns = header_genome[:len(df_genome.columns)]
df_genome = df_genome[df_genome["chrom"].isin(chromosomes_of_interest)]
df_genome.end = df_genome.end.astype(int)
df_genome.start = df_genome.start.astype(int)
df_gene_names = df_genome
df_genome = df_genome.drop(['?1', 'strand', 'name', 'cdsStart', 'cdsEnd', '?2', 'exonCount', 'exonStarts', 'exonEnds'], axis=1)


df.to_csv('ataq_bedtools.bed', sep ='\t', index=False, header=False)
df_genome.to_csv('genome_bedtools.bed', sep = '\t', index=False, header=False)


print("---------Comparing ATAQ data to genome data---------")
f = ataq_filename.split("/")
if len(f) >= 3:
	f = f[2]
elif len(f) == 1:
	f = f[0]
else:
	f = f[1]
#f = ataq_filename
os.system("mkdir genome_comparison_files")
new_filename = "./genome_comparison_files/"+f+".bed"
print("New filename = "+new_filename)
print("Bedtools Intersect")
os.system("bedtools intersect -a ataq_bedtools.bed -b genome_bedtools.bed -f 1.0 -wa -wb >> '{0}'".format(new_filename))

print("Completed successfully.")
