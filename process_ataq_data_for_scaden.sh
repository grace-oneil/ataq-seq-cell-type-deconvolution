#!/bin/bash

#1: ataq_filename
#2: genome_reference_filename
#3: metadata_filename
#4: sample_name
#5: number_of_training_samples
#6: number_of_prediction_samples
#7: number_of_cells_per_sample
#8: chromosomes_of_interest

echo "-----Creating cell type name file-----"
python3 generate_cell_type_names.py $4

echo "-----Comparing ATAQ data to genome-----"
python3 compare_ataq_w_genome.py $1 $2 $8 2> err.txt

echo "-----Generating cell type data for sample-----"
./generate_cell_type_data_for_sample.sh $1 $3 $4 "$4_cell_type_names.csv" $8 2> err.txt

echo "-----Generating cell gene matrices-----"
./generate_all_cell_gene_matrices.sh $4 "$4_cell_type_names.csv" $8 2> err.txt

echo "-----Preprocessing cell gene matrices-----"
python3 preprocess_cell_gene_matrix.py "$4_cell_type_names.csv" $2 $4 $8 2> err.txt

echo "-----Generating simulated bulk samples-----"
python3 generate_bulk_samples.py "$4_gene_cell_matrix.h5ad" $4 $5 $6 $7 2> err.txt

echo "-----Completed successfully-----"
