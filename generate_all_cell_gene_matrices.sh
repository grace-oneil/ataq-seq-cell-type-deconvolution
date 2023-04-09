#!/bin/bash

while read line
do
	echo $line
	python3 generate_cell_gene_matrix.py $1 $line $3
done < $2
