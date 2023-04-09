#!/bin/bash

while read line
do
	echo $line
	python3 generate_cell_type_data_for_sample.py $1 $2 $3 $line $5
done < $4
