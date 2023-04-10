import os
import pandas as pd
import sys
import collections

sample_name = ''.join(sys.argv[1])

df = pd.read_csv("metadata.csv", sep=',')
df.drop(["1","2","3","4"], axis=1, inplace=True)
df = df[df["sample"]==sample_name]
cell_types = [t.replace(' ','') for t in df["cell_type"].unique().tolist()]
cell_types = [t.strip('\"') for t in cell_types]
cell_types = sorted(cell_types)
cell_type_df = pd.DataFrame(cell_types)
output_filename = sample_name + "_cell_type_names.csv"
cell_type_df.to_csv(output_filename, sep="\t", header=None, index=False)
print("Completed successfully")
