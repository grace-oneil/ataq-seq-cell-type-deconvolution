# ATAQ-seq Cell Type Deconvolution using Deep Learning
Cell type deconvolution of ATAQ-seq data using the Scaden model from Menden et al.

## Process raw scATAQ-seq data for use with Scaden model
`./process_ataq_data_for_scaden <raw_data_filename> <genome_filename> <metadata_filename> <sample_name> <cell_types_filename> <num_samples_training> <num_samples_prediction> <num_cells_per_sample> <chr_of_interest>`

### Arguments
`<raw_data_filename>`
Raw ATAQ-seq data file for a single sample (e.g. 'GSE184462_RAW/GSM5589344_adipose_omentum_SM-ADYHB_rep1_fragments.bed.gz').

`<genome_filename>` 
Genome reference file. The genome reference used here is GRCh38 (e.g. 'ucsc_genome_data.bed').

`<metadata_filename>`
Metadata file with single cell names and the corresponding sample they are taken from and the cell type they are labeled as (e.g. 'metadata.csv').

`<sample_name>` 
Name of the sample to be processed (e.g. 'adipose_omentum_SM-ADYHB').

`<cell_types_filename>` 
File file all possible individual cell types within the chosen sample (e.g. 'adipose_omentum_SM-ADYHB_cell_type_names.csv').

`<num_samples_training>` 
Number of bulk samples to be generated for training the Scaden model.

`<num_samples_prediction>` 
Number of bulk samples to be generated for prediction with the Scaden model.

`<num_cells_per_sample>`
Number of cells to be used in each bulk sample.

`<chr_of_interest>`
OPTIONAL. Specific chromosomes of interest. Formatted with '_' between chromosome names (e.g. 'chr1_chr5').

## Train, predict, and evaluate Scaden model with ATAQ-seq data
`./ataq_scaden <sample_name>`

### Arguments
`<sample_name>`
Sample name should be the same as the one used in the above command.

# References
### Model
Menden, K., M. Marouf, S. Oller, A. Dalmia, D.S. Magruder, K. Kloiber, P. Heutink, and S. Bonn. 2020. Deep learning-based cell composition analysis from tissue expression profiles. Sci Adv. 6:eaba2619. 
https://github.com/KevinMenden/scaden

### Data
Zhang, K., J.D. Hocker, M. Miller, X. Hou, J. Chiou, O.B. Poirion, Y. Qiu, Y.E. Li, K.J. Gaulton, A. Wang, S. Preissl, and B. Ren. 2021. A single-cell atlas of chromatin accessibility in the human genome. Cell. 184:5985-6001 e5919. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184462
