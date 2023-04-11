# ATAQ-seq Cell Type Deconvolution using Deep Learning
Cell type deconvolution of ATAQ-seq data using the Scaden model from Menden et al.

## Process raw scATAQ-seq data for use with Scaden model
`./process_ataq_data_for_scaden <raw_data_filename> <genome-filename> <metadata-filename> <sample-name> <num-samples-training> <num-samples-prediction> <num-cells-per-sample> <chr-of-interest>`

### Arguments
`<raw-data-filename>`
Raw ATAQ-seq data file for a single sample (e.g. 'GSE184462_RAW/GSM5589344_adipose_omentum_SM-ADYHB_rep1_fragments.bed.gz').

`<genome-filename>` 
Genome reference file. The genome reference used here is GRCh38 (e.g. 'ucsc_genome_data.bed').

`<metadata-filename>`
Metadata file with single cell names and the corresponding sample they are taken from and the cell type they are labeled as (e.g. 'metadata.csv').

`<sample-name>` 
Name of the sample to be processed (e.g. 'adipose_omentum_SM-ADYHB').

`<num-samples-training>` 
Number of bulk samples to be generated for training the Scaden model.

`<num-samples-prediction>` 
Number of bulk samples to be generated for prediction with the Scaden model.

`<num-cells-per-sample>`
Number of cells to be used in each bulk sample.

`<chr-of-interest>`
OPTIONAL. Specific chromosomes of interest. Formatted with '_' between chromosome names (e.g. 'chr1_chr5'). If interested in all chromosomes, leave this argument empty.

This script takes approximately one hour to run, dependent on the size of the raw ATAQ-seq data.

## Train, predict, and evaluate Scaden model with ATAQ-seq data
`./ataq_scaden <sample-name>`

### Arguments
`<sample-name>`
Sample name should be the same as the one used in the above command.

This script takes approximately 10 minutes to run, dependent on samples sizes and number of cells per sample.

# References
### Model
Menden, K., M. Marouf, S. Oller, A. Dalmia, D.S. Magruder, K. Kloiber, P. Heutink, and S. Bonn. 2020. Deep learning-based cell composition analysis from tissue expression profiles. Sci Adv. 6:eaba2619. 
https://github.com/KevinMenden/scaden

### Data
Zhang, K., J.D. Hocker, M. Miller, X. Hou, J. Chiou, O.B. Poirion, Y. Qiu, Y.E. Li, K.J. Gaulton, A. Wang, S. Preissl, and B. Ren. 2021. A single-cell atlas of chromatin accessibility in the human genome. Cell. 184:5985-6001 e5919. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184462
