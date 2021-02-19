# scGPS

## Introduction

Identify gene pair signature from single cells.


## Identification of scGPS using your own data

### Data

Training data could be provided by one scRNA-seq profile or a sequence of profiles, where rows are genes and columns are cells. Replace the paths of your scRNA-seq profiles in './data/train_data_path.txt'.

Labels could be provided by two ways:

1) Provide a label file for each training profile. Replace the paths of all the label files in './data/train_label_path.txt' with the same order as the corresponding training profile path in './data/train_data_path.txt'.

3) Provide one label for each training profile. If all the cells in the same scRNA-seq profile have the same label, list the label indicators (0/1) for each profile in './data/label_list.txt' with the same order as the corresponding training profile path in './data/train_data_path.txt'.

Label 0 and 1 represent normal and disease, respectively. 

### Usage

python scGPS.py gene_path label_flag [separator] [max_size] [pvalue_thresh]

  -- gene_path: 









