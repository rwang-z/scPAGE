# scPAGE: Transfer gene signatures from single cells to bulk data

## Introduction

Identify gene pair signature from single cells.

<img width="452" alt="image" src="https://user-images.githubusercontent.com/57746198/146361412-4087f10d-9130-467b-b125-ab2af7a68b8e.png">



## Identification of scGPS using your own data

### Data

Training data could be provided by one scRNA-seq profile or a sequence of profiles, where rows are genes and columns are cells. The files will be merged for scGPS identification. Replace the paths of your scRNA-seq profiles in './data/train_data_path.txt'.

Labels could be provided by two ways:

1) Provide a label file for each training profile. Replace the paths of all the label files in './data/train_label_path.txt' with the same order as the corresponding training profile path in './data/train_data_path.txt'.
2) Provide one label for each training profile. If all the cells in the same scRNA-seq profile have the same label, list the label indicators (0/1) for each profile in './data/label_list.txt' with the same order as the corresponding training profile path in './data/train_data_path.txt'.

Label 0 and 1 represent normal and disease, respectively. 

Example scRNA-seq data is provided for demonstration with paths listed in './data/train_data_path_default.txt'. There two profile files for normal and disease, respectively, each including 100 cells. Label list is provided in ‘./data/label_list_default.txt’.

A list of leukemia genes for profile filtering is provided in ‘./data/ NCBI_leukemia_mm_gene.txt’.



### Usage

```
python scGPS.py gene_path label_flag [separator] [max_size] [pvalue_thresh]
```

- gene_path: path of a gene list file to filter the profile. Genes listed in column 'Symbol'.
- label_flag: use label file ('path') or label indicator ('list') for each training profile
- separator: separator used in loading profile matrix. Default '\t'.
- max_size: maximum signature size. Default 50.
- pvalue_thresh: p-value threshold to filter DEPs. Default 0.05.


### Examples

- Using your own data and label files:

```
python scGPS.py './data/NCBI_leukemia_mm_gene.txt' 'path'
```

- Example data is provided for demonstration. There two profile files for normal and disease, respectively, each including 100 cells. Identify scGPS using example data by running:

```
python scGPS.py
```

### Results

Gene pairs selected for scGPS, with corresponding tags and AUC after adding the pair to the scGPS.

scGPS identified from the example data saved in './result/training_results_default.csv'.

scGPS identified from your own data saved in ./result/training_results_test.csv'.

scGPS in the paper provided in './result/pair_list.csv'.

A curve of AUC with respect to different signature size is saved in './figure/training_auc_test.pdf' for your own data or in './figure/training_auc_default.pdf' for the example data.


## Validate scGPS on test data

### Data
RNA-seq profile, where rows are genes and columns are samples.

Labels listed in a file with column name ‘Label’.

Example test data is provided:

Profile matrix in './data/GSE119299_exp_matrix.txt'.

Label in './data/GSE119299_label.txt'.


### Usage

```
python prediction.py profile_path label_path num_pair signature_path [separator]
```

- profile_path: path of data matrix for validation. Profile matrix: row -- gene, column -- cell, tab-seperated.
- label_path: path of data label. Labels listed in column 'Label'.
- num_pair: number of gene pairs to use.
- signature_path: path of file including gene pairs (signature) to test. If not provided or provided with '', scGPS identified from the example data is used.
- separator: separator used in loading profile matrix. Default '\t'.

### Examples

- Test scGPS identified from all single cells in the paper:

```
python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' 30 './result/pair_list.csv' ' '
```

- Test scGPS identified from the example data:
```
python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' 30 '' ' '
```

- Test scGPS identified from your own data:
```
python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' 30 './result/training_results_test.csv' ' '
```
