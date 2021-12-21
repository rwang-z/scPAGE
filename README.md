# scPAGE: Transfer gene signatures from single cells to bulk data

## Introduction

scPAGE is a method to identify gene expression signatures from single cells which were transferable and applicable to bulk data analysis. 

For more details, please read our paper: Improving bulk RNA-seq classification of acute myeloid leukemia by gene signature transfer from single cells. 

![image](https://user-images.githubusercontent.com/57746198/146364178-942adbfc-dd81-47c7-bb1e-4ba1869555e3.png)


## Prerequisite
Python 3.6

Python packages: pandas, numpy, matplotlib, statsmodels, fisher, sklearn, time, tqdm, collections, sys

Platform: Linux


## Identification of scGPS from single cells

### Data

Expression profile of single cells, e.g., scRNA-seq data, and corresponding labels are required as input of the method.

Normalization of the expression data is not required.


- **Example data**

An example training set including the expression level of 200 single cells is provided in './data/casemay10_example.txt' (100 cases) and './data/ctrlmay10_example.txt' (100 controls) for demonstration. The whole dataset used in the paper is available in Gene Expression Omnibus (GEO) with the accession number GSE128423 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128423). 

The mouse luekemia-related genes retrieved from NCBI to filter the expression profiles are listed in './data/NCBI_leukemia_mm_gene.txt'.



- **To use other datasets for scGPS identification, please prepare the required files as follows:**

**1. Expression profiles**

The training expression data could be provided by either a single file or a sequence of files. Each file contains an expression matrix where rows are genes and columns are cells. 

Replace the paths in the column 'path' of './data/train_data_path.txt' with that of your expression profiles for training. The files will be merged later for scGPS identification.

**2. Labels**

Label 0 and 1 represent normal and disease, respectively. 

Labels could be provided in two ways according to the expression profiles:

(1) If the cells/samples in the same expression profile belong to different classes, a separate label file (tab-separated) is required for each expression profile, in which the labels of the cells/samples are listed with the column name 'Label' in the same order as the cell/samples in the corresponding expression profile. List the paths of the label files in './data/train_label_path.txt' with the column name 'path'. 

(2) If all of the cells/samples in the same expression profile belong to the same class, users could use a list of labels where each element indicates the label for an expression profile. The label list should be stored with the column name 'hint_list' in './data/train_label_list.txt'. 

In both of the methods, please make sure that the paths of the label files or the labels in the list are sequenced in the same order as the corresponding expression profiles listed in './data/train_data_path.txt'.

**3. Gene list**

A file (tab-separated) containing a column 'Symbol' of a list of genes used to filter the expression profile. We recommend to use the related genes of the disease. 

Mouse leukemia-related genes used for the example data are listed in './data/NCBI_leukemia_mm_gene.txt'. 

We also provide a list of human immune genes retrieved from InnateDB in 'innateDB_immune_genes_human.txt'. 


### Usage

```
> python scGPS.py data_flag gene_path label_flag [separator] [max_size] [pvalue_thresh]
```

- data_flag:

  - 'example': using the example data for training, loading files in './data/train_data_path_default.txt'

  - 'user-provided': using the data provided by user, loading files in './data/train_data_path.txt'

- gene_path: the path of the gene list file used to filter the expression profiles. 

- label_flag:

  - 'path': using separate label files containing the labels of the samples/cells for the expression profiles and list the paths in './data/train_label_path.txt'.

  - 'list': using a list of labels (stored in './data/train_label_list.txt') where each element indicates the label for an expression profile.

- separator: separator used in loading profile matrix. Default '\t'.

- max_size: maximum signature size. Default 50.

- pvalue_thresh: p-value threshold to filter the differentially expression gene-pairs (DEPs). Default 0.05.


### Examples

- Identify scGPS using the example data:

```
> python scGPS.py 'example'
```

- Using your own data and label files (using the mouse leukemia-related genes to filter the expression profiles):

  - Using separate files to indicate labels
  
  ```
  > python scGPS.py 'user-provided' './data/NCBI_leukemia_mm_gene.txt' 'path'
  ```
  
  - Using a label list
  
  ```
  > python scGPS.py 'user-provided' './data/NCBI_leukemia_mm_gene.txt' 'list'
  ```


### Results

The program will save the results in a file including the gene pairs selected for scGPS, the corresponding tags ('1' for 'disease-positive' and '-1' for 'normal-positive') and AUC after adding the pair to the scGPS.

- The results of the example data will be saved in './result/training_results_default.txt.

- The results of the user-provided data will be saved in ./result/training_results.txt.

- scGPS of mouse leukemia identifed in the paper is provided in './result/pair_list.txt'.

A curve demonstrating AUC with respect to different signature sizes is saved in './figure/training_auc.pdf' for your own data or in './figure/training_auc_default.pdf' for the example data.



## Validate scGPS on test data

### Data

1. A gene expression profile, where rows are genes and columns are samples.

2. A file containing a column 'Label' listing the labels of the samples in the expression profile.

Three example bulk RNA-seq datasets of mouse leukemia are provided:

- Expression profiles:

  - './data/GSE74690_exp_matrix.txt'.

  - './data/GSE78691_exp_matrix.txt'. 

  - './data/GSE119299_exp_matrix.txt'

- Labels:

  - './data/GSE74690_label.txt'

  - './data/GSE78691_label.txt'

  - './data/GSE119299_label.txt'


### Usage

```
> python prediction.py profile_path label_path output_path signature_flag num_pair [separator]
```

- profile_path: the path of expression profile for validation.


- label_path: the path of the label file.

- output_path: the path of the file to output the predicted scores of the test samples.

- signature_flag: the scGPS used for prediction

  - 'example': using the scGPS identified from the example training data (stored in './result/training_results_default.txt').

  - 'signature': using the scGPS identified in the paper (stored in './result/pair_list.txt').
  
  - 'user-provided': using the scGPS identified from the data provided by the user (stored in './result/training_results.txt').

- num_pair: the signature size of scGPS used for prediction. num_pair = 30 when using the scGPS identified in the paper. 

- separator: separator when loading data matrix, default '\t'.


### Examples

- Test the scGPS identified from the example data:

```
> python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'example' 30
```

- Test the scGPS identified from all single cells in the paper:

```
> python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'signature' 30
```

- Test the scGPS identified from the data provided by users:

```
> python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'user-provided' 30
```


### Results

The program will output the test AUC of the scGPS on the validation set, along with a file containing the labels and the predicted scores of the test samples. 

