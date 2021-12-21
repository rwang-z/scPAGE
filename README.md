# scPAGE: Transfer gene signatures from single cells to bulk data

## Introduction

scPAGE is a method to identify **single-cell gene pair signatures (scGPSs)** which were transferable and applicable to bulk data analysis. 

For more details, please read our paper: **Improving bulk RNA-seq classification of acute myeloid leukemia by gene signature transfer from single cells**. 

Workflow of the mehtod:

![image](https://user-images.githubusercontent.com/57746198/146364178-942adbfc-dd81-47c7-bb1e-4ba1869555e3.png)


## Prerequisite
**Python 3.6**

**Python packages**: pandas, numpy, matplotlib, statsmodels, fisher, sklearn, time, tqdm, collections, sys

**Platform**: Linux


## Identification of scGPS from single cells

### Data

Expression profile of single cells (e.g., scRNA-seq data), the corresponding labels and a gene list to filter the expression profile are required as input of the method.

Normalization of the expression data is not required.


- **Example data**

An example training set including the expression level of 200 single cells is provided in **'./data/casemay10_example.txt'** (100 cases) and **'./data/ctrlmay10_example.txt'** (100 controls) for demonstration. The whole dataset used in the paper is available in Gene Expression Omnibus (GEO) with the accession number **GSE128423** (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128423). 



- **Prepare your own data for scGPS identification:**

**1. Expression profiles**

The training expression data could be provided by either a single file or a sequence of files. Each file contains an expression matrix where rows are genes and columns are cells. 

Replace the paths in the column **'path'** of **'./data/train_data_path.txt'** with that of your expression profiles for training. The files will be merged later for scGPS identification.

**2. Labels**

Label **0** and **1** represent normal and disease, respectively. 

Labels could be provided in two ways according to the expression profiles:

- If the cells/samples in the same expression profile belong to different classes, a separate label file (tab-separated) is required for each expression profile, where the labels of the cells/samples are listed in the column **'Label'** in the same order as the cell/samples in the corresponding expression profile. List the paths of the label files in the column **'path'** of **'./data/train_label_path.txt'**.

- If all of the cells/samples in the same expression profile belong to the same class, users could use a list of labels where each element indicates the label of the whole expression profile. The label list should be stored in the column **'hint_list'** of **'./data/train_label_list.txt'**. 

In both of the methods, please make sure that the paths of the label files or the labels in the list are sequenced in the same order as the corresponding expression profiles listed in './data/train_data_path.txt'.

**3. Gene list**

A tab-separated file containing a column **'Symbol'** of genes used to filter the expression profile. We recommend to use the related genes of the disease. Users could prepare their own gene lists of interest or use the following two gene lists provided along with the code:

-	-	The mouse leukemia-related genes retrieved from NCBI are listed in **'./data/NCBI_leukemia_mm_gene.txt'**. This gene list is used to filter the example data.. 

-	The human immune genes retrieved from InnateDB are listed in **'innateDB_immune_genes_human.txt'**. 


### Usage

```
> python scGPS.py data_flag gene_path label_flag [separator] [max_size] [pvalue_thresh]
```

- **data_flag**:

  - 'example': using the example data for training, loading expression files in './data/train_data_path_default.txt'

  - 'user-provided': using the data provided by user, loading expression files in './data/train_data_path.txt'

- **gene_path**: the path of the gene list file used to filter the expression profiles. 

- **label_flag**:

  - 'path': using separate label files where each file contains the labels of the samples/cells for one expression profile. Please replace the label file paths in './data/train_label_path.txt'.

  - 'list': using a list of labels (stored in './data/train_label_list.txt') where each element indicates the label for an expression profile.

- **separator**: separator used in loading profile matrix. Default '\t'.

- **max_size**: maximum signature size. Default 50.

- **pvalue_thresh**: FDR-corrected p-value threshold to filter the differentially expression gene-pairs (DEPs). Default 0.05.


### Examples

- To identify scGPS using the example data:

```
> python scGPS.py 'example'
```

- To use your own data (e.g., using the mouse leukemia-related genes to filter the expression profiles):

  - Using separate files to indicate labels
  
  ```
  > python scGPS.py 'user-provided' './data/NCBI_leukemia_mm_gene.txt' 'path'
  ```
  
  - Using a label list
  
  ```
  > python scGPS.py 'user-provided' './data/NCBI_leukemia_mm_gene.txt' 'list'
  ```


### Results

The program will save the results including the gene pairs selected for scGPS, the corresponding tags ('1' for 'disease-positive' and '-1' for 'normal-positive') and AUC after adding the pair to the scGPS.

- The results of the example data are saved in **'./result/training_results_default.txt'**.

- The results of the user-provided data are saved in **'./result/training_results.txt'**.

- scGPS of mouse leukemia identifed in the paper is provided in **'./result/pair_list.txt'**, which could be tested on external datasets.

A curve demonstrating AUC with respect to different signature sizes is saved in **'./figure/training_auc.pdf'** for your own data or in **'./figure/training_auc_default.pdf'** for the example data.



## Validate scGPS on test data

### Data

1. A gene expression profile, where rows are genes and columns are samples.

2. A file containing a column 'Label' listing the labels of the samples in the expression profile.

Make sure that the samples and labels are sequenced in the same order in the two files.

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

- **profile_path**: the path of the expression profile for validation.

- **label_path**: the path of the label file.

- **output_path**: the path of the file to output the predicted scores of the test samples.

- **signature_flag**: the scGPS used for prediction

  - 'example': using the scGPS identified from the example training data (stored in './result/training_results_default.txt').

  - 'signature': using the scGPS identified in the paper (stored in './result/pair_list.txt').
  
  - 'user-provided': using the scGPS identified from the data provided by the user (stored in './result/training_results.txt').

- **num_pair**: the signature size of scGPS used for prediction. num_pair = 30 when using the scGPS identified in the paper. 

- **separator**: separator when loading data matrix, default '\t'.


### Examples

- Test the scGPS identified from the example data on the provided dataset GSE119299:

```
> python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'example' 30
```

- Test the scGPS identified from the whole single-cell dataset in the paper:

```
> python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'signature' 30
```

- Test the scGPS identified from the data provided by users (after the training step described above):

```
> python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'user-provided' 30
```


### Results

The program will output the test AUC of the scGPS on the validation set, along with a file (path defined by the users) containing the labels and the predicted scores of the test samples. 

