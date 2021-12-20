# coding: utf-8
# To identify scGPS from single-cell RNA-seq data

### Data preparation: see https://github.com/rwangsunshine/scPAGE
# Paths of expression profiles for training listed in: 
#       -- demonstration data: './data/train_data_path_default.txt'
#       -- user-provided data: './data/train_data_path.txt'
# Labels (paths) listed in './data/train_label_path.txt' (paths) or './data/train_label_list.txt' (indicators)

### Usage
# python scGPS.py data_flag gene_path label_flag [separator] [max_size] [pvalue_thresh]

### Input arguments
# data_flag: 
#       -- 'example': using the example data for training, loading files in './data/train_data_path_default.txt'
#       -- 'user-provided': using the data provided by user, loading files in './data/train_data_path.txt'
# gene_path: path of a gene list file to filter the profile. 
#       The file './data/NCBI_leukemia_mm_gene.txt' including the mouse leukemia-related genes retrieved from NCBI is used when identifying scGPS from the example data. 
#       We also provide a list of human immune genes retrieved from InnateDB in './data/innateDB_immune_genes_human.txt'.
#       Other genes could be used by providing the path of the file (tab-separated). Genes shouled be listed with the column name 'Symbol'. 
# label_flag: 
#       -- 'path': using label files containing the labels of samples/cells in the expression profiles. 
#                   List the paths of the files in './data/train_label_path.txt' under the column 'path' with the same order as the training profile files listed in './data/train_data_path.txt'.
#       -- 'list': using a list of labels where each element indicates the label for a training profile. The labels should be listed in './data/train_label_list.txt'.
# separator: separator used in loading profile matrix. Default '\t'.
# max_size: maximum signature size. Default 50.
# pvalue_thresh: p-value threshold to filter DEPs. Default 0.05.

### Examples
### using the example data
# python scGPS.py 'example'

### using the user provided data (using the mouse leukemia-related genes to filter the expression profiles as example)
### using label files
# python scGPS.py 'user-provided' './data/NCBI_leukemia_mm_gene.txt' 'path'

### using label list
# python scGPS.py 'user-provided' './data/NCBI_leukemia_mm_gene.txt' 'list'

import sys
import numpy as np
import load_data
import utils
import step_forward

def main():
    ######### input arguments processing #########
    arg_length = len(sys.argv)
    if arg_length == 1:
        print('Please provide at least one argument to indicate the data used for training')
        exit()
    
    args = sys.argv[1:]
    file_path_flag = args[0]

    # default setting
    print('Training on the %s data' % file_path_flag)
    gene_path = './data/NCBI_leukemia_mm_gene.txt'
    label_flag = 'list'
    separator = '\t'
    pair_limit = 50
    pvalue_threshold = 0.05
    res_path = './result/training_results_default.txt'

    if file_path_flag == 'user-provided':
        # using the data provided by the user
        res_path = './result/training_results.txt'
        gene_path = args[1]
        label_flag = args[2]
        if len(args) > 3:
            separator = args[3]
            print('Setting the separator of loading expression profile to %s' % separator)
        if len(args) > 4:
            pair_limit = int(args[4])
            print('Setting the maximum number of pairs to %d' % pair_limit)
        if len(args) > 5:
            pvalue_threshold = float(args[5])
            print('Setting the threshold of differentially expressed gene-pairs to %.2f' % pvalue_threshold)
    
    print('Genes to filter the profiles: %s' % gene_path)
    print('Using the %s of labels' % label_flag)

    ######### load data #########
    print('Loading data...')
    train_matrix_path,train_label_path,hint_list = load_data.load_path(label_flag,file_path_flag)
    train_matrix, train_label = load_data.train_data_preprocessing(train_matrix_path,train_label_path,hint_list,gene_path,separator,label_flag)
    train_label = np.array(train_label)
    num_gene = train_matrix.shape[1]
    num_cell = train_matrix.shape[0]

    ######### identify DEPs #########
    a,b,c,d = utils.get_value(num_gene, num_cell, train_matrix, train_label)
    pvalue_matrix = utils.fisher_index_mat(a,b,c,d,num_gene)
    j,k = np.where(pvalue_matrix <= pvalue_threshold)
    index = np.array([j,k],dtype = np.uint16)
    num_pair = len(j)
    print('%d DEPs are selected by the pvalue threshold' % num_pair)

    ######### select gene pairs by step forward algorithm #########
    res_df = step_forward.get_optimized_pair(train_matrix,train_label,pvalue_threshold,index,num_pair,a,b,c,d,pair_limit,file_path_flag)

    ######### save results to csv #########
    res_df.to_csv(res_path,sep='\t',index=False)
    print('scGPS saved to %s' % res_path)

if __name__ == "__main__":
   main()
