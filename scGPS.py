# coding: utf-8
# To identify scGPS from single-cell RNA-seq data

### training data and labels
# put paths of your training profiles in './data/train_data_path.txt'
# training profile: rows -- genes, columns -- cells
# labels: 0 -- normal, 1 -- disease
# labels of training data could be provided in two ways:
#   -- 'path': the labels of each training profile is provided in a file, 
#              list all the path of label files in './data/train_label_path.txt'
#              in the same order with training profiles
#   -- 'list': if the samples/cells in the same training profile have the same label, 
#              then a indicator (0/1) of label for each training profile could be provided by a list,
#              put the list in './data/label_list.txt'

### Input arguments
# -- 0: path of gene list to filter the profile, genes in column 'Symbol'
# -- 1: use label path ('path') or indicator ('list') for each train matrix
# -- 2: separator in loading profile matrix, default '\t'
# -- 3: maximum signature size, default 50
# -- 4: p-value threshold to filter DEPs, default 0.05

### using default example data
# python scGPS.py

### using provided data
# python scGPS.py './data/NCBI_leukemia_mm_gene.txt' 'list' ' '
# python scGPS.py './data/NCBI_leukemia_mm_gene.txt' 'path' ' ' 50 0.8

import sys
import numpy as np
import load_data
import utils
import step_forward

def main():
    # input arguments processing
    if len(sys.argv) > 1:
        # using provided data
        args = sys.argv[1:]
        gene_path = args[0]
        label_flag = args[1]
        if len(args) > 2:
            separator = args[2]
        else:
            separator = '\t'
        if len(args) > 3:
            pair_limit = int(args[3])
        else:
            pair_limit = 50
        if len(args) > 4:
            pvalue_threshold = float(args[4])
        else:
            pvalue_threshold = 0.05
        file_path_flag = 'defined'
    else:
        # using example data
        gene_path = './data/NCBI_leukemia_mm_gene.txt'
        label_flag = 'list'
        separator = ' '
        pair_limit = 50
        pvalue_threshold = 0.05
        file_path_flag = 'default'

    # load data
    print('Loading data...')
    train_matrix_path,train_label_path,hint_list = load_data.load_path(label_flag,file_path_flag)
    train_matrix, train_label = load_data.train_data_preprocessing(train_matrix_path,train_label_path,hint_list,gene_path,separator,label_flag)
    train_label = np.array(train_label)
    num_gene = train_matrix.shape[1]
    num_cell = train_matrix.shape[0]

    # identify DEPs
    a,b,c,d = utils.get_value(num_gene, num_cell, train_matrix, train_label)
    pvalue_matrix = utils.fisher_index_mat(a,b,c,d,num_gene)
    j,k = np.where(pvalue_matrix <= pvalue_threshold)
    index = np.array([j,k],dtype = np.uint16)
    num_pair = len(j)
    print('%d DEPs are selected by the pvalue threshold' % num_pair)

    # select gene pairs by step forward algorithm
    res_df = step_forward.get_optimized_pair(train_matrix,train_label,pvalue_threshold,index,num_pair,a,b,c,d,pair_limit)

    # save results to csv
    if file_path_flag == 'defined':
        res_df.to_csv('./result/training_results_test.csv',sep='\t',index=False)
    elif file_path_flag == 'default':
        res_df.to_csv('./result/training_results_default.csv',sep='\t',index=False)
    
    return res_df


if __name__ == "__main__":
   main()
