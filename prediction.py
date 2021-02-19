# coding: utf-8
# To validate scGPS on external data

# Input arguments
# -- 0: path of data matrix for validation
#       profile matrix: row -- gene, column -- cell, tab-seperated
# -- 1: path of data label
#       label file: a list of labels in column 'Label' 
# -- 2: number of gene pairs to use
# -- 3: path of selected gene pairs and corresponding tags
#       gene pairs provided in column 'pairs', tag provided in column 'tag'
#       if not provided or given '', use the default scGPS in './result/training_results_default.csv'
# -- 4: separator when loading data matrix, default '\t'


### use the scGPS identified from the example data, gene pairs in './result/training_results_default.csv'
# python prediction.py './data/GSE78691_exp_matrix.txt' './data/GSE78691_label.txt' 30 '' ' '

### use the scGPS identified from all single cells:
# python prediction.py './data/GSE78691_exp_matrix.txt' './data/GSE78691_label.txt' 30 './result/pair_list.csv' ' '

### use the identified scGPS in './result/training_results_test.csv'
# python prediction.py './data/GSE78691_exp_matrix.txt' './data/GSE78691_label.txt' 30 './result/training_results_test.csv' ' '

import sys
import pandas as pd
import numpy as np
import utils

def main():
    args = sys.argv[1:]
    test_data_path = args[0]
    test_label_path = args[1]
    num_pair = int(args[2])

    if len(args) > 3:
        gene_pair_path = args[3]
        if gene_pair_path == '':
            gene_pair_path = './result/training_results_default.csv'
    else:
        gene_pair_path = './result/training_results_default.csv'

    if len(args) > 4:
        separator = args[4]
    else:
        separator = '\t'
    
    
    # load gene paris and tags
    gene_pairs_res = pd.read_table(gene_pair_path,sep='\t')
    gene_pairs = list(gene_pairs_res['pairs'])[0:num_pair]
    gene_pairs = utils.convert_pair_str_2_list(gene_pairs)
    print('%d gene pairs used:' % num_pair)
    print(gene_pairs)
    pair_tags = list(gene_pairs_res['tag'])[0:num_pair]

    # load test data
    test_data = pd.read_table(test_data_path,sep=separator,index_col=[0])
    test_data = utils.profile_preprocessing(test_data)
    print("Number of genes in the test data: %s" % test_data.shape[1])
    print("Number of samples in the test data: %s" % test_data.shape[0])
    test_label_df = pd.read_table(test_label_path)
    test_label = np.array(list(test_label_df['Label']))
    drop_gene = utils.get_test_drop_gene(test_data,gene_pairs)
    test_auc = utils.test_bulk_sc(drop_gene,test_data,test_label,gene_pairs,pair_tags)
    print('Test AUC is %.4f' % test_auc)
    return test_auc


if __name__ == "__main__":
   main()

