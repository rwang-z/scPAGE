# coding: utf-8
# To validate scGPS on external data

### Details see https://github.com/rwangsunshine/scPAGE

### Usage
# python prediction.py profile_path label_path output_path signature_flag num_pair [separator]

# Input arguments
# profile_path: the path of expression profile for validation. profile matrix: row -- gene, column -- cell, tab-seperated
# label_path: the path of the label file -- a list of labels in column 'Label'.
# output_path: the path of the file to output the predicted scores of the test samples.
# signature_flag: the scGPS used for prediction
#       -- 'example': using the scGPS identified from the example training data (stored in './result/training_results_default.txt')
#       -- 'signature': using the scGPS identified in the paper (stored in './result/pair_list.txt')
#       -- 'user-provided': using the scGPS identified from the data provided by the user (stored in './result/training_results.txt')
# num_pair: the signature size of scGPS used for prediction 
# separator: separator when loading data matrix, default '\t'

### Examples
### Test the scGPS identified from the example data, gene pairs in './result/training_results_default.txt':
# python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'example' 30

### Test the scGPS identified from all single cells in the paper:
# python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'signature' 30

### Test the scGPS identified from the data provided by users:
# python prediction.py './data/GSE119299_exp_matrix.txt' './data/GSE119299_label.txt' './result/GSE119299_pred.txt' 'user-provided' 30 

import sys
import pandas as pd
import numpy as np
import utils

def main():
    ######### input arguments processing #########
    args = sys.argv[1:]
    test_data_path = args[0]
    test_label_path = args[1]
    output_file = args[2]
    signature_flag = args[3]
    num_pair = int(args[4])
    separator = '\t'
    if len(args) > 5:
        separator = args[5]
    
    training_data = signature_flag
    if signature_flag == 'example':
        gene_pair_path = './result/training_results_default.txt'
    elif signature_flag == 'user-provided':
        gene_pair_path = './result/training_results.txt'
    elif signature_flag == 'signature':
        gene_pair_path = './result/pair_list.txt'
        training_data = 'whole single-cell'
        num_pair = 30
    print('The test expression profile: %s' % test_data_path)
    print('The labels of test data: %s' % test_label_path)
    print('Using the scGPS identified from the %s data' % training_data)
    
    # load gene paris and tags
    gene_pairs_res = pd.read_table(gene_pair_path,sep='\t')
    gene_pairs = list(gene_pairs_res['pairs'])[0:num_pair]
    gene_pairs = utils.convert_pair_str_2_list(gene_pairs)
    print('%d gene pairs used for prediction:' % num_pair)
    print(gene_pairs)
    pair_tags = list(gene_pairs_res['tag'])[0:num_pair]

    # load test data
    test_data = pd.read_table(test_data_path,sep=separator, index_col=[0])
    test_data = utils.profile_preprocessing(test_data)  # samples by genes matrix
    samples = test_data.index
    print("Number of genes in the test data: %s" % test_data.shape[1])
    print("Number of samples in the test data: %s" % test_data.shape[0])
    test_label_df = pd.read_table(test_label_path)
    test_label = np.array(list(test_label_df['Label']))
    drop_gene = utils.get_test_drop_gene(test_data,gene_pairs)
    test_auc, pred_score = utils.test_bulk_sc(drop_gene,test_data,test_label,gene_pairs,pair_tags)
    print('Test AUC is %.4f' % test_auc)
    predictions = pd.DataFrame({'sample': samples, 'label': test_label, 'pred_score': pred_score})
    predictions.to_csv(output_file, sep='\t', index=False)
    print('Results saved to %s' % output_file)

if __name__ == "__main__":
   main()

