# # # coding: utf-8
import numpy as np
import statsmodels.stats.multitest as sm
from fisher import pvalue_npy
import fisher
from sklearn import metrics
import time
from tqdm import tqdm
import collections
import pandas as pd

def profile_preprocessing(data):
    # merge duplicate genes
    data = data.T
    # check duplicate genes
    genes = list(data.keys())
    print('Detecting duplicate genes...')
    duplicated_gene = [item for item, count in collections.Counter(genes).items() if count > 1]
    print('Processing duplicate genes')
    for gene in duplicated_gene:
        if pd.notna(gene):
            if gene in data.keys():
                dup_cols = data[gene]
                data = data.drop(gene,axis=1)
                mean_col = dup_cols.mean(axis=1)
                data[gene] = mean_col
    return data

def overlap(a,b):
    overlap_element = list(set(a).intersection(set(b)))
    count = len(overlap_element)
    return count, overlap_element

def pairconvert(data, index):
    sub = np.array(data.iloc[:,index[1,:]]) - np.array(data.iloc[:,index[0,:]]) > 0    
    sub = sub*2-1
    return sub

def tag_assignment(index,i,a,b,c,d):
    a_sub = a[index[0,i],index[1,i]]
    b_sub = b[index[0,i],index[1,i]]
    c_sub = c[index[0,i],index[1,i]]
    d_sub = d[index[0,i],index[1,i]]
    sum_ad = a_sub + d_sub
    sum_bc = b_sub + c_sub
    if sum_ad > sum_bc:
        tag= -1  # 'normal-positive'
    else:
        tag = 1   # 'disease-positive'
    return tag

def train_rank_ratio(data,label,index,a,b,c,d):
    rankdata = pairconvert(data,index)  # -1/+1 for each pair
    num = len(index[0,:])
    tag_list = []  # flag for disease, 1: +, -1: -
    for i in range(0,num):  # for all selected pairs
        tag = tag_assignment(index,i,a,b,c,d)
        tag_list.append(tag)

    tag_list = np.array(tag_list)
    rankdata = rankdata*tag_list > 0  # > 0 --> disease (True, 1); else: health (False, 0)
    num_pair = rankdata.shape[1]
    rankpredt_more = np.squeeze(np.sum(rankdata, axis = 1))/num_pair
    auc = metrics.roc_auc_score(label, rankpredt_more)
    return auc, tag_list

def gene_relative_rank(length,n_sm,x_train):
    t1 = time.time()
    result1 = np.zeros([length,length,n_sm],dtype=np.bool)
    result1_ = np.zeros([length,length,n_sm],dtype=np.bool)
    for i in tqdm(range(0,n_sm)):
        time.sleep(0.05)
        time_l = time.time()
        x = np.tile(x_train.iloc[i][:], (length,1))
        sub = x - x.T
        result1[:,:,i] = (sub > 0)
        result1_[:,:,i] = (sub < 0)
        #print('loop {:d} : {:5f}'.format(i, time.time() - time_l))
    return result1,result1_

def get_value(length,n_sm,x_train,y_train):
    result1,result1_ = gene_relative_rank(length,n_sm,x_train)
    ##Figure out the number of greater or lower rank pairs among all the samples
    con_label = np.squeeze(1 - y_train)
    y_train = np.squeeze(y_train)
    a = np.sum(result1*con_label,axis = 2)
    b = np.sum(result1_*con_label,axis = 2)
    # Case
    c = np.sum(result1*y_train, axis = 2)
    d = np.sum(result1_*y_train, axis = 2)
    return a,b,c,d

def fisher_index_mat(a,b,c,d,num_gene):
    a_ = a.reshape((-1,1))
    a_ = np.squeeze(a_)
    a_ = a_.astype(np.uint)
    b_ = b.reshape((-1,1))
    b_ = np.squeeze(b_)
    b_ = b_.astype(np.uint)
    c_ = c.reshape((-1,1))
    c_ = np.squeeze(c_)
    c_ = c_.astype(np.uint)
    d_ = d.reshape((-1,1))
    d_ = np.squeeze(d_)
    d_ = d_.astype(np.uint)
    # fisher exact test
    _, _, twosided = pvalue_npy(a_, b_, c_, d_)
    rejected, pvalue_Bonf, alphacSidak, alphacBonf = sm.multipletests(twosided, alpha=0.05, method='bonferroni',
                                                                      is_sorted=False, returnsorted=False)
    pvalue_Bonf2 = pvalue_Bonf.reshape((num_gene,num_gene))
    length = len(pvalue_Bonf2[0,:])
    # add 1 in triangle matrix to remove duplicated index
    pvalue_matrix = pvalue_Bonf2 + np.triu(np.ones([length,length]))
    return pvalue_matrix

def get_gene_list(gene_pairs):
    gene_list = []
    for pair in gene_pairs:
        gene_list.append(pair[0])
        gene_list.append(pair[1])
    gene_list = list(set(gene_list))
    return gene_list

def get_test_drop_gene(test_data,gene_pairs):
    gene_list = get_gene_list(gene_pairs)
    test_genes = test_data.columns
    drop_gene = [gene for gene in gene_list if gene not in test_genes]
    if len(drop_gene)!=0:
        print("%d genes in scGPS not detected in the test dataset" % len(drop_gene))
        print("They are: %s" % str(drop_gene))
    else:
        print("All genes in scGPS detected in the test dataset")
    return drop_gene

def convert_pair_str_2_list(pair_str):
    pair_list = []
    for pair in pair_str:
        # print('Pair string: %s' % pair)
        str_list_1 = pair.split('\'')
        gene_1 = str_list_1[1]
        gene_2 = str_list_1[3]
        pair_new = [gene_1, gene_2]
        pair_list.append(pair_new)
        # print('Splited pair:')
        # print(pair_new)
    return pair_list

def pair_2_pair_index(drop_gene,pair_list,tags):
    dele=[]
    pair_index = np.array([[pair[0] for pair in pair_list],[pair[1] for pair in pair_list]])  # reshape to 2 by x
    if len(drop_gene) > 0:
        for gene in drop_gene:
            for i in range(pair_index.shape[1]):
                if gene in pair_index[:,i]:
                    dele.append(i)
        dele = list(set(dele))
        new_index = np.delete(pair_index,dele,axis=1)
        new_tags = np.delete(tags,dele,axis=0)
    else:
        new_index = pair_index
        new_tags = tags
    return new_index, new_tags

def pairconvert_test(data,gene_pair_index):
    # convert expression matrix (data, sample by gene) to 1/-1
    sub = np.array(data.iloc[:][gene_pair_index[1,:]]) - np.array(data.iloc[:][gene_pair_index[0,:]]) > 0  # gj - gi > 0 (gi, gj in gene_pair_index)
    sub = sub * 2 - 1
    return sub # sample by number of gene pairs

def test_rank_ratio(test_data,test_label,pair_index,tag):
    rankdata = pairconvert_test(test_data,pair_index)
    rankdata = rankdata * tag > 0
    num_pair = rankdata.shape[1]
    rankpredt_more = np.squeeze(np.sum(rankdata, axis = 1))/num_pair
    auc = metrics.roc_auc_score(test_label, rankpredt_more)
    return auc, rankpredt_more

def test_bulk_sc(drop_gene,test_data,test_label,gene_pairs,tag):
    pair_index, pair_tags = pair_2_pair_index(drop_gene, gene_pairs, tag)
    test_auc, pred_score = test_rank_ratio(test_data,test_label,pair_index,pair_tags)
    return test_auc, pred_score
