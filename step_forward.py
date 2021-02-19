# coding: utf-8
import utils
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

def sfa_iterative(index_ind,k_pair,x_train,y_train,index,a,b,c,d,num_pair):
    max_auc = 0
    max_ind = 0
    s_index = index_ind
    for i in range(0,num_pair):
        if i in s_index:
            continue
        k_index = list(s_index)
        k_index.append(i)
        rankauc, tag_list = utils.train_rank_ratio(x_train,y_train,index[:,k_index],a,b,c,d)  
        if rankauc > max_auc:
            max_ind = i
            max_auc = rankauc
    max_tag = utils.tag_assignment(index,max_ind,a,b,c,d)
    print("The best training result of %d pairs is: %.4f" % (k_pair,max_auc))
    return max_auc,max_ind,max_tag


def get_optimized_pair(train_x,train_y,pvalue,index,num_pair,a,b,c,d,pair_limit,flag,plot=True):
    # index: index of genes in train_x of each DEP
    k_pair_lst = [i for i in range(1, pair_limit + 1)]  # set of number of features to select
    index_ind = []  # index of pairs selected in index matrix
    selected_pairs = []  # gene name of selected pairs
    tag_list = []
    auc_list = []
    for k_pair in k_pair_lst:
        max_auc,new_index_ind,max_tag = sfa_iterative(index_ind,k_pair,train_x,train_y,index,a,b,c,d,num_pair)
        index_ind.append(new_index_ind)
        new_index = index[:,new_index_ind]
        pair_name = [train_x.columns[new_index[0]],train_x.columns[new_index[1]]]
        selected_pairs.append(pair_name)
        tag_list.append(max_tag)
        auc_list.append(max_auc)
    res_df = pd.DataFrame({'pairs': selected_pairs, 'tag': tag_list, 'auc': auc_list})
    if plot:
        x = [i + 1 for i in range(pair_limit)]
        plt.figure(figsize=(8,5))
        plt.plot(x,auc_list,"-",color="red")
        plt.xlabel("Number of pairs",fontsize=15)
        plt.ylabel("AUC",fontsize=15)
        plt.ylim((0,1.03))
        plt.xlim((0,pair_limit + 1))
        if flag == 'defined':
            plt.savefig('./figure/training_auc_test.pdf')
        elif flag == 'default':
            plt.savefig('./figure/training_auc_default.pdf')
            
    return res_df
