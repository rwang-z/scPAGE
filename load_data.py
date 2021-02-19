# # coding: utf-8
# Load data from a list of data and label paths

import pandas as pd
import utils

class gene_data():
    def __init__(self,sc_dir_list_raw,label_dir_list_raw,filter_dir_raw):
        self.sc_dir_list_raw = sc_dir_list_raw
        self.filter_dir_raw = filter_dir_raw
        self.label_dir_list_raw = label_dir_list_raw
        self.label = []

    def get_matrix(self,hint_list,separator,label_flag):
        # get filter genes
        if self.filter_dir_raw != None:
            filter_path = self.filter_dir_raw
            filter_table = pd.read_table(filter_path)
            gene_filter = filter_table['Symbol']
            print('The total number of related genes: %d' % len(gene_filter))
        
        # merge training matrices from a list of files in sc_dir_list_raw
        sc_table = []
        sc_data_list = self.sc_dir_list_raw
        label_dir_list = self.label_dir_list_raw
        if label_flag == 'list':
            # if case and control cells are provided in seperate files and labels are given by hint_list
            for sc_data_path,hint in zip(sc_data_list,hint_list):
                print('Processing one training profile...')
                this_df = pd.read_table(sc_data_path,sep=separator,index_col=[0])
                sc_profile = utils.profile_preprocessing(this_df)
                sc_table.append(sc_profile)
                print('Shape of the profile: %d, %d' % (sc_profile.shape[0], sc_profile.shape[1]))
                num_cell = sc_profile.shape[0]
                self.label += [hint] * num_cell
        elif label_flag == 'path':
            # if labels are provided in individual files
            for sc_data_path,label_path in zip(sc_data_list,label_dir_list):
                print('Processing one training profile...')
                this_df = pd.read_table(sc_data_path,sep=separator,index_col=[0])
                sc_profile = utils.profile_preprocessing(this_df)
                sc_table.append(sc_profile)
                print('Shape of the profile: %d, %d' % (sc_profile.shape[0], sc_profile.shape[1]))
                label_df = pd.read_table(label_path)
                self.label += list(label_df['Label'])

        if self.filter_dir_raw != None:
            filtered_table = []
            for table in sc_table:
                gene_sc = list(table.columns)
                overlap_count, overlap_gene = utils.overlap(gene_sc, gene_filter)
                print('The number of filtered genes in the profile: %d' % overlap_count)
                filtered_table.append(table.loc[:,overlap_gene]) 
        else:
            filtered_table = sc_table
        
        if len(filtered_table) >= 2:
            sc_table_all = pd.concat(filtered_table,join='inner')
        elif len(filtered_table) == 1:
            sc_table_all = filtered_table[0]
        print('The number of remained genes in the training data:%d' % sc_table_all.shape[1])
        return sc_table_all

    def gene_label(self):
        label_list = [int(label_) for label_ in self.label]
        all_label = pd.DataFrame(label_list,columns=['labels'])
        return all_label


def train_data_preprocessing(train_dir_list,label_dir_list,hint_list,filter_gene,separator,label_flag):
    train_data = gene_data(train_dir_list,label_dir_list,filter_gene)
    train_matrix = train_data.get_matrix(hint_list,separator,label_flag)
    train_label = train_data.gene_label()  # dataframe, {'labels'}
    return train_matrix, train_label


def load_path(label_flag,file_path_flag):
    if file_path_flag == 'defined':
        train_data_file = './data/train_data_path.txt'
    elif file_path_flag == 'default':
        train_data_file = './data/train_data_path_default.txt'
    train_path_table = pd.read_table(train_data_file)
    train_path_list = list(train_path_table['path'])

    if label_flag == 'path':
        if file_path_flag == 'defined':
            train_label_file = './data/train_label_path.txt'
        elif file_path_flag == 'default':
            train_label_file = './data/train_label_path_default.txt'
        train_label_table = pd.read_table(train_label_file)
        label_path_list = list(train_label_table['path'])
        hint_list = []
    elif label_flag == 'list':
        if file_path_flag == 'defined':
            hint_path = './data/label_list.txt'
        elif file_path_flag == 'default':
            hint_path = './data/label_list_default.txt'
        hint_table = pd.read_table(hint_path)
        hint_list = list(hint_table['hint_list'])
        label_path_list = []
        
    return train_path_list,label_path_list,hint_list



















