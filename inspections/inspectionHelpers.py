# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 09:08:48 2021

@author: norab
"""


import pandas as pd
import re
import os
from os.path import isfile, join
from copy import deepcopy
from random import shuffle, seed


def make_filelist(input_files):

    if isinstance(input_files, str):     
        files = [join(input_files, f) for f in os.listdir(input_files) 
                     if isfile(join(input_files, f))]
    return files

def load_uniq_seqs():
    """
    Input: None
    Funtion: Load data. Path is specified in function for functoinality.
    Retun: Pd dataframe with values and gene names as indices.
    """
    
    uniq_seq_dir = "C:/Users/norab/Master/Data/meta_data/9381_uniqseqs.txt"
    uniqseq = pd.read_csv(uniq_seq_dir, "r", delimiter = ":", header = 0) 
    uniqseq.rename(index=lambda s: re.sub('_HUMAN__uniq.*', '', s), 
                   columns = {uniqseq.columns[0]: "uniqseq"}, inplace = True)
    return uniqseq

def load_tot_dist():
    """
    Input: None
    Funtion: Load data. Path is specified in function for functoinality.
    Retun: Pd dataframe with values and gene names as indices.
    """
    totdist_dir = "C:/Users/norab/Master/Data/meta_data/totalDistancesRefined.txt"
    totdist = pd.read_csv(totdist_dir, "r", delimiter = ",", index_col = 0)
    totdist.rename(index=lambda s: re.sub('_HUMAN__full.*', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.rename(index=lambda s: re.sub('C.*trees/', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    
    return totdist


def load_SDRs(level = "all"):
    """ 
    Input: None
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    # if level == "all":
    #     file_path = 'C:/Users/norab/Master/Data/SDR/SDR_all.csv'
    if level == "sub":
        file_path = 'C:/Users/norab/Master/Data/SDR/SDRsub_all.csv'
    elif level == "super":
        file_path == 'C:/Users/norab/Master/Data/SDR/SDRsuper_all.csv'
    elif level == "psuedo":
        file_path = 'C:/Users/norab/Master/Data/SDRnullDist/nullDistSDRsub_23.06.2021_09.14.csv'
        SDRs = pd.read_csv(file_path)
        SDRs.columns = ['gene', 'psuedoSDR']
    elif level == "all":
        file_path = 'C:/Users/norab/Master/Data/SDRnullDist/nullDistSDRsub_23.06.2021_09.14.csv'# Change
    
    
    
    #raw_SDRs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
    #SDRs.columns = ['level', 'gene', 'SDR']
    
    return SDRs
    

def load_SDVs():
    """
    Input: None
    Return: pd DataFrame containing SDV values for super and sub popultions for each tree/gene
    """
    
    file_path = 'E:\Master\SDV\SDV_values_all.csv'
    raw_SDVs = pd.read_csv(file_path, header=0, index_col=0)
    raw_SDVs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
    
    return raw_SDVs

def load_singleSDR(gene_name, load_all = False):
    """ 
    Input: None
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    if load_all: 
        file_path = deepcopy(gene_name)
        gene_name = re.sub('^.*SSDR_', '', gene_name)
        gene_name = re.sub('.csv', '', gene_name)
    else: 
        file_path = 'C:/Users/norab/Master/Data/singleSDRs/SSDR_' + gene_name + '.csv'

        
    SSDRs = pd.read_csv(file_path, header=0, index_col=0)   
    SSDRs.columns = ['pop', f'SSDR_{gene_name}']
    return SSDRs

def load_allSingleSDRs(num_rand_trees = 50):
    
    file_path = 'C:/Users/norab/Master/Data/singleSDRs/'
    file_list = make_filelist(file_path)
    
    seed(num_rand_trees)
    shuffle(file_list)
    SSDR_all = load_singleSDR(file_list.pop(), load_all = True)
    
    n = 0
    for f in file_list: 
        SSDR_current = load_singleSDR(f, load_all = True)
        SSDR_all = SSDR_all.merge(SSDR_current, on = "pop")
        if n >= num_rand_trees: 
            break
        n +=1

    return SSDR_all
    
    
def load_nullDist(level):
    """ 
    Input: None
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    if level == "sub":
        file_path = "C:/Users/norab/Master/Data/SDRnullDist/nullDistSDRsub_23.06.2021_09.14.csv"
    elif level == "super":
        file_path = "C:/Users/norab/Master/Data/SDRnullDist/nullDistSDRsuper_23.06.2021_09.14.csv"
    elif level == "exit":
        return
    else: 
        print("Invalid option. ")
        level = input("type super, sub or exit: ")
        load_nullDist(level)
    
    vals = pd.read_csv(file_path, names = ['gene', 'psuedoSDR'])
    
    return vals


