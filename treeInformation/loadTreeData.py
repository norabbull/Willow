# -*- coding: utf-8 -*-
"""
Created on Sat May 22 10:59:48 2021

@author: norab
"""

import re
import pandas as pd

def load_cd_mat(file):
    """
    Input: 
        Filepath to one cd matrix
    Return: 
        matrix of pairwise distances
    """
    cd = pd.read_csv(file, index_col = 0)
    return cd

def load_uniq_seqs(path):
    """
    Input: None
    Funtion: Load data. Path is specified in function for functoinality.
    Retun: Pd dataframe with values and gene names as indices.
    """ 
    
    uniqseq = pd.read_csv(path, "r", delimiter = ":", header = 0) 
    uniqseq.rename(index=lambda s: re.sub('_HUMAN__uniq.*', '', s), 
                   columns = {uniqseq.columns[0]: "uniqseq"}, inplace = True)
    return uniqseq
#uniq_seq_dir = "C:/Users/norab/MasterDisaster/Data/meta_data/9381_uniqseqs.txt"

def load_tot_dist(path):
    """
    Input: None
    Funtion: Load data. Path is specified in function for functoinality.
    Retun: Pd dataframe with values and gene names as indices.
    """
    
    totdist = pd.read_csv(path, "r", delimiter = ",", index_col = 0)
    totdist.rename(index=lambda s: re.sub('_HUMAN__full.*', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.rename(index=lambda s: re.sub('C.*trees/', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    return totdist
#totdist_dir = "C:/Users/norab/MasterDisaster/Data/meta_data/totalDistancesRefined.txt"

def load_SDRs(path):
    """ 
    Input: None
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    file_path = path
    raw_SDRs = pd.read_csv(file_path, header=0, index_col=0)
    raw_SDRs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
    
    return raw_SDRs
    

def load_SDVs(path):
    """
    Input: None
    Return: pd DataFrame containing SDV values for super and sub popultions for each tree/gene
    """
    
    file_path = path
    raw_SDVs = pd.read_csv(file_path, header=0, index_col=0)
    raw_SDVs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
    return raw_SDVs


def load_GroupwiseSDRs(path):
    
    df = pd.read_csv(path)
    return df
    