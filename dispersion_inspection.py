# -*- coding: utf-8 -*-
"""
Created on Wed May 12 10:05:17 2021

@author: norab
"""

import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Script to find genes by enetring its name from a list of files in a directory


def read_cd_gene(path):
    """
    Input:
        Name: String -  name of gene to find on the form ENSG00000071205___RHG10
        file_dir: Directory with file  of files 
    Function. 
        Read cophenetic distance matrix of a given gene
    Return: 
        
        Distance matrix
    """
    
    cd = pd.read_csv(path, index_col = 0)
    return cd

#%% Test read_cd_gene

directory = 'E:/Master/cophenetic_dists/'
name = 'ENSG00000071205___RHG10___CopD.csv'
full_name = 'E:/Master/cophenetic_dists/ENSG00000086717___PPE1___CopD.csv'

PPE1_cd = read_cd_gene(full_name)
GNAT2_cd= read_cd_gene('E:/Master/cophenetic_dists/ENSG00000134183___GNAT2___CopD.csv')



#%% Load unique sequences and totalDistancesRefined

def load_uniq_seqs():
    """
    Input: None
    Funtion: Load data. Path is specified in function for functoinality.
    Retun: Pd dataframe with values and gene names as indices.
    """
    
    uniq_seq_dir = "C:/Users/norab/MasterDisaster/Data/meta_data/9381_uniqseqs.txt"
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
    totdist_dir = "C:/Users/norab/MasterDisaster/Data/meta_data/totalDistancesRefined.txt"
    totdist = pd.read_csv(totdist_dir, "r", delimiter = ",", index_col = 0)
    totdist.rename(index=lambda s: re.sub('_HUMAN__full.*', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.rename(index=lambda s: re.sub('C.*trees/', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    
    return totdist



def load_SDRs():
    """ 
    Input: None
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    file_path = 'E:\Master\SDR\SDR_values_all.csv'
    raw_SDRs = pd.read_csv(file_path, header=0, index_col=0)
    raw_SDRs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
    
    return raw_SDRs
    

def load_SDVs():
    """
    Input: None
    Return: pd DataFrame containing SDV values for super and sub popultions for each tree/gene
    """
    
    file_path = 'E:\Master\SDV\SDV_values_all.csv'
    raw_SDVs = pd.read_csv(file_path, header=0, index_col=0)
    raw_SDVs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
    
    return raw_SDVs
    

#%% Test loads

uniqseq = load_uniq_seqs()
totdist = load_tot_dist()
SDRs = load_SDRs()
SDVs = load_SDVs()

# Join dfs
result = pd.concat([totdist, uniqseq, SDRs, SDVs], axis=1, join="inner")

# Drop na rows (dont do for analysing different columns)
result.dropna(subset=result.columns, inplace = True)
result.dropna(subset=[result.columns], inplace = True)

#%%


geneNames = totdist.index.values.tolist()
ind = list(range(0,len(geneNames)))
geneNamesDict = dict(zip(ind, geneNames))

uniqseq.rename(index=geneNamesDict)

#%% Plot

result = result.sort_values('uniqseq')
result['log_uniqseq'] = np.log10(result['uniqseq'])

result2 = result.drop(['uniqseq'], axis= 1)

#plt.plot(result2.uniqseq)

# =============================================================================
# Lineplots
# =============================================================================
# Unique seq vs totdist
plt.figure()
result.plot(x='uniqseq', y='totdist', style='o')
plt.xlabel("unique sequences")
plt.ylabel("total non-zero distances in tree")
plt.show()

# Unique seq vs SDR super
 
plt.figure()
result.plot(x='SDR_super', y='uniqseq', style='o', color = "red")
plt.xlabel("SDR super")
plt.ylabel("unique sequences")
plt.title("SDR vs number of unique sequences in gene. No samples removed.")
plt.show()


sns.pairplot(result2)



# All in one plot
result.plot()

# =============================================================================
# Barchart
# =============================================================================

result['uniqseq'].plot.hist(bins=100, figsize=(8,6))
result['SDR_super'].plot.hist(bins=100, figsize=(8,6))

result2.plot.hist(bins=100, alpha = 0.7)


# =============================================================================
# Scatterplots
# =============================================================================

result.plot.scatter(x='uniqueseq', y='SDR_super')
result.plot.hexbin()










#%% Create dict mapping index to gene name

