# -*- coding: utf-8 -*-
"""
Created on Sat May 22 10:59:48 2021
@author: norab
"""

import re
import pandas as pd
from treeInformation.getTreeInfo import *

# def load_groupTypeSets(path):
#     """
#     Input: Filepath to population informtaion
#     Output: Sets 

#     """
#     return createGroupTypeSets(path)
    
class forestInfo:
    
    def __init__(self):

        # SDRs
        self.SDRs = None
        self.groupSDRs = None
        
        # SDVs
        self.SDVs = None
                
        # Other info
        self.uniq_seq = None
        self.nz_dist_frac = None
              
        
    def load_UniqSeqs(self, path):
        """
        Input: None
        Funtion: Load data. Path is specified in function for functoinality.
        Retun: Pd dataframe with values and gene names as indices.
        """ 
        
        uniqseq = pd.read_csv(path, "r", delimiter = ":", header = 0) 
        uniqseq.rename(index=lambda s: re.sub('_HUMAN__uniq.*', '', s), 
                       columns = {uniqseq.columns[0]: "uniqseq"}, inplace = True)
        self.uniq_seq = uniqseq
    #uniq_seq_dir = "C:/Users/norab/MasterDisaster/Data/meta_data/9381_uniqseqs.txt"
    
    def load_NonZeroDistFraction(self, path):
        """
        Input: path to file containing information on total fraction of
        sample comparisons in each tree that is non-zero. 
        Function: Load data from file. 
        """
        
        nz_dist_frac = pd.read_csv(path, "r", delimiter = ",", index_col = 0)
        nz_dist_frac.rename(index=lambda s: re.sub('_HUMAN__full.*', '', s), 
                        columns = {nz_dist_frac.columns[0]: "nz_dist_frac"}, inplace = True)
        nz_dist_frac.rename(index=lambda s: re.sub('C.*trees/', '', s), 
                        columns = {nz_dist_frac.columns[0]: "nz_dist_frac"}, inplace = True)
        
        self.nz_dist_frac = nz_dist_frac
    #totdist_dir = "C:/Users/norab/MasterDisaster/Data/meta_data/totalDistancesRefined.txt"
    
    def load_GroupNonZeroDistFraction(self, path): 
        self.nz_dist_frac = pd.read_csv(path)
    
    def load_ForestSDRs(self, path):
        """ 
        Input: None
        Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
        
        """
        raw_SDRs = pd.read_csv(path, header=0, index_col=0)
        raw_SDRs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
        self.SDRs = raw_SDRs
  

    def load_ForestSDVs(self, path):
        """
        Input: None
        Return: pd DataFrame containing SDV values for super and sub popultions for each tree/gene
        """
        raw_SDVs = pd.read_csv(path, header=0, index_col=0)
        raw_SDVs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
        self.SDVs = raw_SDVs
    
    def load_ForestGroupSDRs(self, path): 
        self.groupSDRs = pd.read_csv(path)
        
    
    def load_AllForestInfo(self, path):
        # Plan: Nake this return infor already in object instead of reading
        # from file. 
        
        return pd.read_csv(path, index_col=0)
    
    def getSuperSDRs(self): return self.SDRs.loc[:,'SDR_super']
    
    def getSubSDRs(self): return self.SDRs.loc[:,'SDR_sub']
    
    def getUniqSeq(self): return self.uniq_seq
    
    def getNonZeroDistFrac(self): return self.nz_dist_frac()
    
    def getSDVs(self): return self.SDVs
    
    def getAllData(self, path): return self.loadAllData(path)
    
    
if __name__ == "__main__":
    
    SDRs_path = 'E:\Master\SDR\SDR_values_all.csv'
    test = forestInfo()
    test = test.loadForestSDRs(path = SDRs_path)

        