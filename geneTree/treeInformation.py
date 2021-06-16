# -*- coding: utf-8 -*-
"""
Created on Fri May 21 08:32:27 2021

@author: norab
"""

import pandas as pd
import re
from os import listdir
from os.path import isfile, join
import csv   
#from SDR_SDV import calculateSDRandSDV

class treeInfo:
    

    def __init__(self):
        
        """
        dist_mat = matrix containing cophenetic distances
        sample_info = indices of sample names 
        pop_info = list of sets with population names.Index 0 = super, index 1 = sub 

        """
        self.dist_mat = None
        self.pop_info = None
        self.super_pops = None
        self.sub_pops = None
        self.sample_info = None
        self.gene_name = None
        self.pop_dists = None
        self.mean_pop_dists = None
        self.mean_type_dists = None
    
    def setup(self, dist_mat_file, pop_info_file):
        self.setDistMat(dist_mat_file)
        self.setPopInfo(pop_info_file)
        self.setGeneName(dist_mat_file)
        self.setSampleInfo()
        
    def setDistMat(self, dist_mat_file):
        self.dist_mat = pd.read_csv(dist_mat_file, index_col = 0)
        
    def setGeneName(self, dist_mat_file):
        """
        Input: 
            file: string filepath to distance matrix-file. 
        Function: 
            Filters out name of gene the tree represents and assign to
            class variable "gene_name".
            Both Ensembl and gene name identifiers included on the form: 
                'ENSG00000000938___FGR'
        """
        subName = re.sub('^.:.*ENS', 'ENS', dist_mat_file)
        self.gene_name = re.sub('___CopD.csv$','', subName)
    
    def setPopInfo(self, pop_info_file):
        """
        file: string filepath to pop type info-file.
              pop types = super and sub
        function: organize information into list of two sets 
                of contained populations in super- and sub pops respectivly
        """
        pop_info = pd.read_csv(pop_info_file, delimiter='\t')
        super_pops = pop_info.loc[pop_info['ClassificationType'] == 'SUPER']
        super_pops = set(super_pops['ClassificationName'])
        sub_pops = pop_info.loc[pop_info['ClassificationType'] == 'SUB']
        sub_pops = set(sub_pops['ClassificationName'])
        
        self.pop_info = [super_pops, sub_pops]   # to be replaced/deleted
        self.super_pops = super_pops
        self.sub_pops = sub_pops
    
 
    def setSampleInfo(self):
        """
        Input:
        dist_mat: dataframe of cophenetic distances between all pairs in tree
        
        Function:
        structure sample_info: dict with created int as index and sample info stored
        in list as value: [gene name, super pop name, sub pop name]
        
        Reason for numerical indices: Can use to iterate triangular matrix in 
        calcTypeDists
        """
        self.sample_info = {}
        
        for ind, sample in enumerate(self.dist_mat.columns):
            subName = re.sub('^...___', '', sample)
            subName = re.sub('___.*$','', subName)
            superName = re.sub('_.*$','', sample)
            
            self.sample_info[ind] = [sample, superName, subName]
     
    def getSampleInfo(self): return self.sample_info
    def getPopInfo(self): return self.pop_info
    def getGeneName(self): return self.gene_name
    def getDistMat(self): return self.dist_mat 
  
    
    #%% 
    
    
    def countPopSamples(sample_info, popType = 'super'):
        """
        TO FINISH AND TEST:
        Count number of samples for each super or sub population
        popType can be 'super' or 'sub' or both
        """
        
        #return sample_info
        if popType == 'super':
            result = pd.DataFrame([sup for name,sup,sub in sample_info.values()], 
                                  columns = ['sampleCount'])
            return result['sampleCount'].value_counts()
        
        elif popType == 'sub':
            result = pd.DataFrame([sub for name,sup,sub in sample_info.values()], 
                                  columns = ['sampleCount'])
            return result['sampleCount'].value_counts()
        
        elif popType == 'all':
            result1 = pd.DataFrame([sup for name,sup,sub in sample_info.values()], 
                                  columns = ['sampleCount'])
            result1 = result1['sampleCount'].value_counts()
            result2 = pd.DataFrame([sub for name,sup,sub in sample_info.values()], 
                                  columns = ['sampleCount'])
            result2 = result2['sampleCount'].value_counts()
            result = result1.append(result2)
            
            return result
        
        else: 
            raise ValueError("popType is unvalid. Options: super, sub")
    
    
    
    #%% Tests
    # dist_mat_file = 'E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv'
    # dist_mat = load_dist_mat(dist_mat_file)
    # sample_info = getSampleInfo(dist_mat)
    # test = countPopSamples(sample_info, popType='all')
    

#%% run


if __name__ == '__main__':

    
    pass





#%% 
# =============================================================================
#     Load data
# =============================================================================
    
    test_gene = 'E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv'
    test_files = {'pop_info': 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv',
                  'dist_mat': test_gene}
    # dist_mat
    dist_mat_file = 'E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv'
    # pop_
    pop_info = setPopInfo('C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv')
    
    # test calcNonZerosForPops
    test = calcNonZerosForPops(dist_mat_file, 
                                pop_info=pop_type_sets, 
                                popType='all',
                                percent = True)
    
    # Small subset for testing run_CalcNonZerosForGroup
    cd_file = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/'
    output_file = 'C:/Users/norab/MasterDisaster/Data/meta_data/group_totdists_subset.csv'
    pop_info_file = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
    # True folders for run_CalcNonZerosForGroup
    cd_dir_redhood = 'E:/Master/cophenetic_dists/'
    output_file_redhood = 'E:/Master/otherTreeData/group_totdists_all.csv'
    pop_info_file = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    

    # test run_calcNonZerosForPops
    run_CalcNonZeroForGroup(cd_dir_redhood,  # Trenger testdir med noen filer 
                            output_file_redhood,
                            pop_info_file)
    



