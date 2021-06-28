# -*- coding: utf-8 -*-
"""
Created on Fri May 21 08:32:27 2021

@author: norab
"""

import pandas as pd
import re
from random import shuffle

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
        self.psuedo_group_sizes = None
        self.random_pops = False
                
    def setup(self, dist_mat_file, pop_info_file):
        self.setDistMat(dist_mat_file)
        self.setPopInfo(pop_info_file)
        self.setGeneName(dist_mat_file)
        self.setSampleInfo()
        
    def setDistMat(self, dist_mat_file):
        self.dist_mat = pd.read_csv(dist_mat_file, index_col = 0, dtype={'a': str})
        
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
        subName = re.sub('^.*ENS', 'ENS', dist_mat_file)
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
    
    def shuffleSampleInfo(self):
        """
        Make randomly defined population groups. 
        To be used for creatin gnull-distribution.
        
        num_groups = number of psuedo-populations to create. 
        
        How: 
            - Create group identifiers
            - Randomly add group identifier to sample name
            - Set psuedoSuperPops and psuedoSubPops groups. This
            infor will be used instead of sample info in calculations.
            (Probably have to expand other functions to allow this.)
            - Add pseudo-info to sample info!
        """
        
        # Create psuedo group indeifiers               

        # psPops = [f'pseudo{num}/{num_pops}' for num in range(num_pops)]
        # num_samples = int(round(len(self.sample_info) / num_pops, 0))
        # psPops = psPops*num_samples
        
        # for ind in range(len(self.sample_info)):
        #     self.sample_info[ind].append(psPops.pop(randrange(len(psPops))))
        
        # Shuffle all population information randmoly. 
        
        # 1. Get info, stor in list
        # 2. Place sample randomly back into sampleInfo
        
        sup = [val[1] for val in self.sample_info.values()]
        sub = [val[2] for val in self.sample_info.values()]
        shuffle(sup)
        shuffle(sub)

        for key, val in self.sample_info.items():
            val[1] = sup.pop()
            val[2] = sub.pop()
            self.sample_info[key] = val
        
        self.randomPops = True
        
        
               
    def setSampleInfoDesc(self, num_pops): 
        
        self.psuedo_pop_sizes.append(num_pops)
                              
    def getSampleInfo(self): return self.sample_info
    def getPopInfo(self): return self.pop_info
    def getGeneName(self): return self.gene_name
    def getDistMat(self): return self.dist_mat 
    def getSampleInfoDesc(self):
        num_def_pops = len(self.sample_info[0]) - 1
        print(f"There are {num_def_pops} defined.")
        if self.psuedo_pop_sizes():
            print("The psuedo-pop sizes are:  {self.psuedo_pop_sizes}.")
            
    
    def countPopSamples(self, sample_info, popType = 'super'):
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
    


if __name__ == '__main__':
  pass



