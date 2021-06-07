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
        self.meanPopDists = None
        self.meanTypeDists = None
        
    
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
        subName = re.sub('^.:.*/ENS', 'ENS', dist_mat_file)
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
    
    def calcPopDists(self):
        """
        Calculate mean inter and intra population distances
        for both classification types (Super and sub).
        These are merged into one function for efficiency reasons.
    
        """
        dist_mat = self.dist_mat.to_numpy()
        
        # Create distance sum placeholders
        # Dict with pop as key and [dist value, number of dist values added] as value
        supWithSums = {pop : [0,0] for pop in self.pop_info[0]}
        supBetSums = {pop : [0,0] for pop in self.pop_info[0]}
        subWithSums = {pop : [0,0] for pop in self.pop_info[1]}   
        subBetSums = {pop : [0,0] for pop in self.pop_info[1]}

        # Iter upper triangular dist_mat
        row_length = len(dist_mat-1)
        row_start = 1
        for col in range(row_length):
            for row in range(row_start, row_length):
                sample1 = col
                sample2 = row
                
                # Get popName
                supPop1 = self.sample_info[sample1][1]
                supPop2 = self.sample_info[sample2][1]
                subPop1 = self.sample_info[sample1][2]
                subPop2 = self.sample_info[sample2][2]
                
                val = dist_mat[sample2][sample1]        # Distance value
    # =============================================================================
    #             Within
    # =============================================================================
                if subPop1 == subPop2:        # Same super pop, same sub pop
                    # Super
                    supWithSums[supPop1] = [i+j for i, j in zip(supWithSums[supPop1], [val, 1])]  
                    # Sub
                    subWithSums[subPop1] = [i+j for i, j in zip(subWithSums[subPop1], [val, 1])] 
    # =============================================================================
    #             Between and within
    # =============================================================================
                elif supPop1 == supPop2:      # Same super pop, different sub-pop
                    # Super within
                    supWithSums[supPop1] = [i+j for i, j in zip(supWithSums[supPop1], [val, 1])]
                    
                    # Sub between
                    subBetSums[subPop1] = [i+j for i, j in zip(subBetSums[subPop1], [val, 1])]
                    subBetSums[subPop2] = [i+j for i, j in zip(subBetSums[subPop2], [val, 1])] 
    # =============================================================================
    #             Between
    # =============================================================================
                else:                           # All different: add all to between - pops            
                    # Super between
                    supBetSums[supPop1] = [i+j for i, j in zip(supBetSums[supPop1], [val, 1])]
                    supBetSums[supPop2] = [i+j for i, j in zip(supBetSums[supPop2], [val, 1])]
                    
                    # Sub between
                    subBetSums[subPop1] = [i+j for i, j in zip(subBetSums[subPop1], [val, 1])]
                    subBetSums[subPop2] = [i+j for i, j in zip(subBetSums[subPop2], [val, 1])]
            row_start += 1
            
        
        
        # TO DO: Write test to check these values are correct for a small matrix
        return {'supWith': supWithSums, 'supBet': supBetSums, 
                'subWith': subWithSums, 'subBet': subBetSums}


        
    def calcMeanPopDists(self):
        """
        Question I repeatedly ask myself: 
            Should I add all and divide by total count or find mean of 
            the individual populations first before adding the means and 
            do a mean of means? 
        """
        pop_dists = self.calcPopDists()
        
        dist_summary = {'supBet':0, 'supWith':0, 'subBet':0, 'subWith':0}
        count_summary = {'supBet':0, 'supWith':0, 'subBet':0, 'subWith':0}
        
        for key, val in pop_dists.items():
            for pop, dist_count in val.items():
                dist_summary[key] += dist_count[0]
                count_summary[key] += dist_count[1]
        
        mean_summary = {key: round(dist_summary[key] / count_summary[key], 8) for key, val
                        in dist_summary.items()}
        
        self.meanPopDists = mean_summary
        
    
    def getMeanPopDists(self): 
        if self.meanPopDists is None: self.calcmeanPopDists()
        return self.meanPopDists
    
    def getSampleInfo(self): 
        if self.sample_info is None: self.setSampleInfo()
        return self.sample_info
    
    def getPopInfo(self): return self.pop_info
    def getGeneName(self): return self.gene_name
    def getDistMat(self): return self.dist_mat 
    
    def calcNonZeroTotal(self, 
                         dist_mat, 
                         percent = True):
        """
        Input:
            dist_mat : DataFrame with cophenetic distances for a tree
    
        Returns:
            totDist: float. total amount of pairdistances being non-zero as 
            percent or count. 
            
        Note: This does the same as what u did in R. Have all values alread. 
        """
        nonZero_row = pd.DataFrame((dist_mat != 0).astype(int).sum(axis=1))
        nonZeroCount = int(nonZero_row.sum(axis=0))
        
        num_entries = (dist_mat.shape[0] * dist_mat.shape[1]) - dist_mat.shape[0]
        nonZeroPercent = nonZeroCount / num_entries
        if percent:    
            return nonZeroPercent
        else: 
            return nonZeroCount
        
    
    def calcNonZerosForSamples(self, 
                               dist_mat, 
                               percent = True):
        """
        Input:
            dist_mat : DataFrame with cophenetic distances for a tree
    
        Returns:
            totDist: dataframe with total amount of pairdistances being non-zero as 
            percent and count per row. 
        """
         
        nonZero_row = pd.DataFrame((dist_mat != 0).astype(int).sum(axis=1), columns = ['nonZero_count'])
        nonZero_row['percent'] = round(nonZero_row['nonZero_count'] / (nonZero_row.shape[0] - 1), 4)
        
        return nonZero_row
    
    # test = calcNonZerosForSamples(MPP5_cd, percent=False)
    # test2 = calcNonZeroTotal(MPP5_cd,)
    
    def calcNonZerosForPops(self, popType='all', percent = True):
        """
        Input: 
            dist_mat : DataFrame with cophenetic distances for a tree
            sample_info: a dict with numeric index values as keys and sample
            infor as value, stored on format [gene_name, superpop, subpop]
            pop_info: Defines populations in popType. 
        Function:
            COunt number of 
        
        Returns:
            Amount of non-zero distances for each subpop
            
        TO DO: 
            create pop_info and pops with functions within this function.
            Have to fix imports and project structure first. 
        """        
        # find total sample cells in matrix
        sampleCount = pd.DataFrame(self.countPopSamples(self.sample_info, popType=popType))
        sampleCount['totalPopComp'] = (sampleCount['sampleCount'] * self.dist_mat.shape[0]) - sampleCount['sampleCount']    
        
        # Make dict counter for pops --- TEST!!!
        if popType=='super':
            sup_sums = {pop : 0 for pop in self.pop_info[0]}   # Placeholder for sum values
            sup_sample_info = dict([val[0], val[1]] for val in self.sample_info.values())
            
        elif popType=='sub':
            sub_sums = {pop : 0 for pop in self.pop_info[1]}   # Placeholder for sum values
            sub_sample_info = dict([val[0], val[2]] for val in self.sample_info.values())
            
        elif popType=='all':
            sup_sums = {pop : 0 for pop in self.pop_info[0]}   # Placeholder for sum values
            sup_sample_info = dict([val[0], val[1]] for val in self.sample_info.values())
            
            sub_sums = {pop : 0 for pop in self.pop_info[1]}   # Placeholder for sum values
            sub_sample_info = dict([val[0], val[2]] for val in self.sample_info.values())
            
        #     sums1 = {pop : 0 for pop in pop_info[0]}   # Placeholder for sum values
        #     sums2 = {pop : 0 for pop in pop_info[1]}
        #     sums = dict(sums1, **sums2)
            
        #     sample_info1 = dict([val[0], val[1]] for val in sample_info.values())
        #     sample_info2 = dict([val[0], val[2]] for val in sample_info.values())
            
        #     sample_info = {**sample_info1, **sample_info2}
            
        #     return sample_info1
        else: 
            raise ValueError("Invalid calssification type. ")
            
        # map nonZero and sample_info together into matrix
        
        nonZeros = self.calcNonZerosForSamples(self.dist_mat, percent = True)  
        if popType == 'all':
            nonZeros['supPop'] = nonZeros.index.map(sup_sample_info)
            nonZeros['subPop'] = nonZeros.index.map(sub_sample_info)
            for key, val in nonZeros.iterrows():
                sup_sums[val[2]] += val[0]
                sub_sums[val[3]] += val[0]
        
            sums = {**sup_sums,**sub_sums}
            
        elif popType == 'super':
            raise ValueError("This option is not yet available")
        else: 
            raise ValueError("This option is not yet available")
            
        # Divide non-zero values for pop on total pop compariosns    
        sampleCount['popSums'] = sampleCount.index.map(sums)
        sampleCount['percentNonZeroForPop'] = sampleCount.popSums / sampleCount.totalPopComp
        sampleCount['gene'] = self.getGeneName()
        
        return sampleCount
    

    #%% But other place. This will instantiate many trees and not be a part of that class. 

    
    #%% 
    
    
    def countPopSamples(sample_info, popType = 'super'):
        """
        Count number of smaples for each super or sub population
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
    



