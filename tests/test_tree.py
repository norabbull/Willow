# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 08:53:43 2021

@author: norab
"""

import pandas as pd
import numpy as np
import pytest
import unittest

from geneTree.treeInformation import treeInfo
from geneTree.treeMetrics import treeMetrics

#%% 

class TestTreeFunctions(unittest.TestCase):
    
    
    def __init__():
        
        pass
        
    def setup_small_tree():
        
        # test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        # pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
        # global test_gene_small
        # global pop_info
        # tree = treeInfo()
        # tree.setup()
        # return tree
        pass

    def test_setGeneName():
        """
        OK.
        """
        test_gene1 = 'E:/Master/cophenetic_dists/ENSG00000188324___OR6C6___CopD.csv'
        test_gene2 = 'E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree1 = treeInfo()
        tree1.setup(test_gene1, pop_info)
        tree2 = treeInfo()
        tree2.setup(test_gene2, pop_info)
        
        assert tree1.getGeneName() == "ENSG00000188324___OR6C6"
        assert tree2.getGeneName() == "ENSG00000000938___FGR"

    
    def test_setPopInfo():
        """
        OK.
        """
        pop_info_file = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
        tree = treeInfo()
        tree.setPopInfo(pop_info_file)
        pop_info = tree.getPopInfo()
        
        sup = pop_info[0]
        sub = pop_info[1]
        
        assert sup == {'EUR', 'EAS', 'SAS', 'AFR', 'AMR'}
        assert sub == {'MXL', 'PUR', 'TSI', 'PEL', 'PJL', 'MSL', 
                       'CHB', 'ASW', 'ESN', 'STU', 'IBS', 'BEB', 
                       'ACB', 'YRI', 'ITU', 'GWD', 'CHS', 'CDX', 
                       'GBR', 'KHV', 'GIH', 'FIN', 'LWK', 'JPT', 
                       'CLM', 'CEU'}
    
    def test_setSampleInfo():
        """
        OK.
        """
        test_gene = 'E:/Master/cophenetic_dists/ENSG00000188324___OR6C6___CopD.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
        tree = treeInfo()
        tree.setup(test_gene, pop_info)
        
        dist_mat = tree.getDistMat()
        sample_info = tree.getSampleInfo()
        
        assert isinstance(sample_info, dict)
        assert isinstance(sample_info[2], list)
        assert sample_info[0] == ['EUR___GBR___HG00261', 'EUR', 'GBR']
        assert sample_info[200] == ['EAS___CHS___HG00409', 'EAS', 'CHS']


    def test_getSampleInfo():
        """
        TO FINISH
        """
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'

        tree = treeInfo()
        tree.setup(test_gene_small, pop_info)
        info = tree.getSampleInfo()
        
        #unfinished
        col = 0
        row_i = 1
        supName1 = info[col][1]
        supName2 = info[row_i][1]
        subName1 = info[col][2]
        subName2 = info[row_i][2]
    
        return info

    #%% 
    
    def test_calcPopDists():
        """
        OK.
        
        """
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeInfo()
        tree.setup(test_gene_small, pop_info)
        pop_dists = tree.calcPopDists() 
        
        AFR_test = pop_dists['supWithSums']['AFR']; AFR_test[0] = round(AFR_test[0], 6)
        EUR_test = pop_dists['supBetSums']['EUR']; EUR_test[0] = round(EUR_test[0], 6)
        LWK_test = pop_dists['subWithSums']['LWK']; LWK_test[0] = round(LWK_test[0], 6)
        CLM_test = pop_dists['subBetSums']['CLM']; CLM_test[0] = round(CLM_test[0], 6)
        
        # Manual calculations retrieved from visual inspection of dist_mat
        supWithSum_AFR = round(0.00236 + 0.00354 + 0.0059 + 0.00472 + 0.00118 + 
                               0.00354 + 0.00236 + 0.00236 + 0.00118 + 0.00118, 6)
        supBetSum_EUR = round(4*0.00236 + 4*0 + 4*(0.00118 + 0.00354 + 2*0.00236), 6)
        subWithSum_LWK = round(0.00236 + 0.0059 + 0.00354, 6)
        subBetSum_CLM = round(0.00472 + 5*0.00236 + 2*0.00118 + 0, 6)

        assert AFR_test == [supWithSum_AFR, 10]  # 5 AFR samples, 10 comparisons
        assert EUR_test == [supBetSum_EUR, 24]   # 4 EUR samples, 24 comparisons
        assert LWK_test == [subWithSum_LWK, 3]   # 3 LWK samples, 3 comparisons
        assert CLM_test == [subBetSum_CLM, 9]    # 1 CLM sample, 9 comparisons 


    def test_calcMeanTypeDists():
        """
        OK.

        """
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeMetrics()
        tree.setup(test_gene_small, pop_info)
        mat = tree.getDistMat()
        pop_dists = tree.getPopDists()
        meanTypeDists = tree.getMeanTypeDists()
        
        # Values read from dist_mat
        subWith = round((0+0.0118)/(6+3), 8)
        subBet = round((0.01416+0.01888+0.01888+0.0472+0.0472) / (9+9+9+24+21), 8) 
        supWith = round((0.02832+0)/(10+6), 8)
        supBet = round((0.0472+0.01888+0.0472) / (25+9+24), 8)
        
        assert supWith == meanTypeDists['supWith']
        assert supBet == meanTypeDists['supBet']
        assert subWith == meanTypeDists['subWith']
        assert subBet == meanTypeDists['subBet']
        

    def test_calcMeanPopDists():
        """
        OK.
        """
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\MIO_5x5.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'    
        tree = treeMetrics()
        tree.setup(test_gene_small, pop_info)
        dist_mat = tree.getDistMat()
        meanPopDists = tree.getMeanPopDists()
        
        # Values read from dist_mat
        supWithAMR = round(0.003 / 1, 6)
        supBetEUR = round((2*0.00108 + 4*0.00054) / 6, 6) 
        subWithGBR = round(0.005 / 1, 6)
        subBetSTU = round((2*0.00108 + 2*0.00054) / 4, 6)
        
        assert supWithAMR == meanPopDists['supWith']['AMR']
        assert supBetEUR == meanPopDists['supBet']['EUR']
        assert subWithGBR == meanPopDists['subWith']['GBR']
        assert subBetSTU == meanPopDists['subBet']['STU']
        

#%% 
    def test_calcSDR():
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\MIO_5x5.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeDispersion()
        tree.setup(test_gene_small, pop_info)
        meanTypeDists = tree.getMeanTypeDists()
        tree.calcSDR()
        treeSubSDR = tree.SDRsub
        treeSupSDR = tree.SDRsuper
        
        # Manual calc: 
        # Mean within superpop dist: 0
        # Mean between superpop dist: 
        mWsup = meanTypeDists['supWith']
        mBsup = meanTypeDists['subBet']
        mWsub = meanTypeDists['subWith']
        mBsub = meanTypeDists['subBet']
        
        subSDR = round(mWsub / mBsub, 6)
        supSDR = round(mWsup / mBsup, 6)        
        
        assert subSDR == treeSubSDR
        assert supSDR == treeSupSDR
        

#%%
    def test_calcSingleSDRs():
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\MIO_5x5.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
        
        tree = treeDispersion()
        tree.setup(test_gene_small, pop_info)
        tree.calcSingleSDRs()
        
        singleSuperSDR = tree.getSingleSuperSDR()
        singleSubSDR = tree.getSingleSubSDR()
        
        treeSubSDR = tree.SDRsub
        treeSupSDR = tree.SDRsuper
        
        # Manual calc: 
        # Mean within superpop dist: 0
        # Mean between superpop dist: 
        mWsup = 0 + 0.005 + 0 + 0 + 0.003
        mBsup = 0.00081 +0.00072 + 0 + 0 + 0.00054
        mWsub = 24*0 + 0.005 + 0.003
        mBsub = 23*0 + 0.00072 + 0.00054 + 0.00081

#%% 

if __name__ == "__main__":
    # Make one instance just to make a samller subset: 
        
    # Specify paths to files
    test_genes = ['E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv',
                 'E:/Master/cophenetic_dists/ENSG00000134183___GNAT2___CopD.csv',
                 'E:/Master/cophenetic_dists/ENSG00000105669___COPE___CopD.csv']
    test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
    pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    

    
# =============================================================================
#     Test claculateions with small dataframe. 
#     Will check correctness by manual calculations
# =============================================================================
    test_files = {'pop_info': pop_info,
                   'dist_mat':test_gene_small}
    
    small_tree = treeInfo()
    small_tree.setup(test_files['dist_mat'], test_files['pop_info'])
    test = small_tree.calcPopDists()
    mat = small_tree.getDistMat()
# =============================================================================
#     Test full set
# =============================================================================
    # Normally, 
    tree1 = treeInfo(test_files1['dist_mat'], test_files1['pop_info'])
    
    mat = tree1.getDistMat()
    
    submat = mat.iloc[0:10,0:10]
    
    submat.to_csv('C:/Users/norab/MasterDisaster/Data/real_tree_data/dist_mat_test/FGR_10x10.csv')
    submat2 = pd.read_csv('C:/Users/norab/MasterDisaster/Data/real_tree_data/dist_mat_test/FGR_10x10.csv', index_col = 0)

    