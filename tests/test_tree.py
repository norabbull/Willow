# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 08:53:43 2021

@author: norab
"""

import pandas as pd
import numpy as np
import pytest
import unittest

from treeAnalysis.treeInformation import treeInfo
from treeAnalysis.dispersionMetrics import treeDispersion

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
        assert treeInfo.getGeneName('E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv') == "ENSG00000000938___FGR"
        assert treeInfo.getGeneName('E:/Master/cophenetic_dists/ENSG00000188324___OR6C6___CopD.csv') == "ENSG00000188324___OR6C6"
    
    
    
    def test_setGroupTypes():
        
        path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
        types = treeInfo.getGroupTypes(path)
        
        sup = types[0]
        sub = types[1]
        
        assert sup == {'EUR', 'EAS', 'SAS', 'AFR', 'AMR'}
        assert sub == {'MXL', 'PUR', 'TSI', 'PEL', 'PJL', 'MSL', 
                       'CHB', 'ASW', 'ESN', 'STU', 'IBS', 'BEB', 
                       'ACB', 'YRI', 'ITU', 'GWD', 'CHS', 'CDX', 
                       'GBR', 'KHV', 'GIH', 'FIN', 'LWK', 'JPT', 
                       'CLM', 'CEU'}
    
    def test_setSampleInfo():
        
        # Create psuedo matrix
        samples = ['EUR___GBR___HG00096', 'EUR___GBR___HG00099', 
                   'EAS___CHB___NA18642','AFR___ASW___NA20332',
                   'AFR___ASW___NA20340', 'AMR___CLM___HG01384']
        
        cd_mat = pd.DataFrame(np.random.randint(0,10,size=(6, 6)), 
                              columns=samples)
        cd_mat.index=samples
        # Test
        info = treeInfo.getSampleInfo(cd_mat)
        
        assert isinstance(info, dict)
        assert isinstance(info[2], list)
        assert info[0] == ['EUR___GBR___HG00096', 'EUR', 'GBR']
        assert info[5] == ['AMR___CLM___HG01384', 'AMR', 'CLM']
        assert info[4][1] == 'AFR'
        
    def test_getSampleInfo():
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
    
    # Ha en setup hvor du initierer et tree eller en skog med test-data! 

    #%% 
    
    def test_calcPopDists():
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeInfo()
        tree.setup(test_gene_small, pop_info)
        dist_mat = tree.getDistMat()
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



#%% 

    def test_calcMeanTypeDists():
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeInfo()
        tree.setup(test_gene_small, pop_info)
        mat = tree.getDistMat()
        #dists, counts = tree.calcMeanPopDists()
        
        supWith_count = 10 + 6  # AFR and EUR
        subWith_count = 3 + 6   # LWK and GBR
        
        #assert supWith_count == counts['supWith']
        #assert subWith_count == counts['subWith']
        
        # TO DO when tired: 
        # manually add dists from dists mat to check if corresponds to
        # values coming from calcMeanPopDists().
        
    def test_calcMeanPopDists():
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'    
        tree = treeInfo()
        tree.setup(test_gene_small, pop_info)
        tree.calcMeanPopDists()
        test = tree.mean_pop_dists
        pop_dists = tree.pop_dists


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

    