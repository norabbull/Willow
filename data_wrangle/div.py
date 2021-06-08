# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 06:45:58 2021

@author: norab
"""

from treeAnalysis.treeInformation import treeInfo
import pandas as pd

test_genes = ['E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000134183___GNAT2___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000105669___COPE___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000164828___SUN1___CopD.csv']

pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    

#%% 

# =============================================================================
#  Make dataframe subset, 10 < 19 dist mat
# =============================================================================
    
test_files1 = {'pop_info': pop_info,
              'dist_mat': test_genes[0]}

tree1 = treeInfo(test_files1['dist_mat'], test_files1['pop_info'])
mat = tree1.getDistMat()

submat = mat.iloc[0:10,0:10]
submat.to_csv('C:/Users/norab/MasterDisaster/Data/real_tree_data/dist_mat_test/FGR_10x10.csv')

# Check if equal: 
submat2 = pd.read_csv('C:/Users/norab/MasterDisaster/Data/real_tree_data/dist_mat_test/FGR_10x10.csv', index_col = 0)
(submat == submat2).eq(True).all().all()


# =============================================================================
# Make dataframe subset, 5 x 5 
# =============================================================================

test_files2 = {'pop_info': pop_info,
              'dist_mat': test_genes[3]}

tree2 = treeInfo()
tree2.setup(test_files2['dist_mat'], test_files2['pop_info'])
mat2 = tree2.getDistMat()

submat2 = mat2.iloc[203:208,203:208]
submat2.to_csv('C:/Users/norab/MasterDisaster/Data/real_tree_data/dist_mat_test/SUN1_5x5.csv')

# Check if equal: 
submat2 = pd.read_csv('C:/Users/norab/MasterDisaster/Data/real_tree_data/dist_mat_test/FGR_10x10.csv', index_col = 0)
(submat == submat2).eq(True).all().all()
