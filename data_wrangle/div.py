# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 06:45:58 2021

@author: norab
"""




# =============================================================================
#  Make dataframe subset
# =============================================================================
    
test_genes = ['E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000134183___GNAT2___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000105669___COPE___CopD.csv']

pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
             
test_files1 = {'pop_info': pop_info,
              'dist_mat': test_genes[0]}

tree1 = treeInfo(test_files1['dist_mat'], test_files1['pop_info'])
mat = tree1.getDistMat()

submat = mat.iloc[0:10,0:10]
submat.to_csv('C:/Users/norab/MasterDisaster/Data/real_tree_data/dist_mat_test/FGR_10x10.csv')

# Check if equal: 
submat2 = pd.read_csv('C:/Users/norab/MasterDisaster/Data/real_tree_data/dist_mat_test/FGR_10x10.csv', index_col = 0)
(submat == submat2).eq(True).all().all()
