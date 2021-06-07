# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 13:42:43 2021

@author: norab
"""

import pandas as pd
import numpy as np
import pytest
import unittest

from treeAnalysis.treeInformation import treeInfo


test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'

test_files = {'pop_info': pop_info,
           'dist_mat':test_gene_small}

tree = treeInfo()
tree.setup(test_files['dist_mat'], test_files['pop_info'])
dist_mat = tree.getDistMat()
sample_info = tree.getSampleInfo()
# Understand iteration of upper matrix
dist_mat = dist_mat.to_numpy()

count = 0

row_length = len(dist_mat-1)
row_start = 1
for col in range(row_length):
    for row in range(row_start, row_length):
        sample1 = col
        sample2 = row_start
        print("-----------------")
        print("row_start: ", row_start)
        print("Sample 1: ",sample1)
        print("Sample 2: ",sample2)
        # # Get popName
        # supPop1 = sample_info[col][1]
        # print("supPop1:", supPop1)
        # supPop2 = sample_info[row_i][1]
        # print("supPop2", supPop2)
        # subPop1 = sample_info[col][2]
        # print("subPop1", subPop1)
        # subPop2 = sample_info[row_i][2]
        # print("subPop2", subPop2)
        val = dist_mat[sample1][sample2]        # Distance value
        print(val)
        print("----------------------")
        count += 1
        
    row_start += 1


9 + 8 + 7 + 6 + 5 + 4 + 3 + 2 + 1 

row_length = len(dist_mat-1)
for col in range(len(dist_mat)-1):
    row_start = 1
    for row in range(row_start, row_length):
        sample1 = col
        sample2 = row_start
        
        # Get popName
        supPop1 = self.sample_info[sample1][1]
        supPop2 = self.sample_info[sample2][1]
        subPop1 = self.sample_info[sample1][2]
        subPop2 = self.sample_info[sample2][2]
        
        val = dist_mat[sample2][sample1]        # Distance value
                
                
row_start = 3
row_length = 7
for row in range(row_start, row_length):
    #sample1 = 0
    sample2 = row
    #print("s1",sample1)
    print("s2",sample2)
    print("row",row)
    
    
liste = [0.4000005, 5.1, 3.2]
liste[0] = round(liste[0], 4)

def play():
    k = 1
    r = 2
    
    return k, r

i, j = play()
