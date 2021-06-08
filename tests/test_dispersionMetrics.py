# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 12:21:00 2021

@author: norab
"""

import pandas as pd
import numpy as np
import pytest
import unittest

from treeAnalysis.treeInformation import treeInfo
from treeAnalysis.dispersionMetrics import treeDispersion

#%% 

class TestDispersionMetrics(unittest.TestCase):
    
    
    def __init__():
        
        pass
    
        
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
    def test_calcPopSDRs():
        test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\MIO_5x5.csv'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeDispersion()
        tree.setup(test_gene_small, pop_info)
        meanPopDists = tree.getMeanPopDists()
        tree.calcSDR()
        treeSubSDR = tree.SDRsub
        treeSupSDR = tree.SDRsuper
        
        # Manual calc: 
        # Mean within superpop dist: 0
        # Mean between superpop dist: 
        mWsup = 0 + 0.005 + 0 + 0 + 0.003
        mBsup = 0.00081 +0.00072 + 0 + 0 + 0.00054
        mWsub = 24*0 + 0.005 + 0.003
        mBsub = 23*0 + 0.00072 + 0.00054 + 0.00081

        
        