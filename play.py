# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 13:42:43 2021

@author: norab
"""

import pandas as pd
import numpy as np
import os
from os.path import join, isfile
from geneTree.treeInformation import treeInfo
from geneTree.treeMetrics import treeMetrics
from datetime import datetime


test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'

test_files = {'pop_info': pop_info,
           'dist_mat':test_gene_small}

tree = treeMetrics()
tree.setup(test_files['dist_mat'], test_files['pop_info'])
dist_mat = tree.getDistMat()
sample_info = tree.getSampleInfo()
# Understand iteration of upper matrix
dist_mat = dist_mat.to_numpy()

