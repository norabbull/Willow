# -*- coding: utf-8 -*-
"""
Created on Wed May 12 10:05:17 2021

@author: norab
"""


import pandas as pd

# Script to find genes by enetring its name from a list of files in a directory


def read_cd_gene(name, file_dir):
    """
    Input:
        Name: String -  name of gene to find on the form ENSG00000071205___RHG10
        file_dir: Directory with file  of files 
    Function. 
        Read cophenetic distance matrix of a given gene
    Return: 
        Distance matrix
    """
    full_path = file_dir+name
    cd = pd.read_csv(full_path, index_col = 0)
    return cd

directory = 'E:/Master/cophenetic_dists/'
name = 'ENSG00000071205___RHG10___CopD.csv'

RHG10_cd = read_cd_gene(name, directory)