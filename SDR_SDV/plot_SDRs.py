# -*- coding: utf-8 -*-
"""
Created on Wed May 12 09:20:29 2021

@author: norab
"""

#%% Imports

import pandas as pd
from os import listdir
from os.path import isfile, join
import numpy as np
import re
import csv   
import matplotlib.pyplot as plt


#%% 

def read_SDRs(file_path):
    
    """
    Input: String - path to csv file containing SDR values
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    
    return pd.read_csv(file_path, header=0)

def plt_SDR(data, plot_col1, plot_col2):
    """
    Input: 
            Data: DataFrame - columns with SDR values to plot
            plot_col1: String - column to plot
            plot_col2: String - column to plot
            
    Function: plots scatterplot
            
    """
    
    # y = np.linspace(0, 1.5, 100)
    # ax = plt.gca()
    # data.plot(x =plot_col1, color = 'red', ax=ax)
    # data.plot(x =plot_col2, color = 'blue', ax=ax)
    super_vals = data[plot_col1].sort_values()
    sub_vals = data[plot_col2].sort_values()
    super_vals.plot(use_index = False)
    sub_vals.plot(use_index = False)
    #return super_vals
   
    
    
    #plt.show()

#%% Test plt_SDR

k = plt_SDR(raw_SDRs, 'SDR_super', 'SDR_sub')
    
    
    


#%% 


if __name__ == '__main__':
    
    
    SDR_file_path = 'E:\Master\SDRs\SDR_values_all.csv'
    raw_SDRs = read_SDRs(SDR_file_path)
    
    # Change 0- values to 1- values for more correct representation: 
    #raw_SDRs.loc[raw_SDRs.SDR_super == 0, 'SDR_super'] = 1
    #raw_SDRs.loc[raw_SDRs.SDR_sub == 0, 'SDR_super'] = 1
    # OBS! This is updated in calculateSDR and will not be a porblem next time.
    
    plt_SDR(raw_SDRs,'SDR_super', 'SDR_sub' )
    
