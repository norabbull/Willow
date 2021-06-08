# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 13:13:50 2021

@author: norab
"""



import unittest
from datetime import datetime
from treeAnalysis.run_calculations import RunStuff 

#%% 

class TestDispersionMetrics(unittest.TestCase):
   
    def __init__():
        pass
        #RunStuff.__init__(self)
    
        
    def test_make_filelist():
        
        pass
        # psuedo_dir = sys.path
        # file_list = ['ENSG00000000419___DPM1___CopD.csv', 
        #              'ENSG00000000938___FGR___CopD.csv', 
        #              'ENSG00000000971___CFAH___CopD.csv', 
        #              'ENSG00000001036___FUCO2___CopD.csv']
        # correct = ['C:/shockwave/supernova/ENSG00000000419___DPM1___CopD.csv', 
        #              'C:/shockwave/supernova/ENSG00000000938___FGR___CopD.csv', 
        #              'C:/shockwave/supernova/ENSG00000000971___CFAH___CopD.csv', 
        #              'C:/shockwave/supernova/ENSG00000001036___FUCO2___CopD.csv']
        
   
        
    
    def test_runCalcSDR():
        
        input_files = 'E:\\Master\\test_runs\\cophenetic_dists_selection'
        pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
        SDRsuper = 'E:\\Master\\test_runs\\SDRsuper_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
        SDRsub = 'E:\\Master\\test_runs\\SDRsub_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
        save_unprocessed = 'E:\\Master\\test_runs\\unprocessed_files.csv'
        
        run = RunStuff()
        run.run_calcSDR(input_files = input_files, 
                             pop_info = pop_info,
                             SDRsuper_output = SDRsuper,
                             SDRsub_output = SDRsub,
                             unprocessed_files = save_unprocessed)

        
#%% Runs

if __name__ == '__main__':
    
    TestDispersionMetrics.test_runCalcSDR()

      
        
