# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 13:13:50 2021
@author: norab
"""


import unittest
from datetime import datetime
from geneTree.treeRun import RunStuff 

#%% 

class TestRunStuff(unittest.TestCase):
   
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
        pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
        SDRsuper = 'E:\\Master\\test_runs\\SDRsuper_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
        SDRsub = 'E:\\Master\\test_runs\\SDRsub_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
        save_unprocessed = 'E:\\Master\\test_runs\\unprocessed_files.csv'
        
        run = RunStuff()
        run.run_calcSDR(input_files = input_files, 
                             pop_info = pop_info,
                             SDRsuper_output = SDRsuper,
                             SDRsub_output = SDRsub,
                             unprocessed_files = save_unprocessed)

    def test_runCalcSingleSDR():
        
        input_files = 'E:\\Master\\test_runs\\cophenetic_dists_selection'
        pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
        SSDR_output_dir = 'E:\\Master\\test_runs\\singleSDR_test\\'
        unprocessed_files = 'E:\\Master\\test_runs\\singleSDR_test\\unprocessed_testRun17.06.21_16.00.csv'
        
        run = RunStuff()
        run.run_calcSingleSDRs(
            input_files = input_files, 
            pop_info = pop_info,
            SSDR_output_dir = SSDR_output_dir,
            unprocessed_files = unprocessed_files)
        
        
#%% Runs

if __name__ == '__main__':
    
    #TestRunStuff.test_runCalcSDR()
    
    TestRunStuff.test_runCalcSingleSDR()