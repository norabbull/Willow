# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 13:13:50 2021

@author: norab
"""

import unittest
import pytest
from treeAnalysis.run_calculations import RunStuff 

class TestMethods(unittest.TestCase):
    
    def test_make_filelist(self):
        
        psuedo_dir = sys.path
        file_list = ['ENSG00000000419___DPM1___CopD.csv', 
                     'ENSG00000000938___FGR___CopD.csv', 
                     'ENSG00000000971___CFAH___CopD.csv', 
                     'ENSG00000001036___FUCO2___CopD.csv']
        correct = ['C:/shockwave/supernova/ENSG00000000419___DPM1___CopD.csv', 
                     'C:/shockwave/supernova/ENSG00000000938___FGR___CopD.csv', 
                     'C:/shockwave/supernova/ENSG00000000971___CFAH___CopD.csv', 
                     'C:/shockwave/supernova/ENSG00000001036___FUCO2___CopD.csv']
        
    
        self.assertRaises(FileNotFoundError, make_filelist(psuedo_dir))
            
        
        result == correct
        
    

runTests.test_make_filelist()


