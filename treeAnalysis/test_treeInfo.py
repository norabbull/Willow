# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 08:53:43 2021

@author: norab
"""

import pandas as pd
import numpy as np
import pytest
import unittest


class testTreeDataFunctions(unittest.TestCase):

    def setup_module():
        global group_types
        group_types = getTreeInfo()
        
    
    
        
    def test_setGeneName():
        
        assert getGeneName('E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv') == "ENSG00000000938___FGR"
        assert getGeneName('E:/Master/cophenetic_dists/ENSG00000188324___OR6C6___CopD.csv') == "ENSG00000188324___OR6C6"
    
    
    
    def test_setGroupTypes():
        
        path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
        types = getGroupTypes(path)
        
        sup = types[0]
        sub = types[1]
        
        assert sup == {'EUR', 'EAS', 'SAS', 'AFR', 'AMR'}
        assert sub == {'MXL', 'PUR', 'TSI', 'PEL', 'PJL', 'MSL', 
                       'CHB', 'ASW', 'ESN', 'STU', 'IBS', 'BEB', 
                       'ACB', 'YRI', 'ITU', 'GWD', 'CHS', 'CDX', 
                       'GBR', 'KHV', 'GIH', 'FIN', 'LWK', 'JPT', 
                       'CLM', 'CEU'}
    
    def test_setSampleInfo():
        
        # Create psuedo matrix
        samples = ['EUR___GBR___HG00096', 'EUR___GBR___HG00099', 
                   'EAS___CHB___NA18642','AFR___ASW___NA20332',
                   'AFR___ASW___NA20340', 'AMR___CLM___HG01384']
        
        cd_mat = pd.DataFrame(np.random.randint(0,10,size=(6, 6)), 
                              columns=samples)
        cd_mat.index=samples
        # Test
        info = getSampleInfo(cd_mat)
        
        assert isinstance(info, dict)
        assert isinstance(info[2], list)
        assert info[0] == ['EUR___GBR___HG00096', 'EUR', 'GBR']
        assert info[5] == ['AMR___CLM___HG01384', 'AMR', 'CLM']
        assert info[4][1] == 'AFR'
        
    def test_getSampleInfo():
        
        
        supName1 = info[col][1]
        supName2 = info[row_i][1]
        subName1 = info[col][2]
        subName2 = infop[row_i][2]
    
        return info
    
    # Ha en setup hvor du initierer et tree eller en skog med test-data! 
    # Lag test- data
    #%% 
    
    def test_calcMeanGroupDists():
        pass
        


def calcMeanGroupDists(cd_mat, 
                       sample_info, 
                       group_types):
    """
    Calculate mean inter and intra population distances
    for both classification types (Super and sub).
    These are merged into one function for efficiency reasons.
    
    cd_mat = matrix containing cophenetic distances
    sample_info = indices of sample names 
    group_types = list of sets with population names. 
                Index 0 = super, index 1 = sub 

    """
    cd_mat = cd_mat.to_numpy()
    
    # Create distance sum placeholders
    supWithSums = {group : [0,0] for group in group_types[0]}
    supBetSums = {group : [0,0] for group in group_types[0]}
    subWithSums = {group : [0,0] for group in group_types[1]}   
    subBetSums = {group : [0,0] for group in group_types[1]}
    
    # Iter upper triangular cd_mat 
    for col in range(len(cd_mat)-1):
        i = 1
        for row in range(len(cd_mat)-1):
            row_i = row + i
            
            # Get groupName
            supName1 = sample_info[col][1]
            supName2 = sample_info[row_i][1]
            subName1 = sample_info[col][2]
            subName2 = sample_info[row_i][2]
            
            val = cd_mat[row_i][col]        # Distance value
# =============================================================================
#             Within
# =============================================================================
            if subName1 == subName2:        # Same super pop, same sub pop
                # Super
                supWithSums[supName1] = [i+j for i, j in zip(supWithSums[supName1], [val, 1])]  
                # Sub
                subWithSums[subName1] = [i+j for i, j in zip(subWithSums[subName1], [val, 1])] 
# =============================================================================
#             Between and within
# =============================================================================
            elif supName1 == supName2:      # Same super pop, different sub-pop
                # Super within
                supWithSums[supName1] = [i+j for i, j in zip(supWithSums[supName1], [val, 1])]
                
                # Sub between
                subBetSums[subName1] = [i+j for i, j in zip(subBetSums[subName1], [val, 1])]
                subBetSums[subName2] = [i+j for i, j in zip(subBetSums[subName2], [val, 1])] 
# =============================================================================
#             Between
# =============================================================================
            else:                           # All different: add all to between - groups            
                # Super between
                supBetSums[supName1] = [i+j for i, j in zip(supBetSums[supName1], [val, 1])]
                supBetSums[supName2] = [i+j for i, j in zip(supBetSums[supName2], [val, 1])]
                
                # Sub between
                subBetSums[subName1] = [i+j for i, j in zip(subBetSums[subName1], [val, 1])]
                subBetSums[subName2] = [i+j for i, j in zip(subBetSums[subName2], [val, 1])]                 
        i += 1
        
    summary = {'supBet':0, 'supWith':0, 'subBet':0, 'subWith':0}
    
    # Calculate means 
    for ind, sum_dict in enumerate([supBetSums, supWithSums, subBetSums, subWithSums]):
        sum_df = pd.DataFrame.from_dict(sum_dict, orient = "index", columns = ['total_distance', 'counts'])        
        sum_df['group_mean'] = sum_df.apply(lambda row: row.total_distance / row.counts if row.counts > 0 else 0, axis = 1)
        summary[list(summary.keys())[ind]] = sum_df

    return summary
#%% 

if __name__ == "__main__":
    
    test_getGeneName()
    test_getGroupTypes()

    