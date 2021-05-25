# -*- coding: utf-8 -*-
"""
Created on Fri May 21 08:32:27 2021

@author: norab
"""

import pandas as pd
from SDR_SDV import *
from treeInformation.loadTreeData import *
import re
from os import listdir
from os.path import isfile, join
import csv   
#from SDR_SDV import calculateSDRandSDV


def getGeneName(cd_mat_path):
    """
    Input: 
        cd_mat_path: path to cd-file. 
    Function: 
        Filters out treename (usually a gene name). 
    Return: 
        treename as string
    """
    subName = re.sub('^.:.*/ENS', 'ENS', cd_mat_path)
    return re.sub('___CopD.csv$','', subName)

# name = getGeneName(file)

def createGroupTypeSets(group_info_path):
    """
    group_info_path:   pandas dataframe 
    return:     list of sets of included classification names in 
                super- and sub pops respectivly
    """
    group_info = pd.read_csv(group_info_path, delimiter='\t')
    super_pops = group_info.loc[group_info['ClassificationType'] == 'SUPER']
    super_pops = set(super_pops['ClassificationName'])
    sub_pops = group_info.loc[group_info['ClassificationType'] == 'SUB']
    sub_pops = set(sub_pops['ClassificationName'])
    
    return [super_pops, sub_pops]


#group_info_path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
#pop_type_sets = createGroupTypeSets(group_info_path)

def getSampleInfo(cd_mat):
    """
    Input:
    cd_mat: dataframe of cophenetic distances between all pairs in tree
    
    Return: 
    sample_info: dict with created int as index and sample info stored
    in list as value: [gene name, super pop name, sub pop name]
    """
    sample_info = {}
    
    for ind, sample in enumerate(cd_mat.columns):
        subName = re.sub('^...___', '', sample)
        subName = re.sub('___.*$','', subName)
        superName = re.sub('_.*$','', sample)
        
        sample_info[ind] = [sample, superName, subName]
    
    return sample_info

#%% 
def calcMeanGroupDists(cd_mat, 
                       sample_info, 
                       group_names):
    """
    Calculate mean inter and intra population distances
    for both classification types (Super and sub).
    These are merged into one function for efficiency reasons.
    
    cd_mat = matrix containing cophenetic distances
    sample_info = indices of sample names 
    group_names = list of sets with population names. 
                Index 0 = super, index 1 = sub 

    """
    cd_mat = cd_mat.to_numpy()
    
    # Create distance sum placeholders
    supWithSums = {group : [0,0] for group in group_names[0]}
    supBetSums = {group : [0,0] for group in group_names[0]}
    subWithSums = {group : [0,0] for group in group_names[1]}   
    subBetSums = {group : [0,0] for group in group_names[1]}
    
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


def calcNonZeroTotal(cd_mat, 
                     percent = True):
    """
    Input:
        cd_mat : DataFrame with cophenetic distances for a tree

    Returns:
        totDist: float. total amount of pairdistances being non-zero as 
        percent or count. 
        
    Note: This does the same as what u did in R. Have all values alread. 
    """
    nonZero_row = pd.DataFrame((cd_mat != 0).astype(int).sum(axis=1))
    nonZeroCount = int(nonZero_row.sum(axis=0))
    
    num_entries = (cd_mat.shape[0] * cd_mat.shape[1]) - cd_mat.shape[0]
    nonZeroPercent = nonZeroCount / num_entries
    if percent:    
        return nonZeroPercent
    else: 
        return nonZeroCount
    

def calcNonZerosForSamples(cd_mat, 
                           percent = True):
    """
    Input:
        cd_mat : DataFrame with cophenetic distances for a tree

    Returns:
        totDist: dataframe with total amount of pairdistances being non-zero as 
        percent and count per row. 
    """
     
    nonZero_row = pd.DataFrame((cd_mat != 0).astype(int).sum(axis=1), columns = ['nonZero_count'])
    nonZero_row['percent'] = round(nonZero_row['nonZero_count'] / (nonZero_row.shape[0] - 1), 4)
    
    return nonZero_row

# test = calcNonZerosForSamples(MPP5_cd, percent=False)
# test2 = calcNonZeroTotal(MPP5_cd,)

def calcNonZerosForGroups(cd_mat_path,
                         popTypeSet, 
                         classificationType='all', 
                         percent = True):
    """
    Input: 
        cd_mat : DataFrame with cophenetic distances for a tree
        sample_info: a dict with numeric index values as keys and sample
        infor as value, stored on format [gene_name, superpop, subpop]
        popTypeSet: Defines populations in classificationType. 
    Function:
        COunt number of 
    
    Returns:
        Amount of non-zero distances for each subgroup
        
    TO DO: 
        create popTypeSet and groups with functions within this function.
        Have to fix imports and project structure first. 
    """
    cd_mat = load_cd_mat(cd_mat_path)   
    sample_info = getSampleInfo(cd_mat)
    
    # find total sample cells in matrix
    sampleCount = pd.DataFrame(countGroupSamples(sample_info, classificationType=classificationType))
    sampleCount['totalPopComp'] = (sampleCount['sampleCount'] * cd_mat.shape[0]) - sampleCount['sampleCount']    
    
    # Make dict counter for pops
    if classificationType=='super':
        sup_sums = {pop : 0 for pop in popTypeSet[0]}   # Placeholder for sum values
        sup_sample_info = dict([val[0], val[1]] for val in sample_info.values())
        
    elif classificationType=='sub':
        sub_sums = {pop : 0 for pop in popTypeSet[1]}   # Placeholder for sum values
        sub_sample_info = dict([val[0], val[2]] for val in sample_info.values())
        
    elif classificationType=='all':
        sup_sums = {pop : 0 for pop in popTypeSet[0]}   # Placeholder for sum values
        sup_sample_info = dict([val[0], val[1]] for val in sample_info.values())
        
        sub_sums = {pop : 0 for pop in popTypeSet[1]}   # Placeholder for sum values
        sub_sample_info = dict([val[0], val[2]] for val in sample_info.values())
        
    #     sums1 = {pop : 0 for pop in popTypeSet[0]}   # Placeholder for sum values
    #     sums2 = {pop : 0 for pop in popTypeSet[1]}
    #     sums = dict(sums1, **sums2)
        
    #     sample_info1 = dict([val[0], val[1]] for val in sample_info.values())
    #     sample_info2 = dict([val[0], val[2]] for val in sample_info.values())
        
    #     sample_info = {**sample_info1, **sample_info2}
        
    #     return sample_info1
    else: 
        raise ValueError("Invalid calssification type. ")
        
    # map nonZero and sample_info together into matrix
    
    nonZeros = calcNonZerosForSamples(cd_mat, percent = True)  
    if classificationType == 'all':
        nonZeros['supPop'] = nonZeros.index.map(sup_sample_info)
        nonZeros['subPop'] = nonZeros.index.map(sub_sample_info)
        for key, val in nonZeros.iterrows():
            sup_sums[val[2]] += val[0]
            sub_sums[val[3]] += val[0]
    
        sums = {**sup_sums,**sub_sums}
        
    elif classificationType == 'super':
        raise ValueError("This option is not yet available")
    else: 
        raise ValueError("This option is not yet available")
        
    # Divide non-zero values for pop on total pop compariosns    
    sampleCount['popSums'] = sampleCount.index.map(sums)
    sampleCount['percentNonZeroForPop'] = sampleCount.popSums / sampleCount.totalPopComp
    sampleCount['gene'] = getGeneName(cd_mat_path)
    
    return sampleCount


def run_CalcNonZeroForGroup(cd_path, 
                            output_path,
                            group_info_path):
    """
    Input: input path for cd files
    Output: file path to store values. 
        Three columns: gene name, pop and nonZero group value. 
    Return: None

    """
    
    # Population information
    pop_type_sets = createGroupTypeSets(group_info_path)
    
    # List of files
    if isinstance(cd_path, str):     
        file_list = [join(cd_path, f) for f in listdir(cd_path) if isfile(join(cd_path, f))]
    else:            # Continue from disruption
        file_list = cd_path.copy()

    ind = 0
    print("Files to process:\n", file_list)
    
    try:
        while file_list:
            cd_file = file_list.pop().strip()
            print("File processed: ", cd_file)
            print("Number: ", ind)
    
            result = calcNonZerosForGroups(cd_mat_path = cd_file, 
                                popTypeSet=pop_type_sets, 
                                classificationType='all',
                                percent = True)
            with open(output_path, 'a') as f:
                result.to_csv(output_path, mode = 'a', index_label = 'pop', header=f.tell()==0)
            ind += 1
    except: 
        print("Disrupted. Something wrong.. ")

#%% 


def countGroupSamples(sample_info, classificationType = 'super'):
    """
    Count number of smaples for each super or sub population
    classificationType can be 'super' or 'sub' or both
    """
    
    #return sample_info
    if classificationType == 'super':
        result = pd.DataFrame([sup for name,sup,sub in sample_info.values()], 
                              columns = ['sampleCount'])
        return result['sampleCount'].value_counts()
    
    elif classificationType == 'sub':
        result = pd.DataFrame([sub for name,sup,sub in sample_info.values()], 
                              columns = ['sampleCount'])
        return result['sampleCount'].value_counts()
    
    elif classificationType == 'all':
        result1 = pd.DataFrame([sup for name,sup,sub in sample_info.values()], 
                              columns = ['sampleCount'])
        result1 = result1['sampleCount'].value_counts()
        result2 = pd.DataFrame([sub for name,sup,sub in sample_info.values()], 
                              columns = ['sampleCount'])
        result2 = result2['sampleCount'].value_counts()
        result = result1.append(result2)
        
        return result
    
    else: 
        raise ValueError("classificationType is unvalid. Options: super, sub")



#%% Tests
# cd_mat_path = 'E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv'
# cd_mat = load_cd_mat(cd_mat_path)
# sample_info = getSampleInfo(cd_mat)
# test = countGroupSamples(sample_info, classificationType='all')


#%% run


if __name__ == '__main__':
    
# =============================================================================
#     Load data
# =============================================================================

    # cd_mat
    cd_mat_path = 'E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv'
    # pop_
    group_type_sets = createGroupTypeSets('C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv')
    
    # test calcNonZerosForGroups
    test = calcNonZerosForGroups(cd_mat_path, 
                                popTypeSet=group_type_sets, 
                                classificationType='all',
                                percent = True)
    
    # Small subset for testing run_CalcNonZerosForGroup
    cd_path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/'
    output_path = 'C:/Users/norab/MasterDisaster/Data/meta_data/group_totdists_subset.csv'
    group_info_path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
    # True folders for run_CalcNonZerosForGroup
    cd_dir_redhood = 'E:/Master/cophenetic_dists/'
    output_path_redhood = 'E:/Master/otherTreeData/group_totdists_all.csv'
    group_info_path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    

    # test run_calcNonZerosForGroups
    run_CalcNonZeroForGroup(cd_dir_redhood,  # Trenger testdir med noen filer 
                            output_path_redhood,
                            group_info_path)
    



