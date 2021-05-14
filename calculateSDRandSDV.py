# -*- coding: utf-8 -*-
"""
Created on Thu May  6 11:33:03 2021

@author: norab

    - Input: 
        For each tree, a matrix of cophenetic distances between all sample pairs
        is given. A sample is a 3'UTR of an individual human, belonging to a 
        super- and a sub-population group.
        
    - Function:
        
        Read matrices from file directory, process one at a time in a loop.
        Calculate SDR and SDV for each tree/gene/matrix. 
        
        
        #####################################################################
        #                              SDR                                  #
        #####################################################################
        
        SDR: subtype diversity ratio.
        
        SDR =   mean within subtype pairwise distance (mWspd) / 
                mean betweem subtype pairwise distance (mBspd)
        
        2 SDR values are calculated for each tree; 
            - Super-population: 5 populations originating from distinct 
            geographical continents. Subtype = 1 super-population, eg. AFR (Africa)
            - Sub-populations: 26 smaller populations/ethnicities originating from 
            distinct geopgraphical areas. Subtype = 1 sub-population, eg. FIN (Finland)
        
        To calculate SDR for a tree:
                - mWspd for all distances within all defined groups are added 
                and divided by number of comparisons.
                -mBspd for all distances between defined groups are added and
                divided by number of comparisons. 
                
        Since the input matrix of cophenetic distances is triangular symmetric, 
        only the lower triangular matrix is iterated. The diagonal (a sample
        compared to itself) is excluded.
        
        Both within- and between pairwise distances for all groups are 
        calculated in the same loop, to avoid looping through the matrices more than 
        once to save time. The function meanGroupDists performs this task,
        which is given as input to CalculateSDR.
        
        Sample names on the form 'AMR___PUR___HG00553' is given in the input matrix. 
        This is used to extract information of: 
            - Which groups are sub- and super populations, respectivly. 
            Sets of these are created with the function "createPopTypeSets".
            - Since numpy arrays are more efficient in calculation, a dictionary
            matching numerical indices and sample names are created in function 
            "storeSampleInfo", to look up sample info for each sample. 
            
            
        #####################################################################
        #                              SDV                                  #
        #####################################################################
        
        SDV: subtype diversity variance
        
        SDV =   sum(SDRgroup - SDRmean)**2 / n - 1, n = number of samples
        This is calculated with pandas function var().
        
        SDRgroup is the caluclated SDR for a single defined group (eg. AFR or GBR)
        SDRmean is the mean of all SDRs calculated. 
        
        This calculation requires SDR calculations for all groups. 
        SDRgroup =  mean within population pairwise distance /
                    mean between populations pairwise distance
        
        With the option "calc_single_groups = True" in meanGroupDists function, 
        these values are calculated and returned instead of the overall tree SDRs.
        The variance of the SDRs for the super- and subpopulations can thus be
        calculated from the returned list of series conatining the values. 
        
        #####################################################################
        
        All caluclations are pipelined into the function "runSDRorSDVpipe".
        By specifying output path (to export calculated values) to either
        SDR_outputfile or SDV_outputfile (otherwise None), either pipe is
        executed. 
        Directory of distance matrices must be given. 
        
        If calculations are disrupted, a list containing all unprocessed files
        are written to a file, so calculation process can be restored at a later
        point. 
    
    Output: 
        
        SDRs: SDR values are written to 
        
        
        
        
        
"""

#%% imports

import pandas as pd
from os import listdir
from os.path import isfile, join
#import numpy as np
import re
import csv   

#%% Import data

def readCDmatrix(file):
    """
    Input: 
        Filepath to one cd matrix
    Return: 
        matrix of pairwise distances
    """
    cd = pd.read_csv(file, index_col = 0)
    return cd

#%% testing

# file = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/ENSG00000144118___RALB___CopD.csv'
# cd_mat = readCDmatrix(file)

#%% getTreeName

def getTreeName(file_path):
    """
    Input: 
        file_path: path to cd-file. 
    Function: 
        Filters out treename (usually a gene name). 
    Return: 
        treename as string
    """
    subName = re.sub('^.:.*/ENS', 'ENS', file_path)
    return re.sub('___CopD.csv$','', subName)

#%% Test getTreeName

# name = getTreeName(file)
    
#%% createPopTypeSets

"""
Import population / grouping info
"""

def createPopTypeSets(pop_info_path):
    """
    pop_info:   pandas dataframe 
    return:     list of sets of included classification names in 
                super- and sub pops respectivly
    """
    pop_info = pd.read_csv(pop_info_path, delimiter='\t')
    super_pops = pop_info.loc[pop_info['ClassificationType'] == 'SUPER']
    super_pops = set(super_pops['ClassificationName'])
    sub_pops = pop_info.loc[pop_info['ClassificationType'] == 'SUB']
    sub_pops = set(sub_pops['ClassificationName'])
    
    return [super_pops, sub_pops]

#%% Test createPopTypeSets

# pop_info_path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
# pop_type_sets = createPopTypeSets(pop_info_path)


#%% storeSampleInfo

def storeSampleInfo(dist_matrix):
    """
    FILL IN WHEN TIRED
    """
    
    sample_info = {}
    
    for ind, sample in enumerate(dist_matrix.columns):
        subName = re.sub('^...___', '', sample)
        subName = re.sub('___.*$','', subName)
        superName = re.sub('_.*$','', sample)
        
        sample_info[ind] = [sample, superName, subName]
    
    return sample_info

#%% Test storeSampleInfo

# sample_info = storeSampleInfo(cd_mat)

#%% meanGroupDists

def meanGroupDists(cd_mat, cd_inds, pop_names):
    """
    Calculate mean inter and intra population distances
    for both classification types (Super and sub).
    
    cd_mat = matrix containing cophenetic distances
    cd_inds = indices of sample names 
    pop_names = list of sets with population names. 
                Index 1 = super, index 2 = sub 

    """
    
    cd_mat = cd_mat.to_numpy() 
    # Create value placeholders
    subWithSums = {pop : [0,0] for pop in pop_names[1]}   # Placeholder for sum values
    subBetSums = {pop : [0,0] for pop in pop_names[1]}   # Placeholder for sum values
    read_cd_file = {pop : [0,0] for pop in pop_names[0]}   # Placeholder for sum values
    supBetSums = {pop : [0,0] for pop in pop_names[0]}   # Placeholder for sum values
    

    # Approach: only iterat through upper triangular matcd_mat 
    for col in range(len(cd_mat)-1):
        i = 1
        for row in range(len(cd_mat)-1):
            
            # Keep track of iteration of triangular mat
            row_i = row + i
            
            # Get popName
            supName1 = cd_inds[col][1]
            subName1 = cd_inds[col][2]
            
            supName2 = cd_inds[row_i][1]
            subName2 = cd_inds[row_i][2]
            
            val = cd_mat[row_i][col]
# =============================================================================
#             Within
# =============================================================================
            if subName1 == subName2:   # Same super pop, same sub pop
                # Add to both, super and sub will be equal 
                
                # Super
                read_cd_file[supName1] = [i+j for i, j in zip(read_cd_file[supName1], [val, 1])]  # Will already exist, as dict is initialized
                # Sub
                subWithSums[subName1] = [i+j for i, j in zip(subWithSums[subName1], [val, 1])] 
# =============================================================================
#             Between and within
# =============================================================================
            elif supName1 == supName2: # Same super pop, different sub-pop
                # Add to super pop within, add to sub pop between
                # Note: subpops will be unequal as the past test failed. 
                
                # Super within
                read_cd_file[supName1] = [i+j for i, j in zip(read_cd_file[supName1], [val, 1])]
                
                # Sub between (add to both)
                subBetSums[subName1] = [i+j for i, j in zip(subBetSums[subName1], [val, 1])]
                subBetSums[subName2] = [i+j for i, j in zip(subBetSums[subName2], [val, 1])] 
                
# =============================================================================
#             Betwee
# =============================================================================
            else: # All different: add all to between - groups
            
                # Super between
                supBetSums[supName1] = [i+j for i, j in zip(supBetSums[supName1], [val, 1])]
                supBetSums[supName2] = [i+j for i, j in zip(supBetSums[supName2], [val, 1])]
                
                # Sub between
                subBetSums[subName1] = [i+j for i, j in zip(subBetSums[subName1], [val, 1])]
                subBetSums[subName2] = [i+j for i, j in zip(subBetSums[subName2], [val, 1])] 
                
            # Calculate means            
            if row == 20: 
                break
        i += 1
        
    summary = {'supBet':0, 'supWith':0, 'subBet':0, 'subWith':0}
    
    # Calculate means for populations
    for ind, sum_dict in enumerate([supBetSums, read_cd_file, subBetSums, subWithSums]):
        sum_df = pd.DataFrame.from_dict(sum_dict, orient = "index", columns = ['total_distance', 'counts'])        
        sum_df['pop_mean'] = sum_df.apply(lambda row: row.total_distance / row.counts if row.counts > 0 else 0, axis = 1)

        summary[list(summary.keys())[ind]] = sum_df

    return summary


#%%  test meanGroupDists
# mean_dists = meanGroupDists(cd_mat, cd_inds, pop_names)
#%% CalculateSDR
def calculateSDR(total_group_dists, calc_single_groups = False):
    """
    Input: 
        total_group_dists: Numpy arrayList of dataframes containing info of total group distances and
        number of comparisons.
        calc_single_groups: If SDR for all single groups should be calculated. This is 
            required to calculate SDVs. 
    Function: 
        Calculates group mean for either full group (calc_single_groups = False)
            or every single group (.. = True.)
    Return:
        Float: overall mean of distances for classification type        
    """

    
    if calc_single_groups:
        SDR_subPops = total_group_dists['subWith']['pop_mean'].divide(total_group_dists['subBet']['pop_mean'])
        SDR_supPops = total_group_dists['supWith']['pop_mean'].divide(total_group_dists['supBet']['pop_mean'])
        
        SDR_subPops = SDR_subPops.fillna(0)
        SDR_supPops = SDR_supPops.fillna(0)
        
        SDR_subPops = SDR_subPops.replace(0,1)
        SDR_supPops = SDR_supPops.replace(0,1)
        
        SDR_subPops = SDR_subPops.rename('SDR')
        SDR_supPops = SDR_supPops.rename('SDR')
        #SDR_supPops.rename(columns = {'pop_mean': 'SDR'}, inplace = True)
        
        return [SDR_supPops, SDR_subPops]

  
    tot_means = {}
    # Calculate for overall group
    for popType, summary_df in total_group_dists.items():
        tot_means[popType] = sum(summary_df['total_distance'])/sum(summary_df['counts'])

    if not (tot_means['supBet'] == 0):
        superSDR = round(tot_means['supWith'] / tot_means['supBet'], 4)
    else: 
        superSDR = 1
    
    if not (tot_means['subBet'] == 0):
        subSDR = round(tot_means['subWith'] / tot_means['subBet'],4)
    else: 
        subSDR = 1
        
    
    
    return [superSDR, subSDR] # Obs, order of these are 

#%% Test calcSDR

test_SDR = calculateSDR(all_means, calc_single_groups=True)
    
#%% SDV 2nd edition

def calculateSDV(SDRs):
    """
    Input: 
        SDR values for all populations
    Return: 
        SDVs for sup and sub pops (ints in list)
    """    
    return [round(classType.var(), 4) for classType in SDRs]    

#%% Test SDV calc

SDVs = calculateSDV(test_SDR)
          
#%% runSDRorSDVpipe


def runSDRorSDVpipe(cd, pop_info, unprocessed_files, SDR_outputfile=None, SDV_outputfile=None):
    
    """
    Input: 
        cd_files: path to files containing cophenetic distance matrices
        pop_info: path to population classes information file
        output_path: path to file SDR values are to be written to. 
        unprocessed_files_path: path to file for saving list of unproccessed file upon disruption
        
    Function: 
        Pipeline of function for SDR calculations. 
        Writes SDR values to file with specified gene name.
        
        If function is disrupted, a list of the remaining files to 
        calculate SDR for from the filelist is written to a file. 
    
    Output: 
        Writes SDR to file.         
    """
    # Population information
    pop_type_sets = createPopTypeSets(pop_info)
    
    # List of files
    if isinstance(cd, str):
        file_list = [join(cd, f) for f in listdir(cd) if isfile(join(cd, f))]
    else: 
        file_list = cd.copy()

    ind = 0     # Only for print
    print(file_list)
    
    try:
        while file_list:
            
            # SDR
            cd_file = file_list.pop()
            
            print("File processed: ", cd_file)
            print("Number: ", ind)

            cd_mat = readCDmatrix(cd_file)              # read cd matrix
            sample_info = storeSampleInfo(cd_mat)       # store sample information
            
            all_means = meanGroupDists(cd_mat,          # calculate all mean dists
                                     sample_info,       #   between pop types
                                     pop_type_sets)
            
            #SDR
            if SDR_outputfile:
                SDRs = calculateSDR(all_means)              # calculate SDR - something strange here
                SDRs.insert(0, getTreeName(cd_file))
                
                with open(SDR_outputfile, 'a', newline='') as f:   # write to file
                    writer = csv.writer(f)
                    writer.writerow(SDRs)

            
            # SDV
            if SDV_outputfile:
                SDRforSDV = calculateSDR(all_means, calc_single_groups=True)
                SDVs = calculateSDV(SDRforSDV)
                SDVs.insert(0,getTreeName(cd_file))

                with open(SDV_outputfile, 'a', newline='') as f:   # write to file
                    writer = csv.writer(f)
                    writer.writerow(SDVs)
            ind += 1       
     
    except: 
        # Write remaining filelist to file
        print("Disrupted.")
        with open(unprocessed_files, 'w', newline='') as f:
               writer = csv.writer(f)
               writer.writerow(file_list)

#%% Run


if __name__ == '__main__': 
    
    #cd_dir_local = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/'
    # Directory of cophenetic distance files
    cd_dir_redhood = 'E:/Master/cophenetic_dists/'
    
    # Directory of population information file
    pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
    # Directory of the file to output SDR values
    SDR_outputfile = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/SDR/SDR_values_all.csv'
    SDR_outputfile_redhood =  'E:/Master/SDR/SDR_values_all.csv'
    
    SDV_outputfile = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/SDV/SDV_values_all.csv'
    SDV_outputfile_redhood =  'E:/Master/SDV/SDV_values_all.csv'
    
    
    # Directory to store list of files that has not been processed because of disruption
    SDR_unprocessed_files = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/SDR/unprocessed_cd_files.csv'
    SDR_unprocessed_files_redhood = 'E:/Master/SDR/unprocessed_cd_files_all.csv'
    
    SDV_unprocessed_files = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/SDV/unprocessed_cd_files.csv'
    SDV_unprocessed_files_redhood = 'E:/Master/SDV/unprocessed_cd_files_all.csv'


    # # If u wanna start form list of unprocessed files, pass this as list of files
    # with open('E:/Master/SDRs/unprocessed_cd_files.csv', newline='') as f:
    #     reader = csv.readlines(f)
    #     cd_dir_redhood_list = list(reader)
        
    # cd_dir_redhood_list = cd_dir_redhood_list[0]
    # ----------------------------------------------------------------------------
    
    
    #runSDRorSDVpipe(cd_dir_local, pop_info, output_file, unprocessed_files)
    runSDRorSDVpipe(cd_dir_redhood, pop_info, SDV_unprocessed_files_redhood, SDV_outputfile = SDV_outputfile_redhood)

    # #---------------------------------------------------------------
    # Create list of all CD files in directory
    path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/'
    all_files = [join(path,f) for f in listdir(path) if isfile(join(path, f))]

    # test only one, else put this in loop
    one_file = all_files[0]
    cd_mat = readCDmatrix(one_file)
    
    # Create sets of population groups
    pop_info_path = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    pop_type_sets = createPopTypeSets(pop_info_path)
    
    # Create dictionary containing indices and sample information
    sample_info = storeSampleInfo(cd_mat)
    
    # Calculate mean distances
    all_means = meanGroupDists(cd_mat, sample_info, pop_type_sets)

    SDRs = calculateSDR(all_means)
    SDVs = calculateSDV(all_means)
    
    SDRs.insert(0, getTreeName(one_file))
    SDVs.insert(0, getTreeName(one_file))
    
    # Read values to file
    SDR_file='C:/Users/norab/MasterDisaster/Data/real_tree_data/SDRs/SDR_values.csv'
    SDV_file='C:/Users/norab/MasterDisaster/Data/real_tree_data/SDVs/SDV_values.csv'

    with open(SDR_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(SDRs)

    with open(SDV_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(SDVs)
    
    
    
    

    


