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
        once to save time. The function calcMeanGroupDists performs this task,
        which is given as input to CalculateSDR.
        
        Sample names on the form 'AMR___PUR___HG00553' is given in the input matrix. 
        This is used to extract information of: 
            - Which groups are sub- and super populations, respectivly. 
            Sets of these are created with the function "createGroupTypeSets".
            - Since numpy arrays are more efficient in calculation, a dictionary
            matching numerical indices and sample names are created in function 
            "getSampleInfo", to look up sample info for each sample. 
            
            
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
        
        With the option "calc_SDRgroupwise = True" in calcMeanGroupDists function, 
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
        
        
    Abbreviations and variable explainations: 
        
        SDR - subtype diversity ratio
        SDV subtype diversity variance
        supPop - super population
        subPop - sub population
        supBet - super between, refers to distance between different popultaions
        supWith - super within, refers to distance within a population
        cd - cophenetic distance 
        mat - matrix
        SDRgroupwise - refers to calculation of SDR for every defined group respectivly
    
        
"""

from os import listdir
from os.path import isfile, join
import csv   
import treeInformation as ti

#file = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/ENSG00000001167___NFYA___CopD.csv'
#cd_mat = readCDgene(file)


#sample_info = getSampleInfo(cd_mat)


#all_means = calcMeanGroupDists(cd_mat, sample_info, group_names)

def calculateSDR(cd_mat, sample_info, group_names, calc_SDRgroupwise = False):
    """
    Input: 
        group_dists: Numpy arrayList of dataframes containing info of total group distances and
        number of comparisons.
        calc_SDRgroupwise: If SDR for all single groups should be calculated. This is 
            required to calculate SDVs. 
    Function: 
        Calculates group mean for either full group (calc_SDRgroupwise = False)
            or every single group (.. = True.)
    Return:
        Float: overall mean of distances for classification type        
    """

    group_dists = ti.calcMeanGroupDists(cd_mat, sample_info, group_names)
    
    if calc_SDRgroupwise:
        SDR_subPops = group_dists['subWith']['group_mean'].divide(group_dists['subBet']['group_mean'])
        SDR_supPops = group_dists['supWith']['group_mean'].divide(group_dists['supBet']['group_mean'])
        
        SDR_subPops = SDR_subPops.fillna(0)
        SDR_supPops = SDR_supPops.fillna(0)
        
        SDR_subPops = SDR_subPops.replace(0,1)
        SDR_supPops = SDR_supPops.replace(0,1)
        
        SDR_subPops = SDR_subPops.rename('SDR')
        SDR_supPops = SDR_supPops.rename('SDR')
        #SDR_supPops.rename(columns = {'group_mean': 'SDR'}, inplace = True)
        
        return [SDR_supPops, SDR_subPops]

  
    tot_means = {}
    # Calculate for overall group
    for groupType, summary_df in group_dists.items():
        tot_means[groupType] = sum(summary_df['total_distance'])/sum(summary_df['counts'])

    if not (tot_means['supBet'] == 0):
        superSDR = round(tot_means['supWith'] / tot_means['supBet'], 4)
    else: 
        superSDR = 1
    
    if not (tot_means['subBet'] == 0):
        subSDR = round(tot_means['subWith'] / tot_means['subBet'],4)
    else: 
        subSDR = 1
        
    
    
    return [superSDR, subSDR] # Obs, order of these are 

#test_SDR = calculateSDR(all_means, calc_SDRgroupwise=True)

def calculateSDV(SDRs):
    """
    Input: 
        SDR values for all populations
    Return: 
        SDVs for sup and sub pops (ints in list)
    """    
    return [round(classType.var(), 4) for classType in SDRs]    

#SDVs = calculateSDV(test_SDR)
          
def runSDRorSDVpipe(cd, 
                    group_info, 
                    unprocessed_files, 
                    SDR_outputfile=None, 
                    SDV_outputfile=None, 
                    calc_SDRgroupwise = False):
    """
    Input: 
        cd_files: path to files containing cophenetic distance matrices
        group_info: path to population classes information file
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
    pop_type_sets = ti.createGroupTypeSets(group_info)
    
    # List of files
    if isinstance(cd, str):     
        file_list = [join(cd, f) for f in listdir(cd) if isfile(join(cd, f))]
    else:            # Continue from disruption
        file_list = cd.copy()

    ind = 0
    print("Files to process:\n", file_list)
    
    try:
        while file_list:
            cd_file = file_list.pop().strip()
            print("File processed: ", cd_file)
            print("Number: ", ind)
    
            cd_mat = ti.readCDgene(cd_file)                # read cd matrix
            sample_info = ti.getSampleInfo(cd_mat)       # get sample info
            
            # all_means = calcMeanGroupDists(cd_mat,          # calculate all mean dists
            #                          sample_info,       # between pop types
            #                          pop_type_sets)
            #SDR
            if SDR_outputfile:
                SDRs = calculateSDR(cd_mat, 
                                    sample_info, 
                                    pop_type_sets, 
                                    calc_SDRgroupwise = calc_SDRgroupwise)
                
                if calc_SDRgroupwise: 
                    SDRs = SDRs[0].append(SDRs[1]).to_frame()
                    gene_name = ti.getGeneName(cd_file)
                    SDRs = SDRs.assign(gene = [gene_name] * SDRs.shape[0])

                    with open(SDR_outputfile, 'a') as f:
                        SDRs.to_csv(SDR_outputfile, mode = 'a', index_label = 'pop', header=f.tell()==0)                

                else:                      
                    SDRs = SDRs.insert(0, ti.getGeneName(cd_file))

                    with open(SDR_outputfile, 'a', newline='') as f:   # write to file    
                        writer = csv.writer(f)
                        writer.writerow(SDRs)
            
            # SDV
            if SDV_outputfile:
                SDRgroupwise = calculateSDR(cd_mat,
                                         sample_info, 
                                         pop_type_sets, 
                                         calc_SDRgroupwise=calc_SDRgroupwise)
                
                SDVs = calculateSDV(SDRgroupwise)
                SDVs.insert(0,ti.getGeneName(cd_file))

                with open(SDV_outputfile, 'a', newline='') as f:   # write to file
                    writer = csv.writer(f)
                    writer.writerow(SDVs)
            ind += 1       
     
    except: 
        # Write remaining filelist to file
        print("Disrupted.")
        print("Last file processed:", cd_file)
        with open(unprocessed_files, 'w', newline='') as f:
               writer = csv.writer(f)
               writer.writerow(file_list)


#%% Run

if __name__ == '__main__': 
    
# =============================================================================
#     # Load data
# =============================================================================
    
    # Distance matrices:
    # Test - locally saved
    cd_dir_local = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/' 
    # True - saved on redhood drive
    cd_dir_redhood = 'E:/Master/cophenetic_dists/'
    
    # Population information file:
    pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    
    # Output path for SDR and SDV values:
    # Test
    SDR_output_test = 'C:/Users/norab/MasterDisaster/Data/SDR/SDR_values_subset_allgroups.csv'
    SDV_output_test = 'C:/Users/norab/MasterDisaster/Data/SDV/SDV_values_subset.csv'
    # True
    SDR_output_redhood =  'E:/Master/SDR/SDR_values_all.csv'
    SDV_output_redhood =  'E:/Master/SDV/SDV_values_all.csv'
    
    # Output path for group SDRs
    # True
    SDR_output_redhood_allGroups =  'E:/Master/SDR/SDR_values_all_groups.csv'
    
    
    # Output path for logging unprocessed cd-files upon disruption
    # Test SDR and SDV
    SDR_unprocessed_test = 'C:/Users/norab/MasterDisaster/Data/SDR/unprocessed_cd_files_subset.csv'
    SDV_unprocessed_test = 'C:/Users/norab/MasterDisaster/Data/SDV/unprocessed_cd_files_subset.csv'
    # True SDR and SDV
    SDR_unprocessed_redhood = 'E:/Master/SDR/unprocessed_cd_files.csv'
    SDV_unprocessed_redhood = 'E:/Master/SDV/unprocessed_cd_files.csv'
    # True SDR groups
    SDR_unprocessed_redhood_allGroups = 'E:/Master/SDR/unprocessed_cd_files_all_groups.csv'

# =============================================================================
#     # Run
# =============================================================================
    
    # Test SDR and SDV
    runSDRorSDVpipe(cd_dir_local, pop_info, SDR_unprocessed_redhood, SDR_outputfile = SDR_output_redhood_allGroups, calc_SDRgroupwise=True)
    # True SDR and SDV
    runSDRorSDVpipe(cd_dir_redhood, pop_info, SDR_unprocessed_redhood, SDR_outputfile = SDR_output_redhood_allGroups, calc_SDRgroupwise=True)
    
    # True SDR groups
    runSDRorSDVpipe(cd_dir_redhood, pop_info, SDR_unprocessed_redhood_allGroups, SDR_outputfile = SDR_output_redhood_allGroups, calc_SDRgroupwise=True)
    # Real SDV
    runSDRorSDVpipe(cd_dir_redhood, pop_info, SDV_unprocessed_redhood, SDV_outputfile = SDV_output_redhood)

    # Disruption: 
    # Open list of unprocessed files if disrupted
    with open('E:/Master/SDR/unprocessed_cd_files_all_groups.csv', newline='') as f:
        cd_dir_redhood_continue = list(f)[0].split(',')

    # Continue run: 
    runSDRorSDVpipe(cd_dir_redhood_continue, pop_info, SDR_unprocessed_redhood_allGroups, SDR_outputfile = SDR_output_redhood_allGroups, calc_SDRgroupwise=True)

    # TO DO: Change homecooked disruption system to better logging-system. 
    


