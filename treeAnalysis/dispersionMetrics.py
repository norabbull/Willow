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
import pandas as pd
from treeAnalysis.treeInformation import treeInfo
#file = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/ENSG00000001167___NFYA___CopD.csv'
#cd_mat = readCDgene(file)

    
    #sample_info = getSampleInfo(cd_mat)
    
    
    #all_means = calcMeanGroupDists(cd_mat, sample_info, group_types
class treeDispersion(treeInfo):
    
    def __init__(self):
        
        
        self.SDRsuper = None
        self.SDRsub = None
        self.treeSDVs = None
        self.popSDRs = None
        
        treeInfo.__init__(self)
        
        
    def calcSDR(self):
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
        
        # Dette blir feil. Du tar gjennomsnitt av gjennomsnitt. 
        # Må bruke dists direkte, legge sammen alle og dele på total count. 
        if not self.mean_type_dists: self.calcMeanTypeDists()

        if self.mean_type_dists['supBet']:
            self.SDRsuper = round(self.mean_type_dists['supWith'] / 
                                  self.mean_type_dists['supBet'], 6)
        else: 
            self.SDRsuper = float("NaN")
        
        if self.mean_type_dists['subBet']:
            self.SDRsub = round(self.mean_type_dists['subWith'] / 
                                self.mean_type_dists['subBet'], 6)
        else:
            self.SDRsub = float("NaN")
        
    
    #test_SDR = calculateSDR(all_means, calc_SDRgroupwise=True)
    
    def calcPopSDRs(self):
        
        frame = pd.DataFrame(self.mean_pop_dists)
        
        for key, val in self.mean_type_dists:
            for pop, mean in val:
                pass
        
        # SDR_subpops = frame['sub']
        # SDR_subops = self.mean_pop_dists['subWith']['group_mean'].divide(group_dists['subBet']['group_mean'])
       # SDR_supPops = group_dists['supWith']['group_mean'].divide(group_dists['supBet']['group_mean'])
        
        # SDR_subPops = SDR_subPops.fillna(0)
        # SDR_supPops = SDR_supPops.fillna(0)
        
        # SDR_subPops = SDR_subPops.replace(0,1)
        # SDR_supPops = SDR_supPops.replace(0,1)
        
        # SDR_subPops = SDR_subPops.rename('SDR')
        # SDR_supPops = SDR_supPops.rename('SDR')
        # #SDR_supPops.rename(columns = {'group_mean': 'SDR'}, inplace = True)
        
        # self.superPopSDRs = SDR_supPops
        # self.subPopSDRs = SDR_subPops
    
    
    def calculateSDV(self, SDRs):
        """
        Input: 
            SDR values for all populations
        Return: 
            SDVs for sup and sub pops (ints in list)
        """    
        return [round(classType.var(), 4) for classType in SDRs]    
    
    #SDVs = calculateSDV(test_SDR)
    
    def getSDRsuper(self): return self.SDRsuper
    def getSDRsub(self): return self.SDRsub


#%% Run

if __name__ == '__main__': 
    pass
    
# =============================================================================
#     # Load data
# =============================================================================

#%% Testing
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



