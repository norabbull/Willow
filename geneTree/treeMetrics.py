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
import numpy as np
from geneTree.treeInformation import treeInfo
#file = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/cophenetic_dists/ENSG00000001167___NFYA___CopD.csv'
#cd_mat = readCDgene(file)

    
    #sample_info = getSampleInfo(cd_mat)
    
    
    #all_means = calcMeanGroupDists(cd_mat, sample_info, group_types
class treeMetrics(treeInfo):
    
    def __init__(self):
        
        self.pop_dists = None
        self.mean_pop_dists = None
        self.mean_type_dists = None
        self.SDRsuper = None
        self.SDRsub = None
        self.SDVsuper = None
        self.SDVsub = None
        self.singleSDRs = None
        
        treeInfo.__init__(self)
    
    def calcPopDists(self):
        """
        Calculate inter and intra population distances
        for both classification types (Super and sub).
        These are merged into one function for efficiency reasons.
    
        """
        dist_mat = self.dist_mat.to_numpy()
        
        # Create distance sum placeholders
        # Dict with pop as key and [dist value, number of dist values added] as value
        supWithSums = {pop : [0,0] for pop in self.pop_info[0]}
        supBetSums = {pop : [0,0] for pop in self.pop_info[0]}
        subWithSums = {pop : [0,0] for pop in self.pop_info[1]}   
        subBetSums = {pop : [0,0] for pop in self.pop_info[1]}

        # Iter upper triangular dist_mat
        row_length = len(dist_mat-1)
        row_start = 1
        for col in range(row_length):
            for row in range(row_start, row_length):
                sample1 = col
                sample2 = row
                
                # Get popName
                supPop1 = self.sample_info[sample1][1]
                supPop2 = self.sample_info[sample2][1]
                subPop1 = self.sample_info[sample1][2]
                subPop2 = self.sample_info[sample2][2]
                
                val = dist_mat[sample2][sample1]        # Distance value
    # =============================================================================
    #             Within
    # =============================================================================
                if subPop1 == subPop2:        # Same super pop, same sub pop
                    # Super
                    supWithSums[supPop1] = [i+j for i, j in zip(supWithSums[supPop1], [val, 1])]  
                    # Sub
                    subWithSums[subPop1] = [i+j for i, j in zip(subWithSums[subPop1], [val, 1])] 
    # =============================================================================
    #             Within and Between
    # =============================================================================
                elif supPop1 == supPop2:      # Same super pop, different sub-pop
                    # Super within
                    supWithSums[supPop1] = [i+j for i, j in zip(supWithSums[supPop1], [val, 1])]
                    
                    # Sub between
                    subBetSums[subPop1] = [i+j for i, j in zip(subBetSums[subPop1], [val, 1])]
                    subBetSums[subPop2] = [i+j for i, j in zip(subBetSums[subPop2], [val, 1])] 
    # =============================================================================
    #             Between
    # =============================================================================
                else:                           # All different: add all to between - pops            
                    # Super between
                    supBetSums[supPop1] = [i+j for i, j in zip(supBetSums[supPop1], [val, 1])]
                    supBetSums[supPop2] = [i+j for i, j in zip(supBetSums[supPop2], [val, 1])]
                    
                    # Sub between
                    subBetSums[subPop1] = [i+j for i, j in zip(subBetSums[subPop1], [val, 1])]
                    subBetSums[subPop2] = [i+j for i, j in zip(subBetSums[subPop2], [val, 1])]
            row_start += 1
            

        # TO DO: Write test to check these values are correct for a small matrix
        self.pop_dists = {'supWith': supWithSums, 'supBet': supBetSums, 
                'subWith': subWithSums, 'subBet': subBetSums}


        
    def calcMeanTypeDists(self):
        """
        Calculates mean distance between all pairwise samples within a
        classification type. A classification type is either 
        Question I repeatedly ask myself: 
            Should I add all and divide by total count or find mean of 
            the individual populations first before adding the means and 
            do a mean of means? 
        """
        if not self.pop_dists: self.calcPopDists()
        
        dist_summary = {'supBet':0, 'supWith':0, 'subBet':0, 'subWith':0}
        count_summary = {'supBet':0, 'supWith':0, 'subBet':0, 'subWith':0}
        
        for key, val in self.pop_dists.items():
            for pop, dist_count in val.items():
                dist_summary[key] += dist_count[0]
                count_summary[key] += dist_count[1]
        
        self.mean_type_dists = {key: round(dist_summary[key] / count_summary[key], 8) 
                                for key, val in dist_summary.items()}
        
        
        
    def calcMeanPopDists(self):
        """
        Calculates mean distance between pairwise samples within a population,
        for each classification type. That is, for both super and sub populations,
        the mean inter- and intra distances are calculated and stored in 
        mean_pop_dists. 
        
        """
        
        if not self.pop_dists: self.calcPopDists()
        
        mean_pop_dists = {}
        
        for key, val in self.pop_dists.items():
            type_means = {}
            
            for pop, dist_count in val.items():
                if dist_count[1]:  
                    type_means[pop] = round(dist_count[0] / dist_count[1], 6)
                else:   # No distance calculated for this pop
                    type_means[pop] = 0
             
            mean_pop_dists[key] = type_means
             
        self.mean_pop_dists = mean_pop_dists
        
        
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
    
    def calcSingleSDRs(self):
        """
        Ratio of mean intra population distance to mean inter population distance
        for each population. 
        
        Example: Africa
            mean(distances between all African samples) /
            mean(distances between all African - non African samples )
            
        Can be used to calculate SDV. 
        
        Returns
        -------
        List of ratios for each population. 
        
        Do: 
            fetch the right mean distances from mean_pop_dists. 
            Super within
            Super between
            Sub within
            Sub between
        """
        
        self.singleSuperSDR = {}
        self.singleSubSDR = {}
        frame = pd.DataFrame(self.mean_pop_dists)
        super_pop_info = self.getSuperPopInfo() # Set with existing population groups        
        sub_pop_info = self.getSubPopInfo()
        
        for pop in super_pop_info: 
            ratio = self.mean_pop_dists['supWith'][pop] / self.mean_pop_dists['supBet'][pop]
            self.singleSuperSDR[pop] = ratio
        
        for pop in super_pop_info: 
            ratio = self.mean_pop_dists['subWith'][pop] / self.mean_pop_dists['subBet'][pop]
            self.singleSubSDR[pop] = ratio

    
    def calculateSDV(self):
        """
        Measure of variability occuring between inidivdual clusters in a tree. 
        
        
        Input: 
            SDR values for all populations
        Return: 
            SDVs for sup and sub pops (ints in list)
        """
        if not self.singleSuperSDR: 
            self.calcSingleSDRs()
        
        # Test if this works
        self.SDVsuper = round(np.array(self.singleSuperSDR.values()).var, 4)
        self.SDVsub = round(np.array(self.singleSubSDR.values()).var, 4)
    
    def calcNonZeroTotal(self, 
                         dist_mat, 
                         percent = True):
        """
        Input:
            dist_mat : DataFrame with cophenetic distances for a tree
    
        Returns:
            totDist: float. total amount of pairdistances being non-zero as 
            percent or count. 
            
        Note: This does the same as what u did in R. Have all values alread. 
        """
        nonZero_row = pd.DataFrame((dist_mat != 0).astype(int).sum(axis=1))
        nonZeroCount = int(nonZero_row.sum(axis=0))
        
        num_entries = (dist_mat.shape[0] * dist_mat.shape[1]) - dist_mat.shape[0]
        nonZeroPercent = nonZeroCount / num_entries
        if percent:    
            return nonZeroPercent
        else: 
            return nonZeroCount
        
    
    def calcNonZerosForSamples(self, 
                               dist_mat, 
                               percent = True):
        """
        Input:
            dist_mat : DataFrame with cophenetic distances for a tree
    
        Returns:
            totDist: dataframe with total amount of pairdistances being non-zero as 
            percent and count per row. 
        """
         
        nonZero_row = pd.DataFrame((dist_mat != 0).astype(int).sum(axis=1), columns = ['nonZero_count'])
        nonZero_row['percent'] = round(nonZero_row['nonZero_count'] / (nonZero_row.shape[0] - 1), 4)
        
        return nonZero_row
    
    # test = calcNonZerosForSamples(MPP5_cd, percent=False)
    # test2 = calcNonZeroTotal(MPP5_cd,)
    
    def calcNonZerosForPops(self, popType='all', percent = True):
        """
        Input: 
            dist_mat : DataFrame with cophenetic distances for a tree
            sample_info: a dict with numeric index values as keys and sample
            infor as value, stored on format [gene_name, superpop, subpop]
            pop_info: Defines populations in popType. 
        Function:
            COunt number of 
        
        Returns:
            Amount of non-zero distances for each subpop
            
        TO DO: 
            create pop_info and pops with functions within this function.
            Have to fix imports and project structure first. 
        """        
        # find total sample cells in matrix
        sampleCount = pd.DataFrame(self.countPopSamples(self.sample_info, popType=popType))
        sampleCount['totalPopComp'] = (sampleCount['sampleCount'] * self.dist_mat.shape[0]) - sampleCount['sampleCount']    
        
        # Make dict counter for pops --- TEST!!!
        if popType=='super':
            sup_sums = {pop : 0 for pop in self.pop_info[0]}   # Placeholder for sum values
            sup_sample_info = dict([val[0], val[1]] for val in self.sample_info.values())
            
        elif popType=='sub':
            sub_sums = {pop : 0 for pop in self.pop_info[1]}   # Placeholder for sum values
            sub_sample_info = dict([val[0], val[2]] for val in self.sample_info.values())
            
        elif popType=='all':
            sup_sums = {pop : 0 for pop in self.pop_info[0]}   # Placeholder for sum values
            sup_sample_info = dict([val[0], val[1]] for val in self.sample_info.values())
            
            sub_sums = {pop : 0 for pop in self.pop_info[1]}   # Placeholder for sum values
            sub_sample_info = dict([val[0], val[2]] for val in self.sample_info.values())
            
        #     sums1 = {pop : 0 for pop in pop_info[0]}   # Placeholder for sum values
        #     sums2 = {pop : 0 for pop in pop_info[1]}
        #     sums = dict(sums1, **sums2)
            
        #     sample_info1 = dict([val[0], val[1]] for val in sample_info.values())
        #     sample_info2 = dict([val[0], val[2]] for val in sample_info.values())
            
        #     sample_info = {**sample_info1, **sample_info2}
            
        #     return sample_info1
        else: 
            raise ValueError("Invalid calssification type. ")
            
        # map nonZero and sample_info together into matrix
        
        nonZeros = self.calcNonZerosForSamples(self.dist_mat, percent = True)  
        if popType == 'all':
            nonZeros['supPop'] = nonZeros.index.map(sup_sample_info)
            nonZeros['subPop'] = nonZeros.index.map(sub_sample_info)
            for key, val in nonZeros.iterrows():
                sup_sums[val[2]] += val[0]
                sub_sums[val[3]] += val[0]
        
            sums = {**sup_sums,**sub_sums}
            
        elif popType == 'super':
            raise ValueError("This option is not yet available")
        else: 
            raise ValueError("This option is not yet available")
            
        # Divide non-zero values for pop on total pop compariosns    
        sampleCount['popSums'] = sampleCount.index.map(sums)
        sampleCount['percentNonZeroForPop'] = sampleCount.popSums / sampleCount.totalPopComp
        sampleCount['gene'] = self.getGeneName()
        
        return sampleCount
    
    
    def getPopDists(self):
        if not self.pop_dists: self.calcPopDists()
        return self.pop_dists
    
    def getMeanTypeDists(self):
        if not self.mean_type_dists: self.calcMeanTypeDists()
        return self.mean_type_dists
        
    def getMeanPopDists(self): 
        if not self.mean_pop_dists: self.calcMeanPopDists()
        return self.mean_pop_dists
        
    def getSingleSuperSDR(self): return self.singleSuperSDR
    def getSingleSubDR(self): return self.singleSubSDR
    def getSDRsuper(self): return self.SDRsuper
    def getSDRsub(self): return self.SDRsuper
    def getSDVsuper(self): return self.SDVsuper
    def getSDVsub(self): return self.SDVsub


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



