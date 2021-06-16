# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 05:48:23 2021

@author: norab
"""


from geneTree.treeMetrics import treeMetrics
import os
from os.path import isfile, join
import csv
from datetime import datetime
import pandas as pd

class RunStuff():
    
    def __init__(self):
        pass

    def make_filelist(self, input_files):

        
        if isinstance(input_files, str):     
            files = [join(input_files, f) for f in os.listdir(input_files) 
                         if isfile(join(input_files, f))]
        return files

    def run_calcSDR(self,
                    input_files, 
                    pop_info,
                    SDRsuper_output,
                    SDRsub_output,
                    unprocessed_files):
        

        if '.csv' not in input_files:
            file_list = self.make_filelist(input_files)

        else: 
            file_list = pd.read_csv(input_files, header = None)
            file_list = list(file_list[0])
            
        ind = 0
        ind_len = len(file_list)
            
        print("Files to procescs:\n", file_list)

        try:
            while file_list:
                dist_mat = file_list.pop().strip()
                
                print("File processed: ", dist_mat)   # TO DO: convert to log  
                print(f"Number: {ind} / {ind_len}")
                
                tree = treeMetrics()
                tree.setup(dist_mat, pop_info)
                tree.calcSDR()
                
                supSDR = tree.getSDRsuper()
                subSDR = tree.getSDRsub()
                
                print("Super SDR: ", supSDR)
                print("Sub SDR: ", subSDR)
                
                ind += 1
                
                supSDR = [tree.getGeneName(), supSDR]
                subSDR = [tree.getGeneName(), subSDR]
        
                with open(SDRsuper_output, 'a', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(supSDR)
                    
                with open(SDRsub_output, 'a', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(subSDR)
                
        except: 
            # Write remaining filelist to file
            print("Disrupted.")
            print("Last file processed:", dist_mat)
            file_list.to_csv(unprocessed_files, index = False, header = False)
            with open(unprocessed_files, 'w', newline='') as f:
                writer = csv.writer(f)
                for f in file_list: 
                    writer.writerow(f)


    def run_calcSingleSDRs(self, 
                           input_files, 
                           pop_info,
                           SSDR_output_dir,
                           unprocessed_files):
        
    """
    Save one file for each gene. 
    File example: 
        
        Level   population    SingleSDR
        Super   AFR           0.4
        Super   EUR           0.5
        Sub     FIN           0.8
        ...     ...           ...
        
    """
        
        if '.csv' not in input_files:
            file_list = self.make_filelist(input_files)

        else: 
            file_list = pd.read_csv(input_files, header = None)
            file_list = list(file_list[0])
            
        ind = 1
        ind_len = len(file_list)
            
        print("Files to procescs:\n", file_list)

        try:
            while file_list:
                dist_file = file_list.pop().strip()
                
                print("File processed: ", dist_file)   # TO DO: convert to log  
                print(f"Number: {ind} / {ind_len}")
                
                tree = treeMetrics()
                tree.setup(dist_file, pop_info)
                tree.calcSingleSDRs()
                
                supSDRs = tree.getSingleSuperSDRs()
                subSDRs = tree.getSingleSubSDRs()
                
                print("Super SDR: ", supSDRs)
                print("Sub SDR: ", subSDRs)
                
                # TEST THIS WORKS
                df_sup = pd.DataFrame([['super', pop, val] for pop, val in supSDRs.items()],
                                   columns=['level', 'pop', 'singleSDR'])
                df_sub = pd.DataFrame([['sub', pop, val] for pop, val in subSDRs.items()],
                                   columns=['level', 'pop', 'singleSDR'])
                
                # Write to file
                df = pd.concat([df_sup, df_sub], ignore_index = True, axis = 0)            
                SSDR_output_file = (SSDR_output_dir + 'SSDR_' + tree.gene_name()) +'.csv'                
                df.to_csv(SSDR_output_file, index = False, header = True)  # Test
                
                ind +=1
                
        except: 
            # Write remaining filelist to file
            print("Disrupted.")
            print("Last file processed:", dist_mat)
            file_list.to_csv(unprocessed_files, index = False, header = False)
            with open(unprocessed_files, 'w', newline='') as f:
                writer = csv.writer(f)
                for f in file_list: 
                    writer.writerow(f)

       
    def run_calcSDV(self, 
                    input_files,
                    pop_info,    
                    SDVsuper_output,
                    SDVsub_output,
                    unprocessed_files):
        """
        Input: 
            input_files: path to files containing cophenetic distance matrices
            pop_info: path to population classes information file
            output_path: path to file SDR values are to be written to. 
            unprocessed: path to file for saving list of unproccessed file upon disruption
            
        Function: 
            Pipeline of function for SDR calculations. 
            Writes SDR values to file with specified gene name.
            
            If function is disrupted, a list of the remaining files to 
            calculate SDR for from the filelist is written to a file. 
        
        Output: 
            Writes SDR to file.         
        """
        print(input_files)
        if '.csv' not in input_files:
            file_list = self.make_filelist(input_files)

        else: 
            file_list = pd.read_csv(input_files, header = None)
            file_list = list(file_list[0])
            
        ind = 1
        ind_len = len(file_list)
            
        print("Files to procescs:\n", file_list)

        try:
            while file_list:
                dist_mat = file_list.pop().strip()
                
                print("File processed: ", dist_mat)   # TO DO: convert to log  
                print(f"Number: {ind} / {ind_len}")
                
                tree = treeMetrics()
                tree.setup(dist_mat, pop_info)
                
                tree.calcSDV()
                
                supSDV = tree.getSDVsuper()
                subSDV = tree.getSDVsub()
                
                print("Super SDV: ", supSDV)
                print("Sub SDV: ", subSDV)
                
                ind += 1
                
                supSDV = [tree.getGeneName(), supSDV]
                subSDV = [tree.getGeneName(), subSDV]
        
                with open(SDVsuper_output, 'a', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(supSDV)
                    
                with open(SDVsub_output, 'a', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(subSDV)
                
        except: 
            # Write remaining filelist to file
            print("Disrupted.")
            print("Last file processed:", dist_mat)
            with open(unprocessed_files, 'w', newline='') as f:
                writer = csv.writer(f)
                for f in file_list: 
                    writer.writerow(f)

    def run_CalcNonZeroForPop(self, 
                              input_files,
                              output_file,
                              pop_info):
        """
        Input: input file for cd files
        Output: file path to store values. 
            Three columns: gene name, pop and nonZero pop value. 
        Return: None
    
        """
        pass
        # # List of files
        # if isinstance(cd_file, str):     
        #     file_list = [join(cd_file, f) for f in listdir(cd_file) if isfile(join(cd_file, f))]
        # else:            # Continue from disruption
        #     file_list = self.cd_file.copy()
    
        # ind = 0
        # print("Files to process:\n", file_list)
        
        # try:
        #     while file_list:
        #         dist_file = file_list.pop().strip()
        #         print("File processed: ", dist__file)
        #         print("Number: ", ind)
                
        #         # make tree
        #         tree = treeInfo(dist_file,
        #                         pop_info)
        
        #         result = tree.calcNonZerosForPops(popType='all', percent = True)
                
        #         with open(output_file, 'a') as f:
        #             result.to_csv(output_file, mode = 'a', index_label = 'pop', header=f.tell()==0)
        #         ind += 1
        # except: 
        #     print("Disrupted. Something wrong.. ")
            

if __name__ == '__main__':
    
# =============================================================================
#     SDR
# =============================================================================
    
    # redhood_input_files =  'E:\\Master\\cophenetic_dists'
    # redhood_input_files_continue1 = 'E:\\Master\\current_run\\unprocessed_genes_09.06_14.40.csv'
    # redhood_input_files_continue2 = 'C:\\Users\\norab\\MasterDisaster\\Data\\runstop_save\\unprocessed_files_SDRcalc_09.06.2021_16.57_update.csv'
    # redhood_input_files_continue3 = 'E:\\Master\\current_run\\unprocessed_files_10.06.2021_13.37.csv'
    # redhood_input_files_continue4 = 'E:\\Master\\current_run\\unprocessed_genes_11.06.21_11.40.csv'
    # redhood_input_files_continue5 = 'E:\\Master\\current_run\\unprocessed_files_12.06.21_11.51.csv'
    
    # pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
    # SDRsuper = 'E:\\Master\\current_run\\SDRsuper_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    # SDRsub = 'E:\\Master\\current_run\\SDRsub_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    # save_unprocessed = 'C:\\Users\\norab\\MasterDisaster\\Data\\runstop_save\\unprocessed_files_SDRcalc_{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))
    
    # run = RunStuff()

    # run.run_calcSDR(input_files = redhood_input_files_continue5, 
    #                 pop_info = pop_info, 
    #                       SDRsuper_output = SDRsuper, 
    #                       SDRsub_output = SDRsub, 
    #                       unprocessed_files = save_unprocessed
    #                       )
    
# =============================================================================
#     SDV
# =============================================================================
    # Testrun
    # redhood_input_files_test =  'E:\\Master\\test_runs\\cophenetic_dists_selection'
    # pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
    # SDVsuper = 'E:\\Master\\test_runs\\SDVsuper_runDate_{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    # SDVsub = 'E:\\Master\\test_runs\\SDVsub_runDate_{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    # save_unprocessed = 'C:\\Users\\norab\\MasterDisaster\\Data\\runstop_save\\unprocessed_files_SDVcalc_{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))

        
    redhood_input_files =  'E:\\Master\\cophenetic_dists'
    pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
    SDVsuper = 'E:\\Master\\current_run\\SDVsuper_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    SDVsub = 'E:\\Master\\current_run\\SDVsub_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    save_unprocessed = 'C:\\Users\\norab\\MasterDisaster\\Data\\runstop_save\\unprocessed_files_SDRcalc_{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))
    
    run = RunStuff()
    
    run.run_calcSDV(input_files = redhood_input_files, 
                          pop_info = pop_info, 
                          SDVsuper_output = SDVsuper, 
                          SDVsub_output = SDVsub, 
                          unprocessed_files = save_unprocessed
                          )
    
