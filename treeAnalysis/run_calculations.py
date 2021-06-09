# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 05:48:23 2021

@author: norab
"""


from treeAnalysis.dispersionMetrics import treeDispersion
import os
from os.path import isfile, join
import csv
from datetime import datetime
import pandas as pd

class RunStuff(treeDispersion):
    
    def __init__(self):
        treeDispersion.__init__(self)

    def make_filelist(self,input_files):
    
        if isinstance(input_files, str):     
            files = [join(input_files, f) for f in os.listdir(input_files) 
                         if isfile(join(input_files, f))]
        return files


    def run_calcSDR(self,input_files, 
                    pop_info,
                    SDRsuper_output,
                    SDRsub_output,
                    unprocessed_files):
        

        if '.csv' not in input_files:
            file_list = self.make_filelist(input_files)
            ind = 0
        else: 
            file_list = pd.read_csv(input_files, header = None)
            file_list = list(file_list[0])
            ind = len(file_list)
            
        print("Files to process:\n", file_list)

        try:
            while file_list:
                dist_mat = file_list.pop().strip()
                
                print("File processed: ", dist_mat)   # TO DO: convert to log
                print("Number: ", ind)
                
                tree = treeDispersion()
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
            with open(unprocessed_files, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(file_list)


    def run_calcPopSDRs():
        pass
        # if calc_SDRpop: 
        #     SDRs = SDRs[0].append(SDRs[1]).to_frame()
        #     gene_name = tree.getGeneName()
        #     SDRs = SDRs.assign(gene = [gene_name] * SDRs.shape[0])
                
        #     with open(SDR_outputfile, 'a') as f:
        #         SDRs.to_csv(SDR_outputfile, mode = 'a', index_label = 'pop', header=f.tell()==0)      
                
    def run_calcSDV(self, 
                    input_files,
                    pop_info,                  
                    unprocessed, 
                    SDV_outputfile=None):
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
        pass
        # List of files (make absolute path name of each file in dir)
        # if isinstance(input_files, str):     
        #     file_list = [join(input_files, f) for f in listdir(input_files) 
        #                  if isfile(join(input_files, f))]
        # else:            
        #     file_list = input_files.copy() # Continue from disruption
    
        # ind = 0
        # print("Files to process:\n", file_list)
        
        # try:
        #     while file_list:
        #         dist_mat = file_list.pop().strip()
        #         print("File processed: ", cd_file)
        #         print("Number: ", ind)
                
        #         tree = treeInfo(dist_mat, pop_info)

                
        #         # SDV
        #         if SDV_outputfile:
        #             SDRpop = calculateSDR(cd_mat,
        #                                      sample_info, 
        #                                      pop_type_sets, 
        #                                      calc_SDRpop=calc_SDRpop)
                    
        #             SDVs = calculateSDV(SDRpop)
        #             SDVs.insert(0,getGeneName(cd_file))
    
        #             with open(SDV_outputfile, 'a', newline='') as f:   # write to file
        #                 writer = csv.writer(f)
        #                 writer.writerow(SDVs)
        #         ind += 1       
         
        # except: 
        #     # Write remaining filelist to file
        #     print("Disrupted.")
        #     print("Last file processed:", cd_file)
        #     with open(unprocessed_files, 'w', newline='') as f:
        #            writer = csv.writer(f)
                   # writer.writerow(file_list)


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
    
    
    redhood_input_files =  'E:\\Master\\cophenetic_dists'
    redhood_input_files_continue1 = 'E:\\Master\\current_run\\unprocessed_genes_09.06_14.40.csv'
    pop_info = 'C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'
    SDRsuper = 'E:\\Master\\current_run\\SDRsuper_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    SDRsub = 'E:\\Master\\current_run\\SDRsub_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    save_unprocessed = 'C:\\Users\\norab\\MasterDisaster\\Data\\runstop_save\\unprocessed_files_SDRcalc_{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))
    
    run = RunStuff()
    
    # run.run_calcSDR(input_files = redhood_input_files, 
    #                      pop_info = pop_info, 
    #                      SDRsuper_output = SDRsuper, 
    #                      SDRsub_output = SDRsub, 
    #                      unprocessed_files = save_unprocessed)
    
    # Run 2, continue form disruption 09.06
    run.run_calcSDR(input_files = redhood_input_files_continue1, 
                         pop_info = pop_info, 
                         SDRsuper_output = SDRsuper, 
                         SDRsub_output = SDRsub, 
                         unprocessed_files = save_unprocessed)
