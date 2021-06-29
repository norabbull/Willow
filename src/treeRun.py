# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 05:48:23 2021

@author: norab
"""

import csv
from datetime import datetime
import pandas as pd
import yaml
import sys
import os
from os.path import isfile, join
from treeMetrics import treeMetrics
import logging
import logging.config
from configparser import ConfigParser

class RunStuff:
    
    def __init__(self, configFilepath):
        
        # Load config arguments
        configFilepath = configFilepath.strip()
        with open(configFilepath, 'r') as c:
            self.config_file = yaml.safe_load(c)
            
        for key, val in self.config_file.items():
            if 'datetime' in val: 
                self.config_file[key] = val.replace('datetime', datetime.now().strftime("%d.%m.%Y_%H.%M"))
        
        # Set function to run
        self.func = self.config_file.get('func')
        
        # Initiate logger 
        self.logger_details = self.read_config(self.config_file['log_config'], self.config_file['LOGGER_SECTION'])
        
    def read_config(self,file_path,section):
        """
        The function to read the configuration parameters from the relevant section in the config file.
        Parameters returned in dictionary. 
        
        :param file_path: the name of the file that will be read.
        :param section: The section of the file has to be read.
        :return section_content: Dictionary with the section parameters.
        :rtype: Dictionary
        
        code taken from: https://medium.com/analytics-vidhya/config-files-and-logging-a7b4cc377fd5
        """
        
        # create parser and read configuration file         
        parser = ConfigParser()
        parser.read(file_path)
                
        section_content = {}
        if parser.has_section(section):
            items = parser.items(section)
            for item in items:
                section_content[item[0]] = item[1]
        else:
            raise Exception('Error in fetching config in read_config method. {0} not found in \
            {1}'.format(section, file_path))

        return section_content
    
    def set_logging_basics(self):
        """
        The function to set the logging information.
        :param config_dict: dictionary with details of configuration
        :return: None
        """
        folder = self.config_dict.get('folder')
        file = self.config_dict.get('file')
        level = self.config_dict.get('level')
        func = self.config_file.get('func')
        self.logging.basicConfig(filename=folder + file
                            , filemode='w+'
                            , format='%(name)s - %(levelname)s - %(message)s'
                            , level=level)

    def make_filelist(self, input_files):

        if isinstance(input_files, str):     
            files = [join(input_files, f) for f in os.listdir(input_files) 
                         if isfile(join(input_files, f))]
        return files
        
    def run_calcSDR(self, random=False):
        
        input_files = self.config_file.get('input_files').strip()
        pop_info = self.config_file.get('pop_info').strip()
        SDRsuper_output = self.config_file.get('SDRsuper_output').strip()
        SDRsub_output = self.config_file.get('SDRsub_output').strip()
        save_unprocessed = self.config_file.get('save_unprocessed').strip()
            
        file_list = self.make_filelist(input_files)
            
        ind = 1
        ind_len = len(file_list)
            
        print("Files to procescs:\n", file_list)
  
        for dist_file in file_list:
            try:         
                dist_file = dist_file.strip()
                
                print("File processed: ", dist_file)   # TO DO: convert to log  
                print(f"Number: {ind} / {ind_len}")
                ind += 1
                
                tree = treeMetrics()
                tree.setup(dist_file, pop_info)
                
                if random: 
                    tree.shuffleSampleInfo()
                
                tree.calcSDR()
                
                supSDR = tree.getSDRsuper()
                subSDR = tree.getSDRsub()
                
                print("Super SDR: ", supSDR)
                print("Sub SDR: ", subSDR)
                print("1st sample: ", tree.getSampleInfo()[0])
                
                supSDR = [tree.getGeneName(), supSDR]
                subSDR = [tree.getGeneName(), subSDR]
        
                with open(SDRsuper_output, 'a', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(supSDR)
                    
                with open(SDRsub_output, 'a', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(subSDR)

            except Exception: 
               
                self.logger.exception(f"File disrupted: {dist_file} ")    
                with open(save_unprocessed, 'w+') as f: 
                    write = csv.writer(f) 
                    write.writerow(dist_file)
                    
                pass
        

       
    def run_calcSingleSDRs(self, random = False):
        
        """
        Save one file for each gene. 
        File example: 
            
            Level   population    SingleSDR
            Super   AFR           0.4
            Super   EUR           0.5
            Sub     FIN           0.8
            ...     ...           ...
            
        """

        input_files = self.config_file.get('input_files')
        pop_info = self.config_file.get('input_pop_info')
        SSDR_output_dir = self.config_file.get('output_SSDR')
        save_unprocessed = self.config_file.get('output_unprocessed')
               
        file_list = self.make_filelist(input_files)

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
                if random: 
                    tree.shuffleSampleInfo()
                    
                tree.calcSingleSDRs()
                
                supSDRs = tree.getSingleSuperSDRs()
                subSDRs = tree.getSingleSubSDRs()

                df_sup = pd.DataFrame([['super', pop, val] for pop, val in supSDRs.items()],
                                   columns=['level', 'pop', 'singleSDR'])
                df_sub = pd.DataFrame([['sub', pop, val] for pop, val in subSDRs.items()],
                                   columns=['level', 'pop', 'singleSDR'])
                
                # Write to file
                df = pd.concat([df_sup, df_sub], ignore_index = True, axis = 0)
                SSDR_output_file = (SSDR_output_dir + 'SSDR_' + tree.getGeneName()) +'.csv'                
                df.to_csv(SSDR_output_file, index = False, header = True,na_rep = "NULL")  # Test
                
                ind +=1
                

        except Exception: 
           
            self.logger.exception(f"File disrupted: {dist_file} ")    
            with open(save_unprocessed, 'w+') as f: 
                write = csv.writer(f) 
                write.writerow(dist_file)
                
            pass

       
    def run_calcSDV(self, random = False):
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
        input_files = self.config_file.get('input_files').strip()
        pop_info = self.config_file.get('input_pop_info').strip()
        SDVsuper_output = self.config_file.get('output_SDRsuper').strip()
        SDVsub_output = self.config_file.get('output_SDRsub').strip()
        save_unprocessed = self.config_file.get('output_unprocessed').strip()
            
        file_list = self.make_filelist(input_files)

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
                
                if random: 
                    tree.shuffleSampleInfo()
                    
                tree.calcSDV()
                
                supSDV = tree.getSDVsuper()
                subSDV = tree.getSDVsub()
                
                print("Super SDV: ", supSDV)
                print("Sub SDV: ", subSDV)
                
                ind += 1
                
                supSDV = [tree.getGeneName(), supSDV]
                subSDV = [tree.getGeneName(), subSDV]
        
                with open(SDVsuper_output, 'w+', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(supSDV)
                    
                with open(SDVsub_output, 'w+', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(subSDV)
                
        except: 
            # Write remaining filelist to file
            print("Disrupted.")
            print("Last file processed:", dist_file)
            
            with open(save_unprocessed, 'w+') as f: 
                write = csv.writer(f) 
                write.writerow(file_list) 

       
    def run_calcTest(self):
        """
        Test function to check functionality of program.
        """

        input_files = self.config_file.get('input_files').strip()
        pop_info = self.config_file.get('input_pop_info').strip()
        measure_output = self.config_file.get('output_measure').strip()
        save_unprocessed = self.config_file.get('output_unprocessed').strip()
            
        file_list = self.make_filelist(input_files)
        
        print("File list: ")
        print(file_list)

        dist_file = file_list.pop().strip()
        print("File processed: ", dist_file)   # TO DO: convert to log  

        tree = treeMetrics()
        tree.setup(dist_file, pop_info)
        print("Gene name: ", tree.getGeneName())
        print("Sample info: ", tree.getSampleInfo())
        
        rand_num = 15.05
        
        with open(measure_output, 'w+', newline='') as f:   # write to file    
            writer = csv.writer(f)
            writer.writerow(rand_num)
        
        

        
    # def run_CalcNonZeroForPop(self):
    #     """
    #     Input: input file for cd files
    #     Output: file path to store values. 
    #         Three columns: gene name, pop and nonZero pop value. 
    #     Return: None
    
    #     """
    #     pass
    #     # # List of files
    #     # file_list = self.make_filelist(input_files)
    
    #     # ind = 0
    #     # print("Files to process:\n", file_list)
        
    #     # try:
    #     #     while file_list:
    #     #         dist_file = file_list.pop().strip()
    #     #         print("File processed: ", dist__file)
    #     #         print("Number: ", ind)
                
    #     #         # make tree
    #     #         tree = treeInfo(dist_file,
    #     #                         pop_info)
        
    #     #         result = tree.calcNonZerosForPops(popType='all', percent = True)
                
    #     #         with open(output_file, 'a') as f:
    #     #             result.to_csv(output_file, mode = 'a', index_label = 'pop', header=f.tell()==0)
    #     #         ind += 1
    #     #         except: 
    #         # # Write remaining filelist to file
    #         # print("Disrupted.")
    #         # print("Last file processed:", dist_file)
            
    #         # with open(save_unprocessed, 'w') as f: 
    #         #     write = csv.writer(f) 
    #         #     write.writerow(file_list) 
    

    def main(self):
        
        logging.info('New session...', extra={'function':self.func})
        
        if self.func == "calcSDR":
            self.run_calcSDR()
        elif self.func == "calcSDV":
            self.run_calcSDV()
        elif self.func == "calcSingleSDRs":
            self.run_calcSingleSDRs()
        elif self.func == "calcNullSDR":
            num_rand_trees = int(self.config_file.get('num_rand_trees'))
            for _ in range(num_rand_trees):
                self.run_calcSDR(random = True)
        elif self.func == "calcTest":
            self.run_calcTest()
            
        else:
            self.logger.error("Not a valid option.")
        
        logging.info('Session finished: ', datetime.now().strftime("%d.%m.%Y_%H.%M"))

if __name__ == '__main__':
    
    # Test on lenovo computer
    #configFilepath = 'E:/Master/jobs/debug_calcSDRnullDist_28.06.yml'
    
    configFilepath = sys.argv[1]
    run = RunStuff(configFilepath)

    run.main()
    
    
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
    #                       save_unprocessed = save_unprocessed
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

        
    # redhood_input_files =  'E:\\Master\\cophenetic_dists'
    # pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
    # SDVsuper = 'E:\\Master\\current_run\\SDVsuper_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    # SDVsub = 'E:\\Master\\current_run\\SDVsub_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
    # save_unprocessed = 'C:\\Users\\norab\\MasterDisaster\\Data\\runstop_save\\unprocessed_files_SDVcalc_{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))
    
    # run = RunStuff()
    
    # run.run_calcSDV(input_files = redhood_input_files, 
    #                       pop_info = pop_info, 
    #                       SDVsuper_output = SDVsuper, 
    #                       SDVsub_output = SDVsub, 
    #                       save_unprocessed = save_unprocessed
    #                       )
    
    
# =============================================================================
#     SingleSDRs
# =============================================================================
    
    # Testrun
    # input_files = 'E:\\Master\\test_runs\\cophenetic_dists_selection'
    # pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
    # SSDR_output_dir = 'E:\\Master\\test_runs\\singleSDR_test\\'
    # save_unprocessed = 'E:\\Master\\test_runs\\singleSDR_test\\unprocessed_testRun17.06.21_16.00.csv'
    
    # run = RunStuff()
    # run.run_calcSingleSDRs(
    #     input_files = input_files, 
    #     pop_info = pop_info,
    #     SSDR_output_dir = SSDR_output_dir
    #     save_unprocessed = save_unprocessed)
    
    # Real run
        
    # redhood_input_files =  'E:\\Master\\cophenetic_dists'
    # continue_files1 = 'C:\\Users\\norab\\MasterDisaster\\Data\\runstop_save\\unprocessed_files_SSDRcalc_17.06.2021_17.27.csv'

    # pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
    # SSDR_output_dir = 'E:\\Master\\current_run\\singleSDRs\\'
    # save_unprocessed = 'C:\\Users\\norab\\MasterDisaster\\Data\\runstop_save\\unprocessed_files_SSDRcalc_{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))
    
    # run = RunStuff()
    
    # run.run_calcSingleSDRs(input_files = continue_files1, 
    #                       pop_info = pop_info, 
    #                       SSDR_output_dir = SSDR_output_dir,
    #                       save_unprocessed = save_unprocessed
    #                       )
    
    