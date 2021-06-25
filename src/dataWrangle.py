# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 06:45:58 2021

@author: norab
"""

from treeAnalysis.treeInformation import treeInfo
import pandas as pd
import os
from os.path import isfile, join
import re


#%% 
test_genes = ['E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000134183___GNAT2___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000105669___COPE___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000164828___SUN1___CopD.csv',
             'E:/Master/cophenetic_dists/ENSG00000164654___MIO___CopD.csv']

pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
    

#%% 

# =============================================================================
#  Make dataframe subset, 10 < 19 dist mat
# =============================================================================
    
test_files1 = {'pop_info': pop_info,
              'dist_mat': test_genes[0]}

tree1 = treeInfo(test_files1['dist_mat'], test_files1['pop_info'])
mat = tree1.getDistMat()

submat = mat.iloc[0:10,0:10]
submat.to_csv('C:/Users/norab/Master/Data/real_tree_data/dist_mat_test/FGR_10x10.csv')

# Check if equal: 
submat2 = pd.read_csv('C:/Users/norab/Master/Data/real_tree_data/dist_mat_test/FGR_10x10.csv', index_col = 0)
(submat == submat2).eq(True).all().all()


# =============================================================================
# Make dataframe subset, 5 x 5 
# =============================================================================

test_files2 = {'pop_info': pop_info,
              'dist_mat': test_genes[4]}

tree2 = treeInfo()
tree2.setup(test_files2['dist_mat'], test_files2['pop_info'])
mat2 = tree2.getDistMat()

submat3 = mat2.iloc[132:137,132:137]
submat3['EUR___GBR___HG00236']['EUR___GBR___HG00159'] = 0.005
submat3['EUR___GBR___HG00159']['EUR___GBR___HG00236'] = 0.005
submat3['AMR___MXL___NA19758']['AMR___MXL___NA19795'] = 0.003
submat3['AMR___MXL___NA19795']['AMR___MXL___NA19758'] = 0.003


submat3.to_csv('C:/Users/norab/Master/Data/real_tree_data/dist_mat_test/MIO_5x5.csv')

# Check if equal: 
submat4 = pd.read_csv('C:/Users/norab/Master/Data/real_tree_data/dist_mat_test/MIO_5x5.csv', index_col = 0)
(submat3 == submat4).eq(True).all().all()



# =============================================================================
#  Make list of unprocces files from SDR file
# =============================================================================

# Interruption 1, 09.06.14.40
processed_genes = pd.read_csv('E:\Master\current_run\SDRsuper_runDate08.06.2021_15.57.csv', index_col = 0)
processed_genes = pd.Series(processed_genes.index) 
# Add ful path

processed_genes.to_csv('E:\\Master\\current_run\\processed_genes_09.06_14.40.csv', index = False)

# all genes
all_genes = make_filelist('E:\Master\cophenetic_dists')
unprocessed_genes = all_genes.copy()

for file in all_genes: 
    name = geneName(file)
    for g in processed_genes: 
        if name in g: 
            unprocessed_genes.remove(file)

unprocessed_genes = pd.Series(unprocessed_genes)
unprocessed_genes.to_csv('E:\\Master\\current_run\\unprocessed_genes_09.06_14.40.csv', index = False)

# 10.06, 13.16 ----------------------------------------------------------------
processed_genes = list(pd.read_csv('E:\Master\current_run\SDRsub_runDate08.06.2021_15.57.csv', index_col = 0, header = None).index)
processed_genes2 = list(pd.read_csv('E:\Master\current_run\SDRsub_runDate09.06.2021_16.57.csv', index_col = 0, header = None).index)
processed_genes3 = list(pd.read_csv('E:\Master\current_run\SDRsub_runDate10.06.2021_03.54.csv', index_col = 0, header = None).index)
processed_genes4 = list(pd.read_csv('E:\Master\current_run\SDRsub_runDate10.06.2021_13.46.csv', index_col = 0, header = None).index)



all_processed_genes = list(set(processed_genes + processed_genes2 + processed_genes3 + processed_genes4))
#all_processed_genes = processed_genes + processed_genes2 + processed_genes3 + processed_genes4
# Add ful path

pd.Series(all_processed_genes).to_csv('E:\\Master\\current_run\\processed_genes_10.06.2021_13.35.csv', index = False)

# all genes
all_genes = list(set(make_filelist('E:\Master\cophenetic_dists')))
unprocessed_genes = all_genes.copy()

for file in all_genes: 
    #name = geneName(file)
    for gene in all_processed_genes: 
        if gene in file: 
            try:
                unprocessed_genes.remove(file)
            except:
                print(f"{file} not removed. ")
            
        

ug = pd.Series(unprocessed_genes)
ug.to_csv('E:\\Master\\current_run\\unprocessed_genes_11.06.21_11.40.csv', index = False)

# Testing stuff
for file in unprocessed_genes: 
    if "ENSG00000110680___CALCA" in file: 
        print("File:", file)

# =============================================================================
# Make list of unprocessed files from unprocessed files list file
# =============================================================================

redhood_input_files_continue2 = 'C:\\Users\\norab\\Master\\Data\\runstop_save\\unprocessed_files_SDRcalc_09.06.2021_16.57.csv' 
liste = pd.read_csv(redhood_input_files_continue2 , header = None, sep = ",",).transpose()
liste.to_csv('C:\\Users\\norab\\Master\\Data\\runstop_save\\unprocessed_files_SDRcalc_09.06.2021_16.57_update.csv' , 
             index = False,
             header = False)
liste.to_csv(file_list, index = False,
             header = False)



# =============================================================================
# Merge SDR calculation files to one
# =============================================================================

# Extract SDRsub and SDRsuper file paths

all_files = make_filelist('E:\\Master\\current_run')
SDRsuper_files = [f for f in all_files if 'SDRsuper' in f]
SDRsub_files = [f for f in all_files if 'SDRsub' in f]

# Load
SDRsuper = pd.DataFrame(columns = ['gene', 'val'])
SDRsub = pd.DataFrame(columns = ['gene', 'val'])
for f in SDRsuper_files: 
    SDRsuper = SDRsuper.append(pd.read_csv(f, names = ['gene', 'val']))
for f in SDRsub_files: 
    SDRsub = SDRsub.append(pd.read_csv(f, names = ['gene', 'val']))

SDRsuper.drop_duplicates(inplace=True)
SDRsub.drop_duplicates(inplace=True)

# Save
SDRsuper.to_csv('E:\\Master\\current_run\\SDRsuper_all.csv', index = False, header = False)
SDRsub.to_csv('E:\\Master\\current_run\\SDRsub_all.csv', index = False, header = False)


# =============================================================================
# Filter genes with SDR <= 0.7, 0.6 and 0.5 - NOT RUN YET CUZ SDV IS RUNNING
# =============================================================================

# Super:

SDRsuper_all = pd.read_csv('C:\\Users\\norab\\Master\\Data\\SDR\\SDRsuper_all.csv',
                           index = False, header = False, names = ['gene', 'val'])
SDRsub_all = pd.read_csv('C:\\Users\\norab\\Master\\Data\\SDR\\SDRsub_all.csv',
                         index = False, header = False, names = ['gene', 'val'])

SDRsuper_07 = SDRsuper_all[SDRsuper_all['val'] <= 0.7]
SDRsub_07 = SDRsub_all[SDRsub_all['val'] <= 0.7]

SDRsuper_06 = SDRsuper_all[SDRsuper_all['val'] <= 0.6]
SDRsub_06 = SDRsub_all[SDRsub_all['val'] <= 0.6]

SDRsuper_05 = SDRsuper_all[SDRsuper_all['val'] <= 0.5]
SDRsub_05 = SDRsub_all[SDRsub_all['val'] <= 0.5]

# Save 
SDRsuper_07.to_csv('C:\\Users\\norab\\Master\\Data\\SDR\\SDRsuper_07.csv', 
                   index = False, header = False)
SDRsub_07.to_csv('C:\\Users\\norab\\Master\\Data\\SDR\\SDRsub_07.csv', 
                   index = False, header = False)
SDRsuper_06.to_csv('C:\\Users\\norab\\Master\\Data\\SDR\\SDRsuper_06.csv', 
                   index = False, header = False)
SDRsub_06.to_csv('C:\\Users\\norab\\Master\\Data\\SDR\\SDRsub_06.csv', 
                   index = False, header = False)
SDRsuper_05.to_csv('C:\\Users\\norab\\Master\\Data\\SDR\\SDRsuper_05.csv', 
                   index = False, header = False)
SDRsub_05.to_csv('C:\\Users\\norab\\Master\\Data\\SDR\\SDRsub_05.csv', 
                   index = False, header = False)

#%%  Make data files ready for inspecton 


# SDR
SDRsuper = pd.read_csv('C:/Users/norab/Master/Data/SDR/SDRsuper_all.csv', header = None, names = ['gene', 'value'])
SDRsub = pd.read_csv('C:/Users/norab/Master/Data/SDR/SDRsub_all.csv', header = None, names = ['gene', 'value'])
SDRpsuedo = pd.read_csv('C:/Users/norab/Master/Data/SDRnullDist/nullDistSDRsub_23.06.2021_09.14.csv', header = None, names = ['gene', 'value'])
SDRsuper.dropna(inplace=True)
SDRsub.dropna(inplace=True)
SDRpsuedo.dropna(inplace=True)

# Add level column
SDRsuper.insert(0, 'level', 'super')
SDRsub.insert(0, 'level', 'sub')
SDRpsuedo.insert(0, 'level', 'psuedo')

SDRall = pd.concat([SDRsub, SDRsuper, SDRpsuedo])

SDRall.to_csv('C:/Users/norab/Master/Data/SDR/SDR_all_shortPsuedo.csv')





#%% Helpers
def make_filelist(input_files):

    if isinstance(input_files, str):     
        files = [join(input_files, f) for f in os.listdir(input_files) 
                     if isfile(join(input_files, f))]
    return files

def geneName(dist_mat_file):
    """
    Input: 
        file: string filepath to distance matrix-file. 
    Function: 
        Filters out name of gene the tree represents and assign to
        class variable "gene_name".
        Both Ensembl and gene name identifiers included on the form: 
            'ENSG00000000938___FGR'
    """
    subName = re.sub('^.:.*ENS', 'ENS', dist_mat_file)
    gene_name = re.sub('___CopD.csv$','', subName)
    
    return gene_name
