# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 14:48:34 2021

@author: norab
"""

import re
import pandas as pd
import matplotlib.pyplot as plt


#%% import data
uniqseq = pd.read_csv("9381_uniqseqs.txt", "r", delimiter = ":", header = 0) 
totdist = pd.read_csv("totalDistancesRefined.txt", "r", delimiter = ",", index_col = 0)


# Rename inds
totdist.rename(index=lambda s: re.sub('_HUMAN__full.*', '', s), columns = {totdist.columns[0]: "totdist"}, inplace = True)
totdist.rename(index=lambda s: re.sub('C.*trees/', '', s), columns = {totdist.columns[0]: "totdist"}, inplace = True)
uniqseq.rename(index=lambda s: re.sub('_HUMAN__uniq.*', '', s), columns = {uniqseq.columns[0]: "uniqseq"}, inplace = True)

#%% Create dict mapping index to gene name

# Join dfs
result = pd.concat([totdist, uniqseq], axis=1, join="inner")

# Drop na rows
result.dropna(subset=[result.columns[0]], inplace = True)
result.dropna(subset=[result.columns[1]], inplace = True)

#%%

geneNames = totdist.index.values.tolist()
ind = list(range(0,len(geneNames)))
geneNamesDict = dict(zip(ind, geneNames))

uniqseq.rename(index=geneNamesDict)

#%% Plot

result = result.sort_values('uniqseq')
#plt.plot(result2.uniqseq)

plt.figure()
result.plot(x='uniqseq', y='totdist', style='o')
plt.xlabel("unique sequences")
plt.ylabel("total non-zero distances in tree")
plt.show()

result[['uniqseq', 'totdist']].corr()
f = plt.figure(figsize=(19, 15))
#plt.matshow(result['totdist'].corr(result['uniqseq']))
#plt.title('Correlation Matrix', fontsize=16)



#%% Checking SDR for one tree

# Check which sample I calculated SDR values for
sample_RALB = result.loc["ENSG00000144118___RALB"]

# Tree: ENSG00000144118___RALB_HUMAN__full_aln.tre

oneSDR = pd.read_csv("C:/Users/norab/MasterDisaster/Data/SDR/test/one_sample_SDR.txt", 
                     "r", delimiter = "\\n",header=None)

df = pd.DataFrame(columns = ['Group', 'SDR'])

with open("C:/Users/norab/MasterDisaster/Data/SDR/test/one_sample_SDR.txt", "r") as file: 
    for line in file:
        line = line.strip()

        if "Group" in line:
            group = line.replace("Group: ", "")

        elif "SDR" in line: 
            sdr = line.replace("SDR:", "")
            sdr = float(sdr.replace(",", "."))
            df = df.append({'Group': group, 'SDR': sdr},
                      ignore_index = True)
        else: 
            pass
            
print(df)

