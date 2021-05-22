# -*- coding: utf-8 -*-
"""
Created on Wed May 12 10:05:17 2021

@author: norab
"""

from inspection.loadTreeData import *

#%% Load single gene matrix

GNAT2_cd= load_cd_mat('E:/Master/cophenetic_dists/ENSG00000134183___GNAT2___CopD.csv')
COPE_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000105669___COPE___CopD.csv')
HSP7C_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000109971___HSP7C___CopD.csv')
S29A4_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000164638___S29A4___CopD.csv')
VAV2_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000160293___VAV2___CopD.csv')
FGR_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv')

#%% Load data

# All data
uniqseq = load_uniq_seqs()
totdist = load_tot_dist()
SDRs = load_SDRs('E:\Master\SDR\SDR_values_all.csv')
SDVs = load_SDVs('E:\Master\SDV\SDV_values_all.csv')
group_SDRs = read_SDR_allgroups('E:\Master\SDR\SDR_values_all_groups.csv')

# Small subset data
group_SDRs = read_SDR_allgroups('C:/Users/norab/MasterDisaster/Data/SDR/SDR_values_subset.csv')

#%% Filtering - full SDR insepction

# Remove entries where SDR = 1
group_SDRs_filter1 = group_SDRs.drop(group_SDRs[group_SDRs.SDR ==1].index)


# Join dfs
# result = pd.concat([totdist, uniqseq, SDRs, SDVs], axis=1, join="inner")
result = pd.merge(uniqseq, group_SDRs_filter1, left_index = True, right_on="gene")
result = pd.merge(totdist,result, left_index = True, right_on="gene")
# Drop na rows (dont do for analysing different columns)
result.dropna(subset=result.columns, inplace = True)
#result.dropna(subset=[result.columns], inplace = True)

# Sort for uniqseq
result2 = result.sort_values('uniqseq')
result['log_uniqseq'] = np.log10(result['uniqseq'])

# Without original uiqseq (only log)
result2 = result.drop(['uniqseq'], axis= 1)

# Without SDR > 1 for super and sub
result3 = result2[result2.SDR_super <= 1]
result3 = result3[result3.SDR_sub <= 1]

# SDR < 0.5 and > 0
result4 = result3[result3.SDR_super < 0.6]
result4 = result4[result4.SDR_sub != 0]

# totdist > 0.39 (this corresponds to 1000 samples)
result5 = result4[result4.totdist > 0.39]
result6 = result5[result5.uniqseq > 150]

#%% Filtering - groupSDR inspection



#%% Rename columns (if necessary for visual purposes)

geneNames = totdist.index.values.tolist()
ind = list(range(0,len(geneNames)))
geneNamesDict = dict(zip(ind, geneNames))
uniqseq.rename(index=geneNamesDict)

#%% Plot

# =============================================================================
# Pairplots
# =============================================================================

# Unique seq vs totdist
plt.figure()
result.plot(x='uniqseq', y='totdist', style='o')
plt.xlabel("unique sequences")
plt.ylabel("total non-zero distances in tree")
plt.show()

# Unique seq vs SDR super
plt.figure()
result.plot(x='SDR_super', y='uniqseq', style='o', color = "red")
plt.xlabel("SDR super")
plt.ylabel("unique sequences")
plt.title("SDR vs number of unique sequences in gene. No samples removed.")
plt.show()

# All against all pairplots
sns.pairplot(result2)
sns.pairplot(result3)
sns.pairplot(result4)
sns.pairplot(result5)
sns.pairplot(result6)

# All in one plot
result.plot()
result4.plot()

# =============================================================================
# Barchart
# =============================================================================

result['uniqseq'].plot.hist(bins=100, figsize=(8,6))
result['SDR_super'].plot.hist(bins=100, figsize=(8,6))

result2.plot.hist(bins=100, alpha = 0.7)


# =============================================================================
# Scatterplots
# =============================================================================

result.plot.scatter(x='uniqueseq', y='SDR_super')


g = sns.catplot(x="gene", y="SDR", hue="pop", data=result)
g.set_xticklabels(rotation=90)

g = sns.catplot(x="gene", y="SDR", hue="pop", size = "log_uniqseq", jitter = False, data=result)
g.set_xticklabels(rotation=90)

g = sns.relplot(x="gene", y="SDR", hue="pop", size="uniqseq", sizes = (1,300),data=result)
g.set_xticklabels(rotation=90)


# Pairplot of SDR groups

sns.pairplot(result, hue='pop')

# Youre onto something here! Do this with all. look for patterns.
#%% analysis of single trees / genes

# =============================================================================
# Heatmap
# =============================================================================

# VAV2 ; totdist = 0,737422 , uniqseq = 244 (OBS! Low), SDR super = 0.3339, SDR sub = 0.3061, also higher variance than others

sns.heatmap(VAV2_cd)


#%% Plot subgroup SDR distributions

# Question to anser:
    # Are there any popultiaons with stronger clustering for multiple genes than others? 
    # Plot SDV against SDR: 
        
        
# Test reading values: 




#%% plot SDRs for populations

group_SDRs.plot.hist(bins = 20)







































