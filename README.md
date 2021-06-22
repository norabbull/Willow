TreeMetrics: - Input: For each tree, a matrix of cophenetic distances between all sample pairs is given. A sample is a 3'UTR of an individual human, belonging to a super- and a sub-population group.

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
