# -*- coding: utf-8 -*-
"""
Created on Fri May 21 09:36:30 2021

@author: norab
"""

  # Re-run: 
    runSDRorSDVpipe(cd_dir_redhood_continue, pop_info, SDR_unprocessed_redhood_allGroups, SDR_outputfile = SDR_output_redhood_allGroups, calc_single_groups=True)
    # ----------------------------------------------------------------------------
    
    # Run
    
    # Test
    runSDRorSDVpipe(cd_dir_local, pop_info, SDR_unprocessed_redhood, SDR_outputfile = SDR_output_redhood_allGroups, calc_single_groups=True)
    
    # Real SDR groups
    runSDRorSDVpipe(cd_dir_redhood, pop_info, SDR_unprocessed_redhood_allGroups, SDR_outputfile = SDR_output_redhood_allGroups, calc_single_groups=True)
    # Real SDV
    runSDRorSDVpipe(cd_dir_redhood, pop_info, SDV_unprocessed_redhood, SDV_outputfile = SDV_output_redhood)