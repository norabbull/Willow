# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 13:14:01 2021

@author: norab
"""


import sys
import os

# Make sure system modules are contained in sys.path
abspath = os.path.abspath('.')
add_paths = ['C:\\Users\\norab\\Master\\Willow']
sys.path = list(set(sys.path + add_paths))
