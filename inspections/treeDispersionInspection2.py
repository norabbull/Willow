# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 13:27:05 2021
@author: norab
Module to inspect data. 
Loads data form files. 
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import re
from plotnine import ggplot, aes, geom_line, geom_violin
from plotnine import *

#%% Selection of plot utlis

# geom_line()   - line. X against y.


#%% SDR - load data

SDRsuper = pd.read_csv('E:\Master\current_run\SDRsuper_all.csv', header = None, names = ['gene', 'value'])
SDRsub = pd.read_csv('E:\Master\current_run\SDRsub_all.csv', header = None, names = ['gene', 'value'])
SDRsuper.dropna(inplace=True)
SDRsub.dropna(inplace=True)

# Add level column
SDRsuper.insert(0, 'Level', 'super')
SDRsub.insert(0, 'Level', 'sub')

SDRall = pd.concat([SDRsub, SDRsuper])

#%% SDRsuper and SDRsub violin plot + boxplot + lines

# =============================================================================
# Simple violin plot: 
# =============================================================================

(
 ggplot(SDRall)
 + aes(y = 'value', x = 'Level', fill = 'Level')
 + geom_violin(scale = "width")
 )

# =============================================================================
# Next level violin plots
# =============================================================================


shift = 0.1

def alt_sign(x):
    "Alternate +1/-1 if x is even/odd"
    return (-1) ** x

m1 = aes(x=stage('Level', after_scale='x+shift*alt_sign(x)'))              # shift outward
m2 = aes(x=stage('Level', after_scale='x-shift*alt_sign(x)'), group='gene')  # shift inward


# Violine + points + lines between sub and super
(ggplot(SDRall, aes('Level', 'value', fill = 'Level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.8, show_legend = False)
 + geom_point(m2, color='none', alpha=0.6, size=1.5, show_legend=False)
 + geom_line(m2, color='gray', size=0.65, alpha=0.6)
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR values for 8782 genes (all)')
)

# Violin + points + boxplot
(ggplot(SDRall, aes('Level', 'value', fill = 'Level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR values for 8782 genes (all)')
)

#%% SDRsuper and SDRsub histograms



#%% Example testing 

from plotnine.data import economics

(
 ggplot(economics)              # what data to use
 + aes(x = "date", y = "pop")   # aes() is used to set th varaibles to use for each axis
 + geom_line()                  # Geometric object to use for drawing
 )  

#%% Psuedo data

np.random.seed(123)
n = 20
mu = (1, 2.3)
sigma = (1, 1.6)

before = np.random.normal(loc=mu[0], scale=sigma[0], size=n)
after = np.random.normal(loc=mu[1], scale=sigma[1], size=n)

df = pd.DataFrame({
    'value': np.hstack([before, after]),
    'when': np.repeat(['before', 'after'], n),
    'id': np.hstack([range(n), range(n)])
})

df['when'] = df['when'].astype(pdtypes.CategoricalDtype(categories=['before', 'after']))
df.head()