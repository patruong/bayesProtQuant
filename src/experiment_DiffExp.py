#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:50:38 2019

@author: ptruong
"""
from triqlerParser import *
from triqlerProcessor import * 
from utils import * 
from parseReformattedPSSS3equDecoy import * # Needs to tidy and place functions in right place!
from corrPlot import * #correlation plotting
from qvalue import *
from diffExp import *

import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy as sp
import subprocess
%matplotlib inline
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})#, 'grid.color': '.1', 'grid.linestyle': u'-'})


if __name__ == "__main__":
    comparisons = getComparisons(fileDir = "triqlerOutput_noShared_largeScaleOptimizations_minSamp15", filePrefix = "proteins", fileSuffix = ".tsv")
    triqlerComps  = getTriqlerComparisons(comparisons, prefix = "proteins", suffix = "tsv")
    spectronautComps = getSpectronautComparison(comparisons)
    
    
    
    # Input directory of triqler outputs    
    diffExpDir = "triqlerOutput_noShared_largeScaleOptimizations_minSamp15"
    FDRTreshold = 0.01
    specie = None
    diffExpComps = processComparisons(triqlerOutputDir = diffExpDir,
                       PSSS3 = df_psss3,
                       comparisons = comparisons,
                       triqlerComps = triqlerComps,
                       spectronautComps = spectronautComps,
                       FDRTreshold = FDRTreshold,
                       specie = specie)
    all_df = meltCompDf(comp_df)
    barplotComps(all_df, specie)
