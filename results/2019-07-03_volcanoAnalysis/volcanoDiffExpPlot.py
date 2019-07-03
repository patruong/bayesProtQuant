#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 11:40:31 2019

@author: ptruong
"""

import pandas as pd
import seaborn as sns
from parseVolcanoResults import *
import matplotlib.pyplot as plt


df = volcanoResultsDict_to_meltedDf(filename = "volcanoResultsDict.txt")

def barplotComps(meltedDf, specie, title, figname = None):    
    """
    Bar diff expression comparison plot between methods.
    """
    import matplotlib.pyplot as plt
    a4_dims = (11.7*1.5, 8.27*1.5)
    #df = mylib.load_data()
    fig, ax = plt.subplots(figsize=a4_dims)
    plt.xticks(rotation=90)
    sns.barplot(x="comparison", y="diffExp", hue="method", data=meltedDf[meltedDf["specie"] == specie], ax = ax).set_title(title)
    if figname:
        fig.savefig(figname)

specie = "ARATH"
barplotComps(df, specie, specie, specie+"_diffExp.png")
