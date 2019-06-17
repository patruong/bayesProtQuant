#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:23:03 2019

@author: ptruong
"""

import pandas as pd

df = pd.read_csv("500-PSSS3-equ decoy_Report.xls", sep ="\t")

# FDR treshold
FDR_treshold = 0.01
df = df[df["PG.Qvalue"].fillna(-1)<FDR_treshold]

FDR_treshold = 0.01
df = df[df["EG.Qvalue"].fillna(-1)<FDR_treshold]

# Cscore treshold
Cscore_treshold = 0.06 
df = df[df["PG.Cscore"].fillna(1)>Cscore_treshold]

# Cscore treshold
Cscore_treshold = 0.06 
df = df[df["EG.Cscore"].fillna(1)>Cscore_treshold]

#df[df["EG.IsDecoy"] == False]

df.to_csv("500-PSSS3-equ decoy_Report_filteredOn_PGQvalue_EGQvalue_PGCscore_EGCscore.xls", sep = "\t")

f = open("500-PSSS3-equ decoy_Report.xls", "r")
emp