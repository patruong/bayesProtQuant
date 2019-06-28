#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 11:40:37 2019

@author: ptruong
"""

import pandas as pd

df = pd.read_pickle("PSSS3_triqlerFormatted_nonShared.pkl")
#df["intenstiy"] = df.loc[(df.proteins.str[:5] != "decoy") & (df["intensity"] == 1)]["intensity"] = 0
df["intensity"] = df.loc[(df.proteins.str[:5] != "decoy") & (df["intensity"] == 1)]["intensity"] = 0
df.loc[(df.proteins.str[:5] != "decoy") & (df["intensity"] == 1), "intensity"] = 0

df.to_csv("PSSS3_triqlerFormatted_nonShared_nonDecoyIntensity1Adj.csv", sep = "\t", index = False)



