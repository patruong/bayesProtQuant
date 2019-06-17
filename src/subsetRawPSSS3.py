#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 14:33:49 2019

@author: ptruong
"""

def subsetPSSS3bySample(rawPSSS3 = "500-PSSS3-equ decoy_Report_nonShared_20190507.xls",
                        sampleList = ["S01", "S02"],
                        filename = "output.xls"):
    df = pd.read_csv("500-PSSS3-equ decoy_Report_nonShared_20190507.xls", sep ="\t")
    df = df[np.isin(df["R.Condition"].str[-3:], ["S02", "S03", "S04", "S10"])]
    df.to_csv("PSSS3_noShared_S02_S03_S04_S10.xls")
    