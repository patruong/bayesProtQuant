#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 16:11:26 2019

@author: ptruong
"""
"""
log2fc
http://rstudio-pubs-static.s3.amazonaws.com/13988_bb11d85b79b2436280de434988558140.html
https://www.quora.com/How-do-you-interpret-Log-Fold-Change-with-Gene-Expression
"""

from proc import *
import numpy as np

def log2fc(a,b):
    """
    fold change a to b
    """
    fc = np.log2(b) - np.log2(a)
    return fc

def get_sample_values(mix, s1, s2):
    """
    s1 and s2 are integers between 1-10 (for samples S01,..., S10)
    mix - mixture from proc.getMixtures()
    """
    samples = np.array([s1,s2])-1    
    comps = []
    for sample in samples:
        comps.append(mix[sample])
    return comps

def get_fc_between_samples(s1,s2):
    """
    Function gets log2fc between samples for all species.
    s1 and s2 are integers between 1-10 ( for samples S01, ..., S10)
    e.g.
    s1 = 2
    s2 = 6
    """
    a_mix, c_mix, h_mix, mix = getMixtures()
    samples = ["S01","S02","S03","S04","S05",
               "S06","S07","S08","S09","S10"]
    mixes = [a_mix, c_mix, h_mix]
    species = ["ARATH", "CAEEL", "HUMAN"]
    fc_dict = dict()
    print("Computing log2fc("+samples[s1-1] +","+samples[s2-1]+") for all species:")
    for i in range(len(mixes)):
        comps = (get_sample_values(mixes[i], s1, s2))
        fc = log2fc(comps[0], comps[1])
        fc_dict[species[i]] = fc
    return fc_dict

if __name__ == "__main__":
    fc_dict = get_fc_between_samples(2,6)

    