#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:20:07 2019

@author: ptruong
"""

import numpy as np
import pandas as pd
from scipy import stats
import scipy as sp
from qvalue import *

def log2fc(a,b):
    """
    fold change a to b
    """
    fc = np.log2(b) - np.log2(a)
    return fc

####
# FUNCTIONS USING at_t, hs_t, ce_t FORMAT e.g. FOR SAMPLE S01
#
#            S01   S01    S01    S01    S01    
#protein a
#protein b
#protein c
#...    
###
    

def ttest(sample1, sample2):
    """
    t-test differential expression.
    """
    if not (sample1.index.tolist() == sample2.index.tolist()):
        raise Exception("matching error!")
    pVals = []
    for i in range(len(sample1)):
        val1 = sample1.iloc[i].values
        val1 = val1[np.logical_not(np.isnan(val1))]
        val2 = sample2.iloc[i].values
        val2 = val2[np.logical_not(np.isnan(val2))]
        pVal = stats.ttest_ind(val1, val2).pvalue
        pVals.append(pVal)
    pVals = np.array(pVals)
    prots = sample1.index
    pVals = pd.DataFrame(pVals, index = prots)#, columns = ["S02"+"_"+"S06"])
    return pVals

# Compute the log2FC
def get_point_estimate_log2fc(sample1, sample2):
    sample1_mu = sample1.mean(axis = 1)
    sample2_mu = sample2.mean(axis = 1)
    log2fc_vals = log2fc(sample1_mu, sample2_mu)
    return log2fc_vals

# Compute -log2(p)
def compute_log2P(sample1, sample2):
    pVals = ttest(sample1, sample2)
    log2P = -np.log2(pVals)
    return log2P

# Compute -log10(p)
def compute_log10P(sample1, sample2):
    pVals = ttest(sample1, sample2)
    log10P = -np.log10(pVals)
    return log10P

# Compute -log(p)
def compute_logP(sample1, sample2, logFunc = np.log10):
    pVals = ttest(sample1, sample2)
    logP = -logFunc(pVals)
    return logP

# Compute -log2(q)
def compute_log2Q(sample1, sample2):
    pVals = compute_log2P(sample1, sample2)    
    qVals = qvalues(pVals)
    qVals = qVals.apply(pd.to_numeric)
    log2Q = -np.log2(qVals["q"])
    return log2Q
