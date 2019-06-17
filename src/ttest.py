#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:44:45 2019

@author: ptruong
"""


"""
alpha = 0.05
tTable = stats.t.ppf(1-(alpha/2), 10) # two-tail t distbution table


hs_s["S03"].mean(axis = 1)
hs_s["S05"].mean(axis = 1)

diff = hs_s["S03"].values - hs_s["S05"].values

muDiff = hs_s["S03"].mean(axis = 1 ) - hs_s["S05"].mean(axis = 1) 
total_sum = muDiff.sum()
squaredDiff = np.power(muDiff, 2)
sse = squaredDiff.sum() 

total_sum = -73
sse = 1131
N = 11

N = len(muDiff)
tValue = (total_sum / N) / (np.sqrt((sse - (np.power(total_sum, 2)/N))/((N-1)*N)))

diffExp = None
if tValue > tTable: #P-value is less than alpha level
    print("Reject the null hypothesis that there is no difference between means - DIFFERENTIALLY EXPRESSED")
    diffExp = True
else:
    print("Cannot Reject null hypothesis - THERE IS NO DIFFERENCE")
    diffExp = False
"""
##########################################################
from scipy import stats
import numpy as np

#values1 = hs_s["S03"]
#values2 = hs_s["S05"]

def computeTvalues(values1, values2):
    """
    Computed student t-test t-values given two data series.
    
    Follow error means we have nan comparisons:
    __main__:2: RuntimeWarning: invalid value encountered in greater
    """     
    diff = values1.values - values2.values
    diffSquared = (values1.values - values2.values)**2
    N = np.shape(diff)[1]
    
    total_sum = (diff.sum(axis = 1))
    sse = (diffSquared.sum(axis = 1))
    
    tValues = (total_sum / N) / (np.sqrt((sse - (np.power(total_sum, 2)/N))/((N-1)*N)))
    return tValues

def checkDiffExp(tValues, alpha = 0.05, twoTailed = True):
    #tValues = computeTvalues(values1, values2)
    alpha = alpha
    twoTailed = twoTailed
    if twoTailed == True:
        tTable = stats.t.ppf(1-(alpha/2), 10) # two-tail t distbution table
    else:
        tTable = stats.t.ppf(1-(alpha), 10) # one-tail t distbution table
    
    nanIdx = np.where(np.isnan(tValues))[0] # Find indices of NaNs
    #print("abs version")
    featDiff = np.array((abs(tValues) > tTable), dtype = object) # Will set nan comparisons to False
    featDiff[nanIdx] = np.NaN # Adjust NaN comparisons to NaN
    return featDiff

def diffExp(values1, values2, alpha = 0.05, twoTailed = True):
    tValues = computeTvalues(values1, values2)
    featDiff = checkDiffExp(tValues, alpha = alpha, twoTailed = twoTailed)
    return featDiff

#diffExp(values1, values2)














