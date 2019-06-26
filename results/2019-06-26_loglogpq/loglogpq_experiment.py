#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:13:58 2019

@author: ptruong
"""
from func import *
from proc import *


def qvalues(pvalues, pi0):
    """
    Takes p-value dataframe as input and output q-value dataframe.
    """

    pcolname = pvalues.columns[0]

    m = float(len(pvalues.transpose().values[0].tolist()))
    assert(m>0)
    pvalues = pvalues.sort_values(by = pcolname)
    
    if pi0 == 1:
        pi0 = 1
    else:
        pi0 = estimatePi0(pvalues.transpose().values[0].tolist())
    print("pi0: " + str(pi0))
    num_p, p_sum = 0, 0.0

    qs = pd.DataFrame(columns = [ 'q' ])

    for index, row in pvalues.iterrows():
        p = row[pcolname]
        num_p += 1
        p_sum += p
        q = pi0*p*m/float(num_p)
        qs.loc[index,'q'] = q 


    qs = qs.iloc[::-1]
    old_q=1.0
    for ix in range(len(qs)):
        q = min(old_q,qs.iloc[ix, 0])
        old_q = q
        qs.iloc[ix, 0] = q
    return qs, pi0


def use_triqler_samples(s1, s2):
    """
    s1, s2 = S01, ..., S10.
    """
    triqlerFile = "tmp"
    protein_id_fdr_treshold = 0.01
    # Read data
    triqler = read_pickle_Triqler(triqlerFile, protein_id_fdr_treshold)
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = False, truncated = False)
    
    #define samples
    sample1 = hs_t[s1]
    sample2 = hs_t[s2]
    return sample1, sample2

def use_spectronaut_samples(s1,s2):
        
    # Read data
    spectronautFile = "tmp"
    protein_id_fdr_treshold = 0.01
    
    spectronaut =  read_pickle_Spectronaut(spectronautFile, protein_id_fdr_treshold, impute = None, global_impute = True)
    at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = False)
    
    #define samples
    sample1 = hs_s[s1]
    sample2 = hs_s[s2]
    return sample1, sample2

#####################
# plot pi0 estimate #
#####################
def plot_loglog_pq(sample1, sample2, pi0, title):
    pVals = ttest(sample1, sample2)
    qVals, pi0 = qvalues(pVals, pi0)
    qVals = qVals.apply(pd.to_numeric)
    f, ax = plt.subplots(figsize=(6, 6))
    ax.loglog(pVals.sort_values(by = 0, ascending = False), qVals)
    ax.plot([0,1],[0,1])
    ax.set_xlabel("p-values")
    ax.set_ylabel("q-values")
    ax.set_title("loglog plot, " + title + " pi0 = " + str(pi0))

##################
# plot histogram #
##################

def plot_pVals_hist(sample1, sample2, n_bins = 100):
    pVals = ttest(sample1, sample2)
    pVals.hist(bins = n_bins)

if __name__ == "__main__":
    sample1, sample2 = use_triqler_samples("S02", "S06")
    #sample1, sample2 = use_spectronaut_samples("S02", "S06")
    plot_loglog_pq(sample1, sample2, pi0 = 1, title = "Human S02 v S06 triqler")
    plot_loglog_pq(sample1, sample2, pi0 = "est", title = "Human S02 v S06 triqler")
    plot_pVals_hist(sample1, sample2, n_bins = 100)



