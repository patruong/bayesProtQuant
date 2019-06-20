#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 19:16:08 2019

@author: ptruong
"""

from scipy import stats
import scipy as sp
from qvalue import *

#ce_s = ce_s.reset_index()


def computeDiffExp(sample1, sample2):
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
    qVals = qvalues(pVals)
    return qVals

qVals = computeDiffExp(at_s["S02"], at_s["S06"])
print((qVals < 0.05).sum())

qVals2 = computeDiffExp(at_t["S02"], at_t["S06"])
print((qVals2 < 0.05).sum())



##########################

triqler = table2Df("proteins.2vs6.tsv", peptideSeperator = True)
triqler = triqler2Numeric(triqler)
triqler = getTriqlerDetailed(triqler)

triqler.sort_values(by = "protein_id_posterior_error_prob", inplace = True)
triqler["FDR"] = triqler["protein_id_posterior_error_prob"].expanding().mean()
triqler = triqler[triqler["FDR"] < 0.01]

#df_triqler = addSpecieProteinColumns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
df_triqler = triqler
runs = getRuns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
conditionRuns = getConditionRuns(runs)
unique_run = conditionRuns.run.unique()
run_id = getConditionRunId(triqler, unique_run[0])


# All runs
at_t_q = df_triqler[df_triqler["specie"] == "ARATH"]
ce_t_q = df_triqler[df_triqler["specie"] == "CAEEL"]
hs_t_q = df_triqler[df_triqler["specie"] == "HUMAN"]

print((hs_t_q.q_value < 0.05).sum())
print((ce_t_q.q_value < 0.05).sum())
print((at_t_q.q_value < 0.05).sum())

print((hs_t_q['diff_exp_prob_0.1'] < 0.05).sum())
print((ce_t_q['diff_exp_prob_0.1'] < 0.05).sum())
print((at_t_q['diff_exp_prob_0.1'] < 0.05).sum())
