#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:16:52 2019

@author: ptruong
"""
# CORRELATION PLOTTING FUNCTIONS.
from scipy.stats.stats import pearsonr
from numpy import log10



def my_corr(row, mix):
    a_mix = np.array([0.5] * 10)
    a_mix[0] += 0.0001
    c_mix = np.array([0.5, 0.25, 0.125, 0.0625, 0.031, 0.0155, 0.008, 0.004, 0.002, 0.0])
    h_mix = np.array([0.0, 0.25, 0.375, 0.4375, 0.469, 0.4845, 0.492, 0.496, 0.498, 0.5])
    a_mix = a_mix/a_mix.sum()
    c_mix = c_mix/c_mix.sum()
    h_mix = h_mix/h_mix.sum()
    mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}
    mymix = mix[row['Species']]
    row = row.loc[repl_cols]
    ixs = (row != 0.0)
    a = np.array(row.loc[ixs])
    b = mymix[ixs]
    return pearsonr(a,b)[0]

def p_r_plot(species,figname=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    N=5000
    for meth in unmeltsamp['Method'].unique():
        sns.distplot(unmeltsamp[(unmeltsamp['Method']==meth) & (unmeltsamp['Species']==species)]['Correlation, R'],
             hist_kws={"cumulative":-1, "linewidth": 3, "alpha": 1,"histtype": "step"},
             bins=N, kde=False, ax=ax, label=meth, norm_hist = True)
    ax.set_xlim(-0.9,1.0)
    ax.set_ylim(0,1.2)
    ax.set_title(species)
    sns.despine(offset=10, trim=True)
    ax.set_ylabel("Number of proteins")
    ax.set_xlabel("Cumulative samples "r"$\leq$" " correlation coefficient")
    plt.legend()
    if figname:
        fig.savefig(figname)