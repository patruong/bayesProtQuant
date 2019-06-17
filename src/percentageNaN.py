#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:46:20 2019

@author: ptruong
"""

def rawSpecies(df, specie = "Caenorhabditis elegans OX=6239"):
    """
    Slices raw data (500-PSSS3-equ decoy_Report.xls) into 
    """
    if specie == "all":
        species = df
    else:
        species = df[df["PG.Organisms"] == specie]
    return species

def samplePercNaN(df, specie = "Caenorhabditis elegans OX=6239", sample = "S01"):
    """
    Get sample percentage NaNs.
    """
    spec = rawSpecies(df = df, specie = specie)
    nonDecoy = spec[spec["EG.IsDecoy"] == False]    
    sample = sample
    sampleDat = nonDecoy[nonDecoy["R.Condition"].str[-3:] == sample]    
    perc_NaN = sampleDat["PG.Quantity"].isna().sum() / (sampleDat["PG.Quantity"].count() + sampleDat["PG.Quantity"].isna().sum())
    return perc_NaN

def allSamplesPercNaN(df, specie = "Caenorhabditis elegans OX=6239"):
    repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    percNaN_list = []
    for i in repl_cols:
        percNaN = samplePercNaN(df, specie = specie, sample = i)
        percNaN_list.append(percNaN)
    percentage_df = pd.DataFrame(np.array([repl_cols, np.around(percNaN_list, 4)]).T,  columns = ["sample", "percentage_NaN"])
    return percentage_df

def computeRawPercNaNTable(df):
    """
    Compute perc NaN for all species in the raw data (df = pd.read_csv("500-PSSS3-equ decoy_Report.xls", sep ="\t"))
    """
    species = df["PG.Organisms"].unique()
    species = np.append(species, "all")
    percentageList = []
    specieList = []
    for i in species:
        percentages = allSamplesPercNaN(df, specie = i)
        percentageList.append(percentages)
        specieList.append(i)
    
    valuesList = []
    samples = percentageList[0]["sample"]
    for i in range(len(percentageList)):
        valuesList.append(percentageList[i]["percentage_NaN"].values)
    percTable = pd.DataFrame(np.array([samples.values] + valuesList).T, columns = ["sample"] + list(species))
    return percTable


def show_values_on_bars(axs):
    def _show_on_single_plot(ax):        
        for p in ax.patches:
            _x = p.get_x() + p.get_width() / 2
            _y = p.get_y() + p.get_height()
            value = '{:.3f}'.format(p.get_height())
            ax.text(_x, _y, value, ha="center") 

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)
        
        
if __name__ == "__main__":
    df = pd.read_csv("500-PSSS3-equ decoy_Report.xls", sep ="\t")
    
    percTable = computeRawPercNaNTable(df)
    
    meltedPercTable = percTable.melt("sample")
    meltedPercTable = meltedPercTable.convert_objects(convert_numeric=True)
    
            
    fig, ax = plt.subplots(figsize=(13*2.5, 8*2.5))
    #sns.violinplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",split=True)
    sns.barplot(data=meltedPercTable, x="sample",y="value", hue="variable")
    ax.set_ylim(0,1.5)
    ax.set_ylabel("Fraction of NaN")
    #ax1= g.axes[0]
    
    ax.axhline(1, ls='--', color = "r")
    show_values_on_bars(ax)
    
    sns.despine(offset=10, trim=True)
    #plt.title("Homo sapiens")  
    fig.savefig("percentageNaN.png")