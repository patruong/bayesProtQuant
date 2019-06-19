# bayesProtQuant

# Table of Contents
1. [Introduction](#introduction)
2. [Project Description](#project-description)
3. [Tools](#tools)
4. [Project Log](#project-log)
5. [Contact](#contact)

## Introduction


[About]()
[Members]()
[Project Log]()

Also here is a link to a page [page2](page_2)
## Project description
## Tools

## Project Log
### 2020-06-20 Differential Expression
Parameters:
S02 vs S06
Triqler parameters: minSamp = 15, log_fold_change = 0.1

Homo Sapiens (S02vS06)


| Spectronaut   | Triqler | Triqler  |
| ------------- | ------------- | ------------- |
| Content Cell  | Content Cell  | Content Cell  |
| Content Cell  | Content Cell  | Content Cell  |

Homo Sapiens (S02vS06)

|Homo Sapiens (S02vS06)|
| Spectronaut (t-test) | Triqler (t-test) | Triqler (reported q-val) |
|----------------------|------------------|--------------------------|
| 4711                 | 6350             | 5887                     |

Caenorhabditis elegans (S02vS06)
| Spectronaut (t-test) | Triqler (t-test) | Triqler (reported q-val) |
|----------------------|------------------|--------------------------|
| 213                  | 2461             | 2161                     |

Arabidopsis thaliana (S02vS06)
| Spectronaut (t-test) | Triqler (t-test) | Triqler (reported q-val) |
|----------------------|------------------|--------------------------|
| 4711                 | 6350             | 5887                     |



### 2019-06-19 top3 most intense peptides.
![](plots/2019-06-19-top3intense/top3IntensePeptide.png)

PG.Quantity is based on a mix of Top3 intense peptides and reproducibility of identification. What about cases without Top3 intense peptides? And how does it account for reproducibility of identification? I guess this should be the C-score (searchScore) in triqler output format.

![](plots/2019-06-19-top3intense/trashGetsIntoTriqler.png)

If we filter on top3 peptides then some trash, like the 1.0 values seen above gets into Triqler. I need to check the searchScore for these to see how much they impact the posteriors.

![](plots/2019-06-19-top3intense/top3IntensitiesHomoSapSpectronaut.png)

Normalized with protein intensity / sum of proteins intensities. Just taking top 3 intensities does not give a good match at all. This confirms that we are sending a lot of trash into triqler. I am still uncertain if triqler searchScore someone manages to adjust this properly.

### 2019-06-18 comparison with default settings.
![](plots/2019-06-18_init/H_sapiens.png)
![](plots/2019-06-18_init/A_thaliana.png)
![](plots/2019-06-18_init/C_elegans.png)

Plotted using:

```console
foo@bar:~$ python main.py comparisonPlot --PSSS3_input "../data/500-PSSS3-raw-reformatted_dropna_dropdup_decoy_nonShared_again.csv" --triqler_input "../data/triqler_output_old/proteins.1vs2.tsv" --fdr_treshold 0.05 --outputSuffix TEST

```

## Contact
