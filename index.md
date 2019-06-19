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

### 2019-06-19 top3 most intense peptides.
![](plots/2019-06-19-top3intense/top3IntensePeptide.png)

PG.Quantity is based on a mix of Top3 intense peptides and reproducibility of identification. What about cases without Top3 intense peptides?

### 2019-06-18 comparison with default settings.
![](plots/2019-06-18_init/H_sapiens.png)
![](plots/2019-06-18_init/A_thaliana.png)
![](plots/2019-06-18_init/C_elegans.png)

Plotted using:

```console
foo@bar:~$ python main.py comparisonPlot --PSSS3_input "../data/500-PSSS3-raw-reformatted_dropna_dropdup_decoy_nonShared_again.csv" --triqler_input "../data/triqler_output_old/proteins.1vs2.tsv" --fdr_treshold 0.05 --outputSuffix TEST

```

## Contact
