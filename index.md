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

Also here is a link to a page [page2](page2)
## Project description
## Tools

## Project Log
### 2020-06-30 Comparison plot for adjusted intensity and parallell plot for non-adjusted intensity.

plots/2019-06-30_compPlot/full/A thaliana_plot.png
![](plots/2019-06-30_compPlot/full/A thaliana_plot.png)

#### full set
![](plots/2019-06-30_compPlot/full/A thaliana_plot.png)
![](plots/2019-06-30_compPlot/full/C elegans_plot.png)
![](plots/2019-06-30_compPlot/full/H sapiens_plot.png)

#### S02-S08
![](plots/2019-06-30_compPlot/S02-S08/A thaliana_plot.png)
![](plots/2019-06-30_compPlot/S02-S08/C elegans_plot.png)
![](plots/2019-06-30_compPlot/S02-S08/H sapiens_plot.png)

#### S03-S08
![](plots/2019-06-30_compPlot/S03-S08/A thaliana_plot.png)
![](plots/2019-06-30_compPlot/S03-S08/C elegans_plot.png)
![](plots/2019-06-30_compPlot/S03-S08/H sapiens_plot.png)

#### Parallell coordinates of non-adjusted intensity
Will do this on adjusted data set later...

[A. Thaliana parallell coordinates](pages/at_t.html)
[C. Elegans parallell coordinates](pages/ce.html)
[H. Sapiens parallell coordinates](pages/hs.html)


### 2020-06-27 fixing the 1.0 intensity error in PSSS3 data set for triqler.

[Here is a link to a short analysis of the 1.0 intensities](pages/check_1_intensity.md)

ToDo: if not "decoy" suffix in triqlerFormatted data: if intensity == 0 => Drop by setting the intensities to 0.0.

### 2020-06-26 loglog plot of p-values vs q-values of triqler results of human S02 v S06

![](results/2019-06-26_loglogpq/pi0_est_hs_t_s02_s06.png)
![](results/2019-06-26_loglogpq/pi0_1_hs_t_s02_s06.png)
![](results/2019-06-26_loglogpq/hs_t_s02_s06_histogram.png)


Converting p-values to q-values should usually lower the amount of significant proteins. Looking at the q-value conversion with pi0 estimates, we can see that at high p- and q-values, there are actually more significant proteins for p-values. The histogram shows that this is because the pi0 estimate is set very low (pi0 = 0.0724). Fixing the pi0 = 1 gives the expected behaviour of less significant proteins at q-value significance.

Higher pi0 gives higher q-values and are thus more conservatives and should give less significance.

These experiment can be generated using code in https://github.com/patruong/bayesProtQuant/tree/master/results/2019-06-26_loglogpq (commit  03f3714)



### 2020-06-20 Differential Expression
Parameters:
S02 vs S06
Triqler parameters: minSamp = 15

**Homo Sapiens (S02vS06)**

| Spectronaut (t-test) | Triqler (t-test) | Triqler (FC_treshold = 0.1) | Triqler (FC_treshold = 1.0)|Triqler (FC_treshold = 1.2) |
|----------------------|------------------|--------------------------|--------------------------|--------------------------| 
| 4711                 | 6350             | 5887                     |   2577                       |   346                       |

**Caenorhabditis elegans (S02vS06)**

| Spectronaut (t-test) | Triqler (t-test) | Triqler (FC_treshold = 0.1) | Triqler (FC_treshold = 1.0)|Triqler (FC_treshold = 1.2) |
|----------------------|------------------|--------------------------|--------------------------|--------------------------| 
| 213                 | 2461             | 2162                     |   886                       |     821                     |

**Arabidopsis thaliana (S02vS06)**

| Spectronaut (t-test) | Triqler (t-test) | Triqler (FC_treshold = 0.1) | Triqler (FC_treshold = 1.0)|Triqler (FC_treshold = 1.2) |
|----------------------|------------------|--------------------------|--------------------------|--------------------------| 
| 39                 | 22             | 606                     |     0                     |      0                    |

Fixed bug in triqlerParser.py causing only half the triqler data sets to be read.

```python
def readFile(filename):
    """
    read in filename and return header and rows of the file as lists.
    """       
    f = open(filename, "r")
    f_header = f.readline().replace("\n","").split("\t")
    f_rows = []
    for i in f:
        f_row = f.readline().replace("\n","").split("\t") #python 3.6 # <---- ERROR
        f_rows.append(f_row)
    return f_header, f_rows
```
Should be 

```python
...
    for i in f:
        f_row = i.replace("\n","").split("\t") #python 3.6 # <---- ERROR
        f_rows.append(f_row)
...
```

New comparison plots below:

![](plots/2019-06-20_fixedParserBug/HS_plot_minSamp15_logFold1.png)
![](plots/2019-06-20_fixedParserBug/CE_plot_minSamp15_logFold1.png)
![](plots/2019-06-20_fixedParserBug/AT_plot_minSamp15_logFold1.png)


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
