
# Analysis of FG.NormalizedMS2PeakArea = 1

Information regarding FG.NormalizedMS2PeakArea = 1.0 from Roland.

Values of 1 can come from two sources:

In Spectronaut, small values for quantities (<1) are set to one, these arise from small noise peaks or from local normalization effects.

In both cases were the signals noise or close to noise.

It mostly arises due to the fact that the dynamic range of MS1 and MS2 are not necessarily the same.

In Spectronaut MS1 and MS2 information is used for identification and it can be that one layer is enough. So the quantitative information of the other layer can be very low.

Import relevant modules.


```python
import pandas as pd 
from utils import *
```

Load the pickled raw data set.


```python
df = pd.read_pickle("../../data/500-PSSS3-equ decoy_Report_nonShared_20190507.pkl")

```

Load triqler formatted data set. (Formatted from the same data set).


```python
df_triqler = pd.read_pickle("../../data/PSSS3_triqlerFormatted_nonShared.pkl")

```

Find proteins where intensity (FG.NormalizedMS2PeakArea) = 1.


```python
proteins = df_triqler[df_triqler["intensity"] == 1].proteins
proteins = proteins[proteins.str[:5] != "decoy"]
proteins = proteins.str.split("_", expand = True)
proteins = proteins.rename(index = str, columns = {0:"protein", 1:"specie"})
unique_proteins = pd.Series(proteins.protein.unique())


```

Show these proteins.


```python
unique_proteins
```




    0        A0A023T4K3
    1        A0A061ACK4
    2        A0A061ACL3
    3        A0A061ACR1
    4        A0A061ACU2
    5        A0A061ACU6
    6        A0A061ACY0
    7        A0A061AD21
    8        A0A061AD39
    9        A0A061AD47
    10       A0A061AE05
    11       A0A061AJK8
    12       A0A061AKY5
    13       A0A061AL89
    14       A0A078BPH9
    15       A0A078BPJ4
    16       A0A078BPM1
    17       A0A0H3W5N0
    18       A0A0K3AQM0
    19       A0A0K3AQN5
    20       A0A0K3AQN9
    21       A0A0K3AQV8
    22       A0A0K3AQY3
    23       A0A0K3AR10
    24       A0A0K3ARF3
    25       A0A0K3ARM6
    26       A0A0K3ARN9
    27       A0A0K3ARR9
    28       A0A0K3ARY0
    29       A0A0K3ARZ7
                ...    
    20675        P06732
    20676        P22352
    20677        P27216
    20678        P55082
    20679        P56381
    20680        Q13233
    20681      Q13555-6
    20682        Q15291
    20683        Q16611
    20684        Q19185
    20685        Q1G3M2
    20686      Q6PJT7-9
    20687        Q7XA86
    20688        Q84UC7
    20689        Q8L4R0
    20690        Q8N8J7
    20691        Q8NG11
    20692        Q8VYK9
    20693        Q9BS34
    20694        Q9H1E5
    20695        Q9LFA4
    20696        Q9MBH1
    20697        Q9NWT8
    20698        Q9STL2
    20699        Q9SZY4
    20700        Q9UAX1
    20701        Q9UBM7
    20702      Q9UPN9-2
    20703        Q9XIM0
    20704        Q9ZVA5
    Length: 20705, dtype: object




```python
query_prot = unique_proteins[0]
query = df[df["PG.ProteinAccessions"] == query_prot]

```


```python
query
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>R.Condition</th>
      <th>R.FileName</th>
      <th>PG.Organisms</th>
      <th>PG.ProteinAccessions</th>
      <th>PG.Cscore</th>
      <th>PG.NrOfStrippedSequencesIdentified</th>
      <th>PG.Qvalue</th>
      <th>PG.Quantity</th>
      <th>EG.StrippedSequence</th>
      <th>EG.IsDecoy</th>
      <th>EG.PrecursorId</th>
      <th>EG.PEP</th>
      <th>EG.Qvalue</th>
      <th>EG.Cscore</th>
      <th>FG.NormalizedMS2PeakArea</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>SFYYLVQDLK</td>
      <td>False</td>
      <td>_SFYYLVQDLK_.2</td>
      <td>1.479293e-02</td>
      <td>0.003878</td>
      <td>2.355977</td>
      <td>7.319080e+05</td>
    </tr>
    <tr>
      <th>1</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>2.210408e-02</td>
      <td>0.005114</td>
      <td>2.073540</td>
      <td>2.299625e+04</td>
    </tr>
    <tr>
      <th>2</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.2</td>
      <td>1.000000e+00</td>
      <td>0.762912</td>
      <td>-2.174942</td>
      <td>4.274121e+04</td>
    </tr>
    <tr>
      <th>3</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>AADFYVR</td>
      <td>False</td>
      <td>_AADFYVR_.2</td>
      <td>9.465682e-03</td>
      <td>0.015673</td>
      <td>1.283862</td>
      <td>3.631738e+04</td>
    </tr>
    <tr>
      <th>4</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>IGALADVNNSKDPDGLR</td>
      <td>False</td>
      <td>_IGALADVNNSKDPDGLR_.3</td>
      <td>4.225811e-01</td>
      <td>0.102376</td>
      <td>0.036027</td>
      <td>2.253101e+04</td>
    </tr>
    <tr>
      <th>5</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>IGALADVNNSKDPDGLR</td>
      <td>False</td>
      <td>_IGALADVNNSKDPDGLR_.2</td>
      <td>1.000000e+00</td>
      <td>0.760221</td>
      <td>-2.151609</td>
      <td>2.341229e+04</td>
    </tr>
    <tr>
      <th>6</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.2</td>
      <td>4.517035e-01</td>
      <td>0.117293</td>
      <td>-0.067812</td>
      <td>1.782665e+04</td>
    </tr>
    <tr>
      <th>7</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.3</td>
      <td>6.079661e-01</td>
      <td>0.219510</td>
      <td>-0.606611</td>
      <td>8.547810e+04</td>
    </tr>
    <tr>
      <th>8</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>9.981995e-01</td>
      <td>0.702184</td>
      <td>-1.814023</td>
      <td>8.090994e+03</td>
    </tr>
    <tr>
      <th>9</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.3</td>
      <td>9.884009e-01</td>
      <td>0.693141</td>
      <td>-1.773851</td>
      <td>2.529676e+04</td>
    </tr>
    <tr>
      <th>10</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000e+00</td>
      <td>0.584152</td>
      <td>-1.466288</td>
      <td>1.000000e+00</td>
    </tr>
    <tr>
      <th>11</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>QELEILYK</td>
      <td>False</td>
      <td>_QELEILYK_.2</td>
      <td>5.828639e-01</td>
      <td>0.268762</td>
      <td>-0.809965</td>
      <td>1.876812e+05</td>
    </tr>
    <tr>
      <th>12</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>IIEDSEIMQEDDDNWPEPDKIGR</td>
      <td>False</td>
      <td>_IIEDSEIMQEDDDNWPEPDKIGR_.3</td>
      <td>7.583292e-01</td>
      <td>0.375077</td>
      <td>-1.108799</td>
      <td>1.882612e+04</td>
    </tr>
    <tr>
      <th>13</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>EATVSESVLSELK</td>
      <td>False</td>
      <td>_EATVSESVLSELK_.2</td>
      <td>1.000000e+00</td>
      <td>0.758897</td>
      <td>-2.142972</td>
      <td>4.276881e+04</td>
    </tr>
    <tr>
      <th>14</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>IGALADVNNSK</td>
      <td>False</td>
      <td>_IGALADVNNSK_.2</td>
      <td>9.720947e-01</td>
      <td>0.773937</td>
      <td>-2.279790</td>
      <td>1.754513e+04</td>
    </tr>
    <tr>
      <th>15</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>IIEDSEIMQEDDDNWPEPDK</td>
      <td>False</td>
      <td>_IIEDSEIMQEDDDNWPEPDK_.2</td>
      <td>1.000000e+00</td>
      <td>0.736726</td>
      <td>-1.985881</td>
      <td>2.863420e+03</td>
    </tr>
    <tr>
      <th>16</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>QELEILYKNEHISFTTGK</td>
      <td>False</td>
      <td>_QELEILYKNEHISFTTGK_.4</td>
      <td>1.000000e+00</td>
      <td>0.766065</td>
      <td>-2.205700</td>
      <td>4.767215e+05</td>
    </tr>
    <tr>
      <th>17</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.15625</td>
      <td>QELEILYKNEHISFTTGK</td>
      <td>False</td>
      <td>_QELEILYKNEHISFTTGK_.3</td>
      <td>1.000000e+00</td>
      <td>0.742685</td>
      <td>-2.022538</td>
      <td>1.182180e+04</td>
    </tr>
    <tr>
      <th>18</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NSEIGFTTHK</td>
      <td>True</td>
      <td>_NSEIGFTTHK_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-0.886623</td>
      <td>3.543555e+04</td>
    </tr>
    <tr>
      <th>19</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NSEIGFTTHK</td>
      <td>True</td>
      <td>_NSEIGFTTHK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.928180</td>
      <td>5.298748e+04</td>
    </tr>
    <tr>
      <th>20</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>VADFYAR</td>
      <td>True</td>
      <td>_VADFYAR_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.296488</td>
      <td>3.413020e+04</td>
    </tr>
    <tr>
      <th>21</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SAVVLEETSESLKR</td>
      <td>True</td>
      <td>_SAVVLEETSESLKR_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.588380</td>
      <td>1.039064e+04</td>
    </tr>
    <tr>
      <th>22</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SAVVLEETSESLKR</td>
      <td>True</td>
      <td>_SAVVLEETSESLKR_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.314687</td>
      <td>7.577624e+05</td>
    </tr>
    <tr>
      <th>23</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>FGHENLFFLEREFGSPR</td>
      <td>True</td>
      <td>_FGHENLFFLEREFGSPR_.4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.283839</td>
      <td>6.206880e+02</td>
    </tr>
    <tr>
      <th>24</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>FGHENLFFLEREFGSPR</td>
      <td>True</td>
      <td>_FGHENLFFLEREFGSPR_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.757634</td>
      <td>7.628955e+03</td>
    </tr>
    <tr>
      <th>25</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>EELQYLIK</td>
      <td>True</td>
      <td>_EELQYLIK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.280042</td>
      <td>4.756001e+04</td>
    </tr>
    <tr>
      <th>26</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>QEEEILLKNIYHSFTTGK</td>
      <td>True</td>
      <td>_QEEEILLKNIYHSFTTGK_.4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.222821</td>
      <td>1.015951e+05</td>
    </tr>
    <tr>
      <th>27</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>QEEEILLKNIYHSFTTGK</td>
      <td>True</td>
      <td>_QEEEILLKNIYHSFTTGK_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-0.119591</td>
      <td>7.760912e+04</td>
    </tr>
    <tr>
      <th>28</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>GKFFHEFSEGEFRPNGLLR</td>
      <td>True</td>
      <td>_GKFFHEFSEGEFRPNGLLR_.4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.395250</td>
      <td>3.577354e+03</td>
    </tr>
    <tr>
      <th>29</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>VGNLANIDASK</td>
      <td>True</td>
      <td>_VGNLANIDASK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.343939</td>
      <td>4.672448e+04</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>52155410</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.2</td>
      <td>3.332570e-06</td>
      <td>0.000132</td>
      <td>2.548576</td>
      <td>1.445004e+05</td>
    </tr>
    <tr>
      <th>52155411</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.3</td>
      <td>1.681911e-03</td>
      <td>0.000467</td>
      <td>2.016687</td>
      <td>1.781270e+04</td>
    </tr>
    <tr>
      <th>52155412</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>2.766263e-09</td>
      <td>0.000060</td>
      <td>2.865193</td>
      <td>6.142457e+04</td>
    </tr>
    <tr>
      <th>52155413</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.3</td>
      <td>1.977536e-10</td>
      <td>0.000057</td>
      <td>2.962228</td>
      <td>8.107272e+04</td>
    </tr>
    <tr>
      <th>52155414</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>2.559424e-09</td>
      <td>0.000060</td>
      <td>2.868159</td>
      <td>2.840809e+04</td>
    </tr>
    <tr>
      <th>52155415</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>QELEILYK</td>
      <td>False</td>
      <td>_QELEILYK_.2</td>
      <td>3.783442e-05</td>
      <td>0.000191</td>
      <td>2.401054</td>
      <td>2.886927e+04</td>
    </tr>
    <tr>
      <th>52155416</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>IIEDSEIMQEDDDNWPEPDKIGR</td>
      <td>False</td>
      <td>_IIEDSEIMQEDDDNWPEPDKIGR_.3</td>
      <td>3.861019e-08</td>
      <td>0.000075</td>
      <td>2.759681</td>
      <td>2.008109e+04</td>
    </tr>
    <tr>
      <th>52155417</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>EATVSESVLSELK</td>
      <td>False</td>
      <td>_EATVSESVLSELK_.2</td>
      <td>3.696719e-02</td>
      <td>0.003889</td>
      <td>0.774860</td>
      <td>2.714550e+04</td>
    </tr>
    <tr>
      <th>52155418</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>IGALADVNNSK</td>
      <td>False</td>
      <td>_IGALADVNNSK_.2</td>
      <td>9.931060e-03</td>
      <td>0.001268</td>
      <td>1.462434</td>
      <td>3.666480e+04</td>
    </tr>
    <tr>
      <th>52155419</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>IIEDSEIMQEDDDNWPEPDK</td>
      <td>False</td>
      <td>_IIEDSEIMQEDDDNWPEPDK_.2</td>
      <td>3.624652e-01</td>
      <td>0.051595</td>
      <td>-1.179870</td>
      <td>8.187363e+03</td>
    </tr>
    <tr>
      <th>52155420</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>QELEILYKNEHISFTTGK</td>
      <td>False</td>
      <td>_QELEILYKNEHISFTTGK_.4</td>
      <td>6.368914e-01</td>
      <td>0.109106</td>
      <td>-1.759308</td>
      <td>4.004741e+05</td>
    </tr>
    <tr>
      <th>52155421</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>11</td>
      <td>0.000646</td>
      <td>190713.62500</td>
      <td>QELEILYKNEHISFTTGK</td>
      <td>False</td>
      <td>_QELEILYKNEHISFTTGK_.3</td>
      <td>1.611248e-01</td>
      <td>0.018855</td>
      <td>-0.347197</td>
      <td>3.279107e+04</td>
    </tr>
    <tr>
      <th>52155422</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NSEIGFTTHK</td>
      <td>True</td>
      <td>_NSEIGFTTHK_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.625697</td>
      <td>8.006405e+04</td>
    </tr>
    <tr>
      <th>52155423</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NSEIGFTTHK</td>
      <td>True</td>
      <td>_NSEIGFTTHK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.655590</td>
      <td>1.387353e+04</td>
    </tr>
    <tr>
      <th>52155424</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>VADFYAR</td>
      <td>True</td>
      <td>_VADFYAR_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.541574</td>
      <td>6.440874e+04</td>
    </tr>
    <tr>
      <th>52155425</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SAVVLEETSESLKR</td>
      <td>True</td>
      <td>_SAVVLEETSESLKR_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.266281</td>
      <td>2.088234e+04</td>
    </tr>
    <tr>
      <th>52155426</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SAVVLEETSESLKR</td>
      <td>True</td>
      <td>_SAVVLEETSESLKR_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.667791</td>
      <td>1.217765e+05</td>
    </tr>
    <tr>
      <th>52155427</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>FGHENLFFLEREFGSPR</td>
      <td>True</td>
      <td>_FGHENLFFLEREFGSPR_.4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.190037</td>
      <td>1.841496e+03</td>
    </tr>
    <tr>
      <th>52155428</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>FGHENLFFLEREFGSPR</td>
      <td>True</td>
      <td>_FGHENLFFLEREFGSPR_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.551600</td>
      <td>1.682023e+04</td>
    </tr>
    <tr>
      <th>52155429</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>EELQYLIK</td>
      <td>True</td>
      <td>_EELQYLIK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.859780</td>
      <td>1.620637e+03</td>
    </tr>
    <tr>
      <th>52155430</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>QEEEILLKNIYHSFTTGK</td>
      <td>True</td>
      <td>_QEEEILLKNIYHSFTTGK_.4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.595791</td>
      <td>1.904184e+04</td>
    </tr>
    <tr>
      <th>52155431</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>QEEEILLKNIYHSFTTGK</td>
      <td>True</td>
      <td>_QEEEILLKNIYHSFTTGK_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-0.315062</td>
      <td>1.340976e+06</td>
    </tr>
    <tr>
      <th>52155432</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>GKFFHEFSEGEFRPNGLLR</td>
      <td>True</td>
      <td>_GKFFHEFSEGEFRPNGLLR_.4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.517616</td>
      <td>1.472539e+03</td>
    </tr>
    <tr>
      <th>52155433</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>VGNLANIDASK</td>
      <td>True</td>
      <td>_VGNLANIDASK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.681926</td>
      <td>6.866233e+03</td>
    </tr>
    <tr>
      <th>52155434</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>VDALADPNNSKIGDGLR</td>
      <td>True</td>
      <td>_VDALADPNNSKIGDGLR_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.507371</td>
      <td>1.992760e+04</td>
    </tr>
    <tr>
      <th>52155435</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>VDALADPNNSKIGDGLR</td>
      <td>True</td>
      <td>_VDALADPNNSKIGDGLR_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.980533</td>
      <td>3.229269e+04</td>
    </tr>
    <tr>
      <th>52155436</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>LQYYLVSDFK</td>
      <td>True</td>
      <td>_LQYYLVSDFK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.511397</td>
      <td>2.168764e+05</td>
    </tr>
    <tr>
      <th>52155437</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>LTAVSESELSVEK</td>
      <td>True</td>
      <td>_LTAVSESELSVEK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.657463</td>
      <td>1.782522e+05</td>
    </tr>
    <tr>
      <th>52155438</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>IIDDSEMEWEDDDNQPEPKIGIR</td>
      <td>True</td>
      <td>_IIDDSEMEWEDDDNQPEPKIGIR_.3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.562214</td>
      <td>5.152218e+04</td>
    </tr>
    <tr>
      <th>52155439</th>
      <td>S500-PSSS3-S01</td>
      <td>G_D180330_S500-PSSS3-S01_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>MIIDSEDIEEDDQPWNEPDK</td>
      <td>True</td>
      <td>_MIIDSEDIEEDDQPWNEPDK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-2.688025</td>
      <td>5.576035e+03</td>
    </tr>
  </tbody>
</table>
<p>1800 rows × 15 columns</p>
</div>



A lot of these proteins in fact has peptide intensities. Lets investigate where intensities are 0.


```python
query[query["FG.NormalizedMS2PeakArea"] == 1]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>R.Condition</th>
      <th>R.FileName</th>
      <th>PG.Organisms</th>
      <th>PG.ProteinAccessions</th>
      <th>PG.Cscore</th>
      <th>PG.NrOfStrippedSequencesIdentified</th>
      <th>PG.Qvalue</th>
      <th>PG.Quantity</th>
      <th>EG.StrippedSequence</th>
      <th>EG.IsDecoy</th>
      <th>EG.PrecursorId</th>
      <th>EG.PEP</th>
      <th>EG.Qvalue</th>
      <th>EG.Cscore</th>
      <th>FG.NormalizedMS2PeakArea</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>10</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>2</td>
      <td>0.000646</td>
      <td>377452.156250</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.584152</td>
      <td>-1.466288</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>2128802</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>0.971182</td>
      <td>0.877314</td>
      <td>-1.756694</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>3193196</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.731580</td>
      <td>-1.305529</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>3193202</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>IGALADVNNSK</td>
      <td>False</td>
      <td>_IGALADVNNSK_.2</td>
      <td>1.000000</td>
      <td>0.913514</td>
      <td>-2.087202</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>4257592</th>
      <td>S500-PSSS3-S10</td>
      <td>G_D180330_S500-PSSS3-S10_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.882508</td>
      <td>-1.837320</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>5321988</th>
      <td>S500-PSSS3-S08</td>
      <td>G_D180330_S500-PSSS3-S08_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.749338</td>
      <td>-1.102499</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>7450780</th>
      <td>S500-PSSS3-S08</td>
      <td>G_D180330_S500-PSSS3-S08_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>0.949748</td>
      <td>0.900993</td>
      <td>-1.883789</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>7450795</th>
      <td>S500-PSSS3-S08</td>
      <td>G_D180330_S500-PSSS3-S08_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>FGHENLFFLEREFGSPR</td>
      <td>True</td>
      <td>_FGHENLFFLEREFGSPR_.4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.485949</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>9579572</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.744226</td>
      <td>-1.351521</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>9579573</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.3</td>
      <td>1.000000</td>
      <td>0.646199</td>
      <td>-1.209836</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>10643968</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>0.982366</td>
      <td>0.922066</td>
      <td>-2.057751</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>13837149</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>0.978127</td>
      <td>0.832960</td>
      <td>-2.124984</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>13837156</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.749722</td>
      <td>-1.772931</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>14901545</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>1.000000</td>
      <td>0.948016</td>
      <td>1.891271</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>14901554</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R04_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.948983</td>
      <td>1.806365</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>17030343</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.3</td>
      <td>0.948465</td>
      <td>0.832369</td>
      <td>-1.417118</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>17030344</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.801456</td>
      <td>-1.343402</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>18094733</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>1.000000</td>
      <td>0.779697</td>
      <td>-2.380094</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>18094742</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>0.956955</td>
      <td>0.758794</td>
      <td>-2.127311</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>19159131</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>AADFYVR</td>
      <td>False</td>
      <td>_AADFYVR_.2</td>
      <td>0.950289</td>
      <td>0.733853</td>
      <td>-1.887230</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>19159135</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.3</td>
      <td>0.937416</td>
      <td>0.726019</td>
      <td>-1.865410</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>20223532</th>
      <td>S500-PSSS3-S08</td>
      <td>G_D180330_S500-PSSS3-S08_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.889777</td>
      <td>-1.878795</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>21287921</th>
      <td>S500-PSSS3-S10</td>
      <td>G_D180330_S500-PSSS3-S10_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>1.000000</td>
      <td>0.923936</td>
      <td>-1.457002</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>21287923</th>
      <td>S500-PSSS3-S10</td>
      <td>G_D180330_S500-PSSS3-S10_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>AADFYVR</td>
      <td>False</td>
      <td>_AADFYVR_.2</td>
      <td>0.945463</td>
      <td>0.957276</td>
      <td>-1.718887</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>21287927</th>
      <td>S500-PSSS3-S10</td>
      <td>G_D180330_S500-PSSS3-S10_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.3</td>
      <td>0.964483</td>
      <td>0.940092</td>
      <td>-1.562498</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>21287928</th>
      <td>S500-PSSS3-S10</td>
      <td>G_D180330_S500-PSSS3-S10_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>0.820509</td>
      <td>0.656198</td>
      <td>-0.931232</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>22352330</th>
      <td>S500-PSSS3-S03</td>
      <td>G_D180330_S500-PSSS3-S03_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>3</td>
      <td>0.000646</td>
      <td>224372.781250</td>
      <td>IGALADVNNSK</td>
      <td>False</td>
      <td>_IGALADVNNSK_.2</td>
      <td>0.823999</td>
      <td>0.376484</td>
      <td>-1.475100</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>23416720</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>1</td>
      <td>0.000646</td>
      <td>830202.687500</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.703447</td>
      <td>-1.285880</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>23416722</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>1</td>
      <td>0.000646</td>
      <td>830202.687500</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.915586</td>
      <td>-2.226666</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>26609901</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>0.923761</td>
      <td>0.922159</td>
      <td>-1.667343</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>26609910</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>0.948613</td>
      <td>0.948630</td>
      <td>-1.890678</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>26609912</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R05_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>IIEDSEIMQEDDDNWPEPDKIGR</td>
      <td>False</td>
      <td>_IIEDSEIMQEDDDNWPEPDKIGR_.3</td>
      <td>1.000000</td>
      <td>0.769534</td>
      <td>-1.065896</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>29803091</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>AADFYVR</td>
      <td>False</td>
      <td>_AADFYVR_.2</td>
      <td>1.000000</td>
      <td>0.946887</td>
      <td>-1.492815</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>29803095</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R03_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.3</td>
      <td>1.000000</td>
      <td>0.888488</td>
      <td>-1.251080</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>30867485</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>0.942934</td>
      <td>0.912783</td>
      <td>-1.570019</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>30867487</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>AADFYVR</td>
      <td>False</td>
      <td>_AADFYVR_.2</td>
      <td>0.852968</td>
      <td>0.897701</td>
      <td>-1.457832</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>30867492</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.926072</td>
      <td>-1.735365</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>30867493</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.3</td>
      <td>1.000000</td>
      <td>0.841734</td>
      <td>-1.236124</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>31931881</th>
      <td>S500-PSSS3-S09</td>
      <td>G_D180330_S500-PSSS3-S09_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>1.000000</td>
      <td>0.832340</td>
      <td>-1.661479</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>34060675</th>
      <td>S500-PSSS3-S08</td>
      <td>G_D180330_S500-PSSS3-S08_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>AADFYVR</td>
      <td>False</td>
      <td>_AADFYVR_.2</td>
      <td>0.956659</td>
      <td>0.952245</td>
      <td>-1.919342</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>34060680</th>
      <td>S500-PSSS3-S08</td>
      <td>G_D180330_S500-PSSS3-S08_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>0.906760</td>
      <td>0.926398</td>
      <td>-1.619079</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>34060697</th>
      <td>S500-PSSS3-S08</td>
      <td>G_D180330_S500-PSSS3-S08_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>EELQYLIK</td>
      <td>True</td>
      <td>_EELQYLIK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-0.977456</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>35125069</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>0.962650</td>
      <td>0.896980</td>
      <td>-1.603957</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>35125071</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>AADFYVR</td>
      <td>False</td>
      <td>_AADFYVR_.2</td>
      <td>0.959751</td>
      <td>0.873907</td>
      <td>-1.456922</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>36189465</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>0.947165</td>
      <td>0.865787</td>
      <td>-1.914139</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>36189472</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>FGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_FGHEFLEFEFRPNGSLR_.4</td>
      <td>1.000000</td>
      <td>0.894660</td>
      <td>-2.123795</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>36189489</th>
      <td>S500-PSSS3-S07</td>
      <td>G_D180330_S500-PSSS3-S07_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>EELQYLIK</td>
      <td>True</td>
      <td>_EELQYLIK_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.575424</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>37253861</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>1.000000</td>
      <td>0.687519</td>
      <td>-1.132478</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>37253874</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>IGALADVNNSK</td>
      <td>False</td>
      <td>_IGALADVNNSK_.2</td>
      <td>0.875654</td>
      <td>0.499425</td>
      <td>-0.858331</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>38318257</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>1.000000</td>
      <td>0.900043</td>
      <td>-2.153550</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>38318266</th>
      <td>S500-PSSS3-S06</td>
      <td>G_D180330_S500-PSSS3-S06_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>0.979426</td>
      <td>0.832022</td>
      <td>-1.763272</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>39382655</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>1</td>
      <td>0.000646</td>
      <td>735648.250000</td>
      <td>AADFYVR</td>
      <td>False</td>
      <td>_AADFYVR_.2</td>
      <td>0.889480</td>
      <td>0.733052</td>
      <td>-1.420900</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>40447049</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>EATVSESVLSELKR</td>
      <td>False</td>
      <td>_EATVSESVLSELKR_.3</td>
      <td>0.945992</td>
      <td>0.616186</td>
      <td>-1.354987</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>40447055</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>NEHISFTTGK</td>
      <td>False</td>
      <td>_NEHISFTTGK_.3</td>
      <td>0.921434</td>
      <td>0.785328</td>
      <td>-1.711489</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>40447059</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>QELEILYK</td>
      <td>False</td>
      <td>_QELEILYK_.2</td>
      <td>1.000000</td>
      <td>0.854983</td>
      <td>-2.048480</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>40447079</th>
      <td>S500-PSSS3-S05</td>
      <td>G_D180330_S500-PSSS3-S05_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>NaN</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>VDALADPNNSKIGDGLR</td>
      <td>True</td>
      <td>_VDALADPNNSKIGDGLR_.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-1.610883</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>41511458</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>1</td>
      <td>0.000646</td>
      <td>667178.000000</td>
      <td>IGALADVNNSK</td>
      <td>False</td>
      <td>_IGALADVNNSK_.2</td>
      <td>0.989736</td>
      <td>0.502940</td>
      <td>-1.490796</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>42575850</th>
      <td>S500-PSSS3-S04</td>
      <td>G_D180330_S500-PSSS3-S04_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>0</td>
      <td>0.000646</td>
      <td>NaN</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>0.950841</td>
      <td>0.643791</td>
      <td>-1.965023</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>44704642</th>
      <td>S500-PSSS3-S03</td>
      <td>G_D180330_S500-PSSS3-S03_MHRM_R02_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>5</td>
      <td>0.000646</td>
      <td>246017.718750</td>
      <td>GKFGHEFLEFEFRPNGSLR</td>
      <td>False</td>
      <td>_GKFGHEFLEFEFRPNGSLR_.4</td>
      <td>0.995765</td>
      <td>0.584585</td>
      <td>-2.150846</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>48962232</th>
      <td>S500-PSSS3-S02</td>
      <td>G_D180330_S500-PSSS3-S02_MHRM_R01_T0</td>
      <td>Caenorhabditis elegans OX=6239</td>
      <td>A0A023T4K3</td>
      <td>0.182511</td>
      <td>8</td>
      <td>0.000646</td>
      <td>97127.828125</td>
      <td>QELEILYKNEHISFTTGK</td>
      <td>False</td>
      <td>_QELEILYKNEHISFTTGK_.4</td>
      <td>0.997542</td>
      <td>0.368788</td>
      <td>-2.107628</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>
</div>



Most of the PG.Quantity values are NaN. But, looking at above table we can see that for e.g. S04:R03 We have PG.Quantity = 377452, even though the FG.NormalizedMS2PeakArea = 0. This is because we have multiple peptide sequences (EG.StrippedSequence) which belong to the protein A0A023T4K3. Spectronaut uses a combination of top 3 intensities combined and Cscore to quantify proteins. 

In cases where all samples and runs have 1.0 intensities the PG.Quantity is simply registered as NaN, meaning that the 1.0 intensities do not induce any bias to any downstream analysis. These low intensities are therefore ignored in PG.Quantity.


The 1.0 intensities are not handled in the triqler data parsing module (Code below is not executed as it is only to show snippets from the triqler data parsing module).


```python

from __future__ import print_function

import sys
import numpy as np
import math
import csv
import os
import itertools
import re
from collections import defaultdict, namedtuple


def getTsvReader(filename):
  # Python 3
  if sys.version_info[0] >= 3:
    return csv.reader(open(filename, 'r', newline = ''), delimiter = '\t')
  # Python 2
  else:
    return csv.reader(open(filename, 'rb'), delimiter = '\t')

TriqlerSimpleInputRowHeaders = "run condition charge searchScore intensity peptide proteins".split(" ")
TriqlerInputRowHeaders = "run condition charge spectrumId linkPEP featureClusterId searchScore intensity peptide proteins".split(" ")
TriqlerInputRowBase = namedtuple("TriqlerInputRow", TriqlerInputRowHeaders)

class TriqlerInputRow(TriqlerInputRowBase):
  def toList(self):
    l = list(self)
    return l[:-1] + l[-1]
  
  def toSimpleList(self):
    l = list(self)
    return l[:3] + l[6:-1] + l[-1]
  
  def toString(self):
    return "\t".join(map(str, self.toList()))

def parseTriqlerInputFile(triqlerInputFile):
  reader = getTsvReader(triqlerInputFile)
  headers = next(reader)
  hasLinkPEPs = "linkPEP" in headers
  getUniqueProteins = lambda x : list(set([p for p in x if len(p.strip()) > 0]))
  intensityCol = 7 if hasLinkPEPs else 4
  seenPeptChargePairs = dict()
  for i, row in enumerate(reader):
    if i % 1000000 == 0:
      print("  Reading row", i)
    
    intensity = float(row[intensityCol])
    if intensity > 0.0: # <-------------------- 1.0 intensities passes through triqler 
      if hasLinkPEPs:
        proteins = getUniqueProteins(row[9:])
        yield TriqlerInputRow(row[0], row[1], int(row[2]), int(row[3]), float(row[4]), int(row[5]), float(row[6]), intensity, row[8], proteins)
      else:
        key = (int(row[2]), row[5])
        if key not in seenPeptChargePairs:
          seenPeptChargePairs[key] = len(seenPeptChargePairs)
        proteins = getUniqueProteins(row[6:])
        yield TriqlerInputRow(row[0], row[1], int(row[2]), (i+1) * 100, 0.0, seenPeptChargePairs[key], float(row[3]), intensity, row[5], proteins)

```

Note the line "if intensity > 0.0:" allowing 1.0 intensities to be passed in triqler. We could set the non decoy intensities to 0.0 and let triqler filter the 1.0 intensities, since they are also ignored in spectronaut.


```python
print("triqler input rows: " + str(len(df_triqler)))

# I don't think decoy intensities are used...
#print("triqler input rows with intensity = 1" + str(len(df_triqler[df_triqler["intensity"] == 1])))
#print("ratio of noise rows to all rows: " + str(len(df_triqler[df_triqler["intensity"] == 1])/len(df_triqler)))

# Therefore it should be
print("triqler input rows with intensity = 1" + str(len(proteins)))
print("ratio of noise rows to all rows: " + str(len(proteins)/len(df_triqler)))
```

    triqler input rows: 52354084
    triqler input rows with intensity = 11131973
    ratio of noise rows to all rows: 0.021621484199780862


2.16% of the input is noise. I am not certain if this is a large or small number in this case.
