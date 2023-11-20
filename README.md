# Code for "The Cost of Deviation: A Generalized Spatial Autoregressive Model"

Documented by Jaekyeong Shin (<Jaekyeong.Shin@Colorado.edu>)

Last updated: November 17, 2023

## Table of contents

- Introduction
- Data cleaning
- Main result

## Introduction

The included codes run all the estimations in my job market paper, "The Cost of Deviation: A Generalized Spatial Autoregressive Model." In the paper, I estimate the peer effect with the generalized spatial autoregressive model I proposed using the dataset collected for Banerjee et al. (2013), which I will refer to as BCDJ.

Users must download two sets of data from <https://web.stanford.edu/~jacksonm/Data.html>.

1. [Raw data](https://www.stanford.edu/~jacksonm/IndianVillagesDataFiles.zip)
1. [The "Diffusion of Microfinance" study data](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21538)

Extract both .zip files in ```/gsar```, then you should have ```/gsar/2010-0760_Data``` and ```/gsar/datav4.0```. Throughout the rest of the steps, we won't make any changes to these folders. All of the processed data will be stored in ```/gsar/processed```.

BCDJ and the GSAR have two main differences in terms of the data-cleaning process. First, isolated individuals are excluded in BCDJ but not in the GSAR. This will alter the label of individuals in village networks. Second, only the household-level attributes are used in BCDJ, but individual-level data are also used in the GSAR. The purpose of the following process is to address these differences.

## Data cleaning

Note: if you don't want to run the data cleaning yourself, the processed data is available here: [download](https://www.dropbox.com/scl/fi/vwdhukftqfk5kbwdml4xy/processed.zip?rlkey=0shrjdiapzccyvj43myyn093r&dl=0). Unzip the file in ```/)gsar``` and skip to [Main result](#main-result).

### Step 1. Run Stata do file

The following two .dta files will be used:

```
/datav4.0/Data/2. Demographics and Outcomes/individual_characteristics.dta
/datav4.0/Data/2. Demographics and Outcomes/household_characteristics.dta
```

The household characteristics contain caste and religion information. The individual characteristics file is used only to correct those two variables, which could have been incorrectly surveyed in the initial household-level survey.

Run ```/gsar/datacleaning_step1.do```.

- Caste
I replace the missing variables with "minority." Most of them are related to Islam and Christianity. The subcastes variables are also dropped. It won't be used for the estimation.

- Categorical variables
Convert categorical variables to dummy variables. The bases are as follows:

1. Caste: Other Backward Castes (OBC).
1. Religion: Hinduism
1. Latrine: other than no access
1. Roof type: other than "tile"
1. Rent status: a house is privately owned or shared
1. Electricity: other than no access

BSS entered only 43 villages. I drop the rest of the villages after running a summary statistics.

Caste data is missing in the first few participating villages. Those observations are also dropped.

The output file is ```hh_cov_cleaned.csv```. This will be separated into individual villages with the network data in the next step.

### Step 2. Run Python script

The script consists of two parts: first, matching the outcome variable (microfinance take-up) to the adjacency matrix provided by BCDJ and, second, reconstructing the

The adjacency matrices provided by BCDJ do not contain the household ID, and isolated households are excluded. In Step 1, however, I included the isolated nodes and dropped some observations missing the caste information. Therefore, the adjacency matrices must be reconstructed accordingly so that the dimensions are matched to the population of each village and the places of the nodes are to the hhid in the covariates.

**Note: no changes are made to the covariate file in Step 2. Therefore, if users want to try different sets of covariates, it must be done in the previous step.**

This Python script will separate the .csv file into the folders labeled with village numbers.

The two required files must be placed in ```/gsar```:

```
network_list.txt: the filenames of each network provided by BCDJ.
vil_used.txt: the list of 43 villages used.
```

These two files are already included in the repository, so no further action is needed if users have not removed them.

For most of the users, it won't be necessary to modify this portion of the code. You can change all of the specifications a researcher would want to in Step 1 or the next section with the main Matlab code.

Run ```/gsar/datacleaning_step2.py``` with your preferred environment for Python. The processed dataset will be stored in ```/gsar/output```.

## Main Result

This Matlab code runs estimation by importing the cleaned dataset. The six specifications will be estimated with the embedded minimizer (fminunc). The code will read the network and characteristics data of each village and store them in cells.

Note: Parallel Computing Toolbox will be used for some highly repetitive tasks. While it can significantly reduce the computation time for most of the modern multi-core processors, it will utilize more memory to allocate the task to each thread, depending on the number of cores. If you don't think it is necessary, replace _parfor_ with _for_.

Run ```/gsar/main_estimation.m```. The code has been tested with Matlab R2023b.
At the end of each run, the result will be recorded in ```/gsar/results.csv```.

The benchmark logit model (Model 1) is also used as the initial point for the SAR models.

### Appendix 3

In Appendix 3 of my job market paper, I compare the methodology and implication of BCDJ with mine. Run ```/gsar/Appendix3/main_comparison.m``` for the result.
