# Code for "The Cost of Deviation: A Generalized Spatial Autoregressive Model"

Documented by Jaekyeong Shin

Last updated: November 17, 2023

## Table of contents

- Introduction
- Data cleaning
- Matlab code

## Introduction

The included codes run all the estimations in my job market paper, "The Cost of Deviation: A Generalized Spatial Autoregressive Model." In the paper, I estimate the peer effect with the generalized spatial autoregressive model proposed by myself, using the dataset collected for Banerjee et al. (2013), which I will refer to as BCDJ.

The data can be downloaded here: <https://web.stanford.edu/~jacksonm/Data.html> Unzip the downloaded file (datav4.0.zip) in ```/1_datacleaning```.

BCDJ and the GSAR have two main differences in terms of the data-cleaning process. First, isolated individuals are excluded in BCDJ but not in the GSAR. This will alter the label of individuals in village networks. Second, only the household-level attributes are used in BCDJ, but individual-level data are also used in the GSAR. The purpose of the following process of to address these differences.

## Data cleaning

### Step 1. Run Stata do file

The following two .dta file will be used:

```
/Datasets/datav4.0/Data/2. Demographics and Outcomes/individual_characteristics.dta
/Datasets/datav4.0/Data/2. Demographics and Outcomes/household_characteristics.dta
```

Only the heads of household will be extracted from the individual data and merged into the household data based on the household id (hhid).

Run ```/1_datacleaning/data_cleaning_hh.do```

- Caste

Replace the "missing" variable with "minority." Most of the missing castes are related with Islams and Christians.

Drop the subcastes variable. It will not be used for the estimation.

- Religion

- Categorical variables

Convert categorical variables to dummy variables. The bases are as follow:

1. Caste: Other Backward Castes (OBC).
1. Religion: Hinduism
1. Latrine:
1. Roof type: everything except "tile"
1. Rent status: a house is privately owned
1. Electricity:

BSS entered only into the following 43 villages: (list here). Drop the rest of the villages after running a summary statistics.

Caste data is missing in the first few participating villages. Drop them.

The output file is ```hh_cov_cleaned.csv```. This will be separated into individual villages in the next step.

### Step 2. Run Python script

The Python script will separate the .csv file into the folders labeled with village numbers.

## Matlab code

This code runs estimation importing the cleaned dataset. A total of five specifications will be estimated with the embedded minimizer (fminunc). The code will read the network and characteristics data of each village and store them in cells.

### Appendix 3

The set of codes in This set of codes compares the
