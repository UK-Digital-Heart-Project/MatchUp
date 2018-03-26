# MatchUp
Matching multiple groups of subjects by covariates. 

## Requisites
R (and RStudio).
No packages needed.

## Installation
Open MatchUp.R  
Enter ```getwd()``` in the console to find which folder you are in.
Put `Matching.csv` in this folder.
This contains some sample data shown below - in this instance with 3 disease groups.

| Unique ID | Group | Sequence number | AGE | SBP | BSA  | SEX | 
|-----------|-------|-----------------|-----|-----|------|-----| 
| 0114/0014 | 1     | 1               | 45  | 134 | 1.23 | 1   | 
| 0115/0061 | 1     | 2               | 43  | 123 | 1.31 | 1   | 
| 0115/0086 | 1     | 3               | 56  | 143 | 1.02 | 0   | 
| 0117/0014 | 2     | 4               | 23  | 123 | 1.34 | 1   | 
| 0117/0016 | 2     | 5               | 33  | 132 | 1.33 | 1   | 
| 0117/0018 | 2     | 6               | 78  | 142 | 1.56 | 0   | 
| 0117/0022 | 3     | 7               | 44  | 111 | 1.44 | 0   | 
| 0117/0023 | 3     | 8               | 54  | 99  | 1.33 | 0   | 
| 0217-/006 | 3     | 9               | 33  | 166 | 1.09 | 1   | 

Run the function by typing: ```MatchUp("Matching.csv", final.group.size = ...)```  
The output is the same CSV file with one extra column appended which denotes whether that subject should be included or not in the matched dataset.

Edit `Matching.csv` as necessary / change `final.group.size` as required.

## Notes
`Matching.csv` must be (a) .csv format, (b) column names for each column, (c) all covariate names in capitals / no non-covariates in capitals, (d) ID columnname = 'Unique.ID' (NB cases)

`Final.group.size` decides the number in each final group. Small samples may lead to non-convergence in the GLM, although the warnings for this have been switched off


## Assumptions:

1. There are three or more non-overlapping clinical groups which need to be matched on the basis of covariates

2. The initial groups may be of unequal sizes

3. You want to match the groups 1:1:1

4. There are no missing data



## Design:

i. Propensity matching approach

ii. Sequential comparison of each group ('target' group) with the remaining groups ('pooled' groups)

iii. Subject in the target group with lowest propensity score is discarded, provided the number of subjects in the target group exceeds a minimum threshold

