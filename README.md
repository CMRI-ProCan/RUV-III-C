# Remove Unwanted Variation III, Complete (RUV-III-C)

Normalisation is an essential step in the analysis of any large omic dataset. Some omics technologies present particular challenges due to the large number of missing measurements. Some of these measurements are missing not at random, where a correct measurement is either zero or below the limit of detection. Others are missing at random, where a correct measurement might be large. As a result, treating missing values as zero creates problems for normalisation, but so does imputing a non-zero value; neither approach is accurate. Treating those values as missing also presents a problem, as most normalisation methods require a complete data matrix. RUV-III-C can normalize multi-omic datasets containing missing values which cannot easily or should not be imputed. 

# Repo Contents

- [data](./data): Example data from Collins, B.C., Hunter, C.L., Liu, Y. *et al.* Multi-laboratory assessment of reproducibility, qualitative and quantitative performance of SWATH-mass spectrometry. *Nat Commun* **8**, 291 (2017). 
- [docker](./docker): Dockerfiles and build scripts related to distribution via docker. 
- [man](./man): Man packages for provided R functions. 
- [R](./R): `R` package code.
- [src](./src): C++ implementations of the core algorithms.
- [tests](./tests): Unit tests

# System Requirements

## Hardware Requirements

There is mimimum requirement to use the `RUV-III-C` package, but you must have sufficient RAM to run the analysis requested. Here ``sufficient'' means enough RAM to load the input data, and also enough for a working set per-thread. This means that memory scales linearly with the number of threads. If you run out of memory, consider lowering the number of threads used. The package will by default use as many threads as there are cores available, although the number of threads can be restricted. 

Four CPU cores and 16GB RAM are much more than is necesarry to run the provided example. 

## Software Requirements

### OS Requirements

This package was developed and tested primarily on Linux. It is compatible with the Linux, Mac and Windows operating systems. 

Several dependency packages must be installed, and a minimum R version of 3.4.0 is suggested. Earlier versions of R are expected to also work, but this is entirely untested. 

# Installation Guide

## Development Version

### Package Dependencies

Users should install the following packages prior to installing `RUVIIIC`, from an `R` terminal:
```
install.packages(c('RSpectra', 'progress', 'Rcpp', 'RcppEigen'))
```

### Package Installation

From an `R` session, type:

```
require(devtools)
install_github('CMRI-ProCan/RUVIIIC')
```

The package should take less than a minute to install.

# Instructions for use

Run this software on your dataset using the same method outline in the following demo. 

# Demo

The following demo runs in a couple of minutes:

```
data(crossLab)
#Design matrix containing information about which runs are technical replicates of each other.
#In this case, random pairings of mass-spec runs analysing the same sample, at different sites.
#Note that we specify no intercept term!
M <- model.matrix(~ grouping - 1, data = peptideData)
#Get out the list of peptides, both HEK (control) and peptides of interest.
peptides <- setdiff(colnames(peptideData), c("filename", "site", "mixture", "Date", "grouping"))
#Reduce the data matrix to only the peptide data
onlyPeptideData <- data.matrix(peptideData[, peptides])
#All the human peptides are potential controls. That is, everything that's not an SIS peptides.
potentialControls <- setdiff(peptides, sisPeptides)
#But we want to use controls that are always found
potentialControlsAlwaysFound <- names(which(apply(onlyPeptideData[, potentialControls], 2, function(x) sum(is.na(x))) == 0))
#Because there are so many potential controls here, we only use 500.
actualControls <- head(potentialControlsAlwaysFound, 500)
#Actually run correction
results <- RUVIII_C(k = 11, Y = onlyPeptideData, M = M, toCorrect = c(sisPeptides, actualControls), controls = actualControls, filename = "results.RData")
```
