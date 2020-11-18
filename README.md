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

There is no minimum requirement to use the `RUV-III-C` package, but you must have sufficient RAM to run the analysis requested. Here ``sufficient'' means enough RAM to load the input data, and also enough for a working set per-thread. This means that memory scales linearly with the number of threads. If you run out of memory, consider lowering the number of threads used. The package will by default use as many threads as there are cores available, although the number of threads can be restricted. 

Four CPU cores and 16GB RAM are much more than is necessary to run the provided example. 

## Software Requirements

### OS Requirements

This package was developed and tested primarily on Linux. It is compatible with the Linux, Mac and Windows operating systems.

Several dependency packages must be installed, and a minimum R version of 3.4.0 is suggested. Earlier versions of R are expected to also work, but this is entirely untested. 

On Linux the package is configured to use the system Basic Linear Algebra Subprogram (BLAS) library, making performance on Linux much better than on Windows. The requirement for a system BLAS library can be dropped by removing the flag `-DEIGEN_USE_BLAS` from src/Makevars, however this will remove the performance improvement compared to Windows. 

# Installation Guide

## Package Dependencies

This package requires:
* R 3.5 or later
* R packages Rcpp, RSpectra and RcppEigen

Users should install the dependency packages prior to installing `RUVIIIC`, from an `R` terminal:
```
install.packages(c('RSpectra', 'progress', 'Rcpp', 'RcppEigen'))
```

The package has been tested with the following versions of these dependencies:
```
R 3.6.3
RSpectra_0.14-0
progress_1.2.0
Rcpp_1.0.1
RcppEigen_0.3.3.5.0
```
However, it is expected that all newer versions of these packages will also work. 

## Package Installation

From an `R` session, type:

```
require(devtools)
install_github('CMRI-ProCan/RUV-III-C')
```
Alternatively, from the console run:
```
git clone https://github.com/CMRI-ProCan/RUV-III-C.git RUV-III-C
R CMD INSTALL RUV-III-C
```

The package should take less than a minute to install.

# Example

We use data from the multi-laboratory proteomics study of Collins et al. (2017). This study analysed the same dilution series of spiked-in peptides in a background of human peptides, at eleven laboratories around the world.
```
data(crossLab)
```
We create a design matrix `M` containing information about which runs are technical replicates of each other. In this case, we have already specified random pairings of mass-spec runs analysing the same sample, at different sites. Note that the `- 1` specifies no intercept term!
```
M <- model.matrix(~ grouping - 1, data = peptideData)
```
Get out the list of all peptides (human and spiked-in). 
```
peptides <- setdiff(colnames(peptideData), c("filename", "site", "mixture", "Date", "grouping"))
```
Reduce the data matrix to only the peptide data.
```
onlyPeptideData <- data.matrix(peptideData[, peptides])
```
All the human peptides are potential negative control variables. That is, everything that's not a spiked-in peptide. But we want to use negative control variables that are found in every sample.
```
potentialControls <- setdiff(peptides, sisPeptides)
potentialControlsAlwaysFound <- names(which(apply(onlyPeptideData[, potentialControls], 2, function(x) sum(is.na(x))) == 0))
```
Because there are so many potential controls here, we only use 500.
```
actualControls <- head(potentialControlsAlwaysFound, 500)
```
Now we actually run the normalization using RUV-III-C. We use a very high value of `k` because we have a large number of very high quality negative control variables in this study. Non-dilution datasets will likely use a smaller value of `k`. 
```
results <- RUVIII_C(k = 11, Y = onlyPeptideData, M = M, toCorrect = c(sisPeptides, actualControls), controls = actualControls)
```
Alternatively, we can use controls that are _not_ found in every sample. 
```
actualControls <- head(potentialControls, 1000)
results <- RUVIII_C_Varying(k = 11, Y = onlyPeptideData, M = M, toCorrect = c(sisPeptides, actualControls), controls = actualControls)
```

# Credits
Rohan Shah, Sean Peters, Qing Zhong. 


