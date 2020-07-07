# RUV-III-C: Application of RUV-III, using only complete (non-missing) data

RUV-III is a normalization method based on replication, factor analysis and negative control variables (variables that are known _a priori_ to be affected only by unwanted variation). RUV-III-C extends this method to the case where some values in the dataset are missing or not recorded. This is important for biological data generated using a number of technologies (e.g. mass-spectrometry based proteomics). 

# Requirements
* R 3.4 or greater
* R packages Rcpp, RSpectra and RcppEigen

# Installation
For the latest stable release from github, run
``` 
devtools::install_github(repo = "CMRI-ProCan/RUV-III-C")
```
or
```
git clone https://github.com/CMRI-ProCan/RUV-III-C.git RUV-III-C
R CMD INSTALL RUV-III-C
```

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
Rohan Shah, Qing Zhong. 
