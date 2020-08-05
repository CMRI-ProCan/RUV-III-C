#' Peptide data from the cross-lab study of Collins et al (2017)
#'
#' A dataset containing data on around 1,000 peptides from the cross-lab study of Collins et al (2017), originally published in Nature Communications. The data consists of a data.frame named \code{peptideData}, and a character vector sisPeptides, listing the peptides that should be analysed. 
#' @format A data.frame containing data on around 1,000 peptides, with peptides that were not detected in a mass-spectrometry acquisition coded as NA. All columns are peptides, except for \code{filename}, \code{site}, \code{mixture}, \code{grouping} and \code{Date}. The peptides include around 1,000 background peptides, and the stable isotope-labelled standard (SIS) peptides which have intensities that change across mass-spectrometry acquisitions. 
#' @usage data(crossLab)
"peptideData"

#' List of stable isotope-labelled standard (SIS) peptides from the cross-lab study of Collins et al (2017)
#' @usage data(crossLab)
"sisPeptides"
