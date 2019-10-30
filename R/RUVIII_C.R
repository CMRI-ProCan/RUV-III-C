#' RUV-III-C
#'
#' Apply RUV-III-C, a variation of RUV-III that only uses non-missing values
#' 
#' RUV-III is a sophisticated method for removing unwanted variation. The key difficulty in removing unwanted variation is distinguishing wanted from unwanted variation. RUV-III solves this by relying on technical replication, and a list of variables (known as negative control variables) which are known a priori to be constant across all observations. Any variation in the negative control variables across the dataset is (by assumption) unwanted. So we can distinguish wanted from unwanted variation, and therefore estimate the unwanted variation and remove it. 
#' 
#' One problem with this approach is the presence of ``missing'' or zero values in certain application domains. For example, in proteomics it will sometimes be the case that a protein or peptide is not detected in a specific technical replicate of a sample, for purely technical reasons relating to data collection. These missing values are often not related to censoring or the limit of detection. Similar problems occur in metabolomics and single-cell transcriptomics. In all these cases, the metabolite, gene or peptide will be recorded as a zero in the data matrix. Where this type of variation occurs between technical replicates (e.g. one records a zero value and one records a non-zero value) is not correctable.  
#'
#' Regardless of the reason for these zeros, and whether they are accurate or not, zero values are not affected by technical variation, which breaks an assumption of the RUV-III model. In the case that a zero value is incorrect, more serious problems occur. The discrepancies between a pair of technical replicates due to zero values will appear to be much larger than the discrepancies due to other (correctable) technical factors. RUV-III will attempt to correct for the larger (uncorrectable) discrepancy, and ignore the correctable technical factors. 
#' 
#' RUV-III-C is a variation of RUV-III that attempts to solve this problem, by applying RUV-III separately to every variable. If variable X is being corrected, we take the rows of the data matrix for which X is non-missing. RUV-III is then applied, and the corrected values of X is retained. The corrected values of all other variables are discarded. Note that when we take a subset of the rows of the data matrix, other columns besides X will still have missing values. These values are replaced with zero in order to apply RUV-III. No additional transformation is applied to the input data matrix. If normalization should be applied on the log-scale, then logged data must be input. 
#'
#' There are two implementations of this function, one in C++ and one in R. Select which version using the \code{version} argument, which must be either "CPP" or "R"
#' 
#' The R implementation of this function takes an optional intermediate filename, used to store the results. If the function is run again, the previously computed results are loaded directly from the intermediate file, without checking if the inputs specified to the function are the same as those used to generate the previous intermediate result. To re-run a computation, the intermediate file must be deleted. If the intermediate filename is set to NULL, no intermediate file is used. The C++ implementation never saves an intermediate file. 
#' 
#' @param k The number of factors of unwanted variation to remove
#' @param ruvInputData The input data matrix. Must be a matrix, not a data.frame. It should contain missing (NA) values, rather than zeros. No additional transformation is applied to the input data. 
#' @param M The design matrix containing information about technical replicates. It should not contain an intercept term!
#' @param toCorrect The names of the variables to correct using RUV-III-C
#' @param filename The intermediate file in which to save the results, in the R version. Set to NULL to not use an intermediate file. The C++ version never saves an intermediate file. 
#' @param controls The names of the control variables which are known to be constant across the observations
#' @param withW Should we keep the estimate of the unwanted variation factors W, for each variable which is corrected?
#' @param batchSize How often should we write to the intermediate file? The default of 1000 implies that results are written to file every 1000 variables. 
#' @param version The version of the underlying code to use. Must be either "CPP" or "R"
#'
#' @return If withW = FALSE, returns a matrix. If withW = TRUE, returns a list with entries named \code{newY} and \code{W}.
#'
#' @examples
#' data(crossLab)
#' #Design matrix containing information about which runs are technical replicates of each other. 
#' #In this case, random pairings of mass-spec runs analysing the same sample, at different sites.
#' #Note that we specify no intercept term!
#' M <- model.matrix(~ grouping - 1, data = peptideData)
#' #Get out the list of peptides, both HEK (control) and peptides of interest.
#' peptides <- setdiff(colnames(peptideData), c("filename", "site", "mixture", "Date", "grouping"))
#' #Reduce the data matrix to only the peptide data
#' onlyPeptideData <- data.matrix(peptideData[, peptides])
#' #All the human peptides are potential controls. That is, everything that's not an SIS peptides.
#' potentialControls <- setdiff(peptides, sisPeptides)
#' #But we want to use controls that are always found
#' potentialControlsAlwaysFound <- names(which(apply(onlyPeptideData[, potentialControls], 2, function(x) sum(is.na(x))) == 0))
#' #Because there are so many potential controls here, we only use 500. 
#' actualControls <- head(potentialControlsAlwaysFound, 500)
#' #Actually run correction
#' \dontrun{results <- RUVIII_C_R(k = 11, ruvInputData = onlyPeptideData, M = M, toCorrect = c(sisPeptides, actualControls), controls = actualControls, filename = "results.RData")}
#' @export
RUVIII_C <- function(k, ruvInputData, M, toCorrect, filename, controls, withW = FALSE, batchSize = 1000, version = "CPP")
{
	if(version == "CPP")
	{
		return(RUVIII_C_CPP(k = k, input = ruvInputData, M = M, toCorrect = toCorrect, controls = controls, withW = withW))
	}
	else if(version == "R")
	{
		return(RUVIII_C_R(k = k, ruvInputData = ruvInputData, M = M, toCorrect = toCorrect, filename = filename, controls = controls, withW = withW, batchSize = batchSize))
	}
	else
	{
		stop("version must be either 'CPP' or 'R'")
	}
}
RUVIII_C_R <- function(k, ruvInputData, M, toCorrect, filename, controls, withW = FALSE, batchSize = 1000)
{
	if(missing(filename))
	{
		stop("Function RUVIII_NM_Varying requires a filename for intermediate results, or NULL if an intermediate file should not be used")
	}
	if(k <= 0)
	{
		stop("Input k, the number of factors of unwanted variation, must be positive")
	}
	if(is.factor(controls))
	{
		stop("Input controls was a factor. Please input a vector of variable names")
	}
	if(is.null(rownames(ruvInputData)))
	{
		stop("Function RUVIII_NM_Varying requires row-names for input ruvInputData")
	}
	if(any(apply(M, 2, sum) == nrow(M)))
	{
		warning("Input design matrix had a column of all ones, indicating an intercept term. This is generally an error!")
	}
	if(any(is.na(ruvInputData[, controls])))
	{
		stop("The negative control variables should never be NA")
	}
	if(k >= nrow(ruvInputData))
	{
		stop("Input k cannot be larger than the number of rows in the input matrix")
	}
	if(k >= length(controls))
	{
		stop("Input k cannot be larger than the number of negative controls")
	}
	#Replace NAs with 0
	ruvInputDataWithoutNA <- ruvInputData
	ruvInputDataWithoutNA[is.na(ruvInputDataWithoutNA)] <- 0

	results <- list()
	results$peptideResults <- results$alphaResults <- results$W <- list()

	#Load previous results set. 
	if(!is.null(filename) && file.exists(filename))
	{
		load(filename)
	}
	pb <- progress_bar$new(total = length(toCorrect), format = "[:bar] :percent :eta")
	#Don't bother writing anything to file, if there are no new computations. 
	newComputation <- FALSE
	if(length(results$peptideResults) != length(toCorrect) || any(sort(names(results$peptideResults)) != sort(toCorrect)))
	{
		for(peptideIndex in 1:length(toCorrect))
		{
			peptide <- toCorrect[peptideIndex]
			if(!(peptide %in% names(results$peptideResults)))
			{
				#row indices of the rows which have actual values
				indices <- which(!is.na(ruvInputData[,peptide]))
				#If tere are no non-missing observations, then we can just mark this peptide as uncorrected / uncorrectable
				if(length(indices) > 0)
				{
					#Submatrix of the entire data matrix
					submatrix <- ruvInputDataWithoutNA[indices, , drop = F]
					#Submatrix of the replicate matrix
					Msubset <- M[indices, , drop = F]
					#Some whole biological replicates are removed in the process. So remove those entire columns. This has consequences for the dimension tests later on. 
					Msubset <- Msubset[, apply(Msubset, 2, function(x) sum(x) > 0), drop=F]
					#You need at least two observations, across two different biological samples, in ordor to make any kind of correction
					if(ncol(Msubset) < 2 || nrow(Msubset) < 2) 
					{
						#Mark this peptide as uncorrected / uncorrectable.
						newComputation <- TRUE
						results$alphaResults[peptide] <- list(c())
						results$peptideResults[peptide] <- list(c())
						if(withW) results$W[peptide] <- list(c())
						next
					}

					#Copied from the RUVIII code
					Msubset <- replicate.matrix(Msubset)
					#Residual with respect to Msubset
					Y0 <- ruv::residop(submatrix, Msubset)
					m <- nrow(submatrix)
					#If the dimensions don't work, we can't correct this variable
					if(min(m - ncol(Msubset), length(controls)) < k)
					{
						newComputation <- TRUE
						results$alphaResults[peptide] <- list(c())
						results$peptideResults[peptide] <- list(c())
						if(withW) results$W[peptide] <- list(c())
						next
					}
					try({
						#Now the RUV-III code. This may throw exceptions, possibly for numerical reasons, hence the try / catch. 
						eigenDecomp <- eigs_sym(Y0 %*% t(Y0), k = min(m - ncol(Msubset), length(controls)), which = "LM")
						fullalpha <- t(eigenDecomp$vectors) %*% submatrix
						alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
						ac <- alpha[,controls, drop = FALSE]
						currentPeptideW <- submatrix[, controls] %*% t(ac) %*% solve(ac %*% t(ac))
						currentResult <- fullalpha[1:k, peptide]

						#Now we extend the corrected results out to a full matrix. 
						adjusted <- vector(mode = "numeric", length = nrow(ruvInputDataWithoutNA))
						adjusted[] <- NA
						adjusted[indices] <- submatrix[,peptide, drop=F] - currentPeptideW %*% cbind(currentResult)
						names(adjusted) <- rownames(ruvInputDataWithoutNA)

						results$alphaResults[[peptide]] <- currentResult
						results$peptideResults[[peptide]] <- adjusted
						if(withW)
						{
							results$W[[peptide]] <- currentPeptideW
						}
						newComputation <- TRUE
					}, silent = TRUE)
				} else
				{
					newComputation <- TRUE
					results$alphaResults[peptide] <- list(c())
					results$peptideResults[peptide] <- list(c())
					if(withW) results$W[peptide] <- list(c())
				}
			}
			#This case only triggers if there is an exception in the try() above. 
			if(!(peptide %in% names(results$peptideResults)))
			{
				newComputation <- TRUE
				results$alphaResults[peptide] <- list(c())
				results$peptideResults[peptide] <- list(c())
				if(withW) results$W[peptide] <- list(c())
			}
			#If we've hit the batch size, and we've actually made a new computation, write out to the temporary file. 
			if(peptideIndex %% batchSize == 0 && newComputation && !is.null(filename))
			{
				save(results, file = filename)
			}
			pb$tick()
		}
		#Exit the progress bar.
		pb$terminate()
		#If we've made any new computations write out to file. 
		if(newComputation && !is.null(filename)) save(results, file = filename)
	}
	if(withW)
	{
		return(list(newY = do.call(cbind, results$peptideResults), W = results$W))
	}
	else
	{
		return(do.call(cbind, results$peptideResults))
	}
}
