#' @include RcppExports.R
#' @title RUV-III-C
#'
#' @description Apply RUV-III-C, a variation of RUV-III that only uses non-missing values
#' 
#' @details RUV-III is a sophisticated method for removing unwanted variation. The key difficulty in removing unwanted variation is distinguishing wanted from unwanted variation. RUV-III solves this by relying on technical replication, and a list of variables (known as negative control variables) which are known a priori to be constant across all observations. Any variation in the negative control variables across the dataset is (by assumption) unwanted. So we can distinguish wanted from unwanted variation, and therefore estimate the unwanted variation and remove it. 
#' 
#' One problem with this approach is the presence of ``missing'' or zero values in certain application domains. For example, in proteomics it will sometimes be the case that a protein or peptide is not detected in a specific technical replicate of a sample, for purely technical reasons relating to data collection. These missing values are often not related to censoring or the limit of detection. Similar problems occur in metabolomics and single-cell transcriptomics. In all these cases, the metabolite, gene or peptide will be recorded as a zero in the data matrix. Where this type of variation occurs between technical replicates (e.g. one records a zero value and one records a non-zero value) is not correctable.  
#'
#' Regardless of the reason for these zeros, and whether they are accurate or not, zero values are not affected by technical variation, which breaks an assumption of the RUV-III model. In the case that a zero value is incorrect, more serious problems occur. The discrepancies between a pair of technical replicates due to zero values will appear to be much larger than the discrepancies due to other (correctable) technical factors. RUV-III will attempt to correct for the larger (uncorrectable) discrepancy, and ignore the correctable technical factors. 
#' 
#' RUV-III-C is a variation of RUV-III that attempts to solve this problem, by applying RUV-III separately to every variable. If variable X is being corrected, we take the rows of the data matrix for which X is non-missing. RUV-III is then applied, and the corrected values of X is retained. The corrected values of all other variables are discarded. Note that when we take a subset of the rows of the data matrix, other columns besides X will still have missing values. These values are replaced with zero in order to apply RUV-III. No additional transformation is applied to the input data matrix. If normalization should be applied on the log-scale, then logged data must be input. 
#'
#' There are two implementations of this function, the preferred C++ version and the original protoype R code. Select which version using the \code{version} argument, which must be either "CPP" or "R"
#' 
#' @param k The number of factors of unwanted variation to remove
#' @param Y The input data matrix. Must be a matrix, not a data.frame. It should contain missing (NA) values, rather than zeros. No additional transformation is applied to the input data. 
#' @param M The design matrix containing information about technical replicates. It should not contain an intercept term!
#' @param toCorrect The names of the variables to correct using RUV-III-C
#' @param controls The names of the control variables which are known to be constant across the observations
#' @param withExtra Should we generate extra information?
#' @param withAlpha Should we generate, per-peptide, the matrix alpha giving the effects of the unwanted factors?
#' @param withW Should we generate the matrices W giving information about the unwanted factors, for every peptide?
#' @param version The version of the underlying code to use. Must be either "CPP" or "R"
#' @param progress Should a progress bar be displayed?
#' @param ... Other arguments for the prototype R code. Supported values are \code{filename} for a checkpoint file, and \code{batchSize} for the frequency with which the checkpoint file is written. 
#'
#' @return If withExtra = FALSE, returns a matrix. If withExtra = TRUE, returns a list with entries named \code{newY}, \code{residualDimensions} and \code{W}.
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
#' potentialControlsAlwaysFound <- names(which(apply(onlyPeptideData[, potentialControls], 2, 
#'     function(x) sum(is.na(x))) == 0))
#' #Actually run correction
#' #Set number of threads for CRAN
#' try(RUVIIIC::omp_set_num_threads(2L), silent=TRUE)
#' \donttest{results <- RUVIII_C(k = 11, Y = log10(onlyPeptideData), M = M, toCorrect = 
#'     colnames(onlyPeptideData), controls = potentialControlsAlwaysFound)}
#' @include RcppExports.R
#' @export
RUVIII_C <- function(k, Y, M, toCorrect, controls, withExtra = FALSE, withW = FALSE, withAlpha = FALSE, version = "CPP", progress = TRUE, ...)
{
	if(nrow(M) != nrow(Y))
	{
		stop("The number of rows in Y and M must be identical")
	}
	if(length(controls) != length(unique(controls)))
	{
		stop("Controls were not unique")
	}
	if(!all(controls %in% colnames(Y)))
	{
		stop("Not all controls were columns of Y")
	}
	if(version == "CPP")
	{
		return(RUVIIIC_CPP(k = k, Y = Y, M = M, toCorrect = toCorrect, controls = controls, withExtra = withExtra, withW = withW, withAlpha = withAlpha, progress = progress))
	}
	else if(version == "R")
	{
		return(RUVIIIC_R(k = k, Y = Y, M = M, toCorrect = toCorrect, controls = controls, withExtra = withExtra, withW = withW, withAlpha = withAlpha, ...))
	}
	else
	{
		stop("version must be either 'CPP' or 'R'")
	}
}
RUVIIIC_R <- function(k, Y, M, toCorrect, filename, controls, withExtra = FALSE, withW = FALSE, withAlpha = FALSE, batchSize = 1000)
{
	if(missing(filename))
	{
		stop("Function RUVIII_C_Varying requires a filename for intermediate results, or NULL if an intermediate file should not be used")
	}
	if(k <= 0)
	{
		stop("Input k, the number of factors of unwanted variation, must be positive")
	}
	if(is.factor(controls))
	{
		stop("Input controls was a factor. Please input a vector of variable names")
	}
	if(is.null(rownames(Y)))
	{
		stop("Function RUVIII_C_Varying requires row-names for input Y")
	}
	if(any(apply(M, 2, sum) == nrow(M)))
	{
		warning("Input design matrix had a column of all ones, indicating an intercept term. This is generally an error!")
	}
	if(any(is.na(Y[, controls])))
	{
		stop("The negative control variables should never be NA")
	}
	if(k >= nrow(Y))
	{
		stop("Input k cannot be larger than or equal to the number of rows in the input matrix")
	}
	if(k > length(controls))
	{
		stop("Input k cannot be larger than the number of negative controls")
	}
	if(nrow(M) != nrow(Y))
	{
		stop("The number of rows in Y and M must be identical")
	}
	if(length(controls) != length(unique(controls)))
	{
		stop("Controls were not unique")
	}
	if(!all(controls %in% colnames(Y)))
	{
		stop("Not all controls were columns of Y")
	}
	#Replace NAs with 0
	YWithoutNA <- Y
	YWithoutNA[is.na(YWithoutNA)] <- 0

	results <- list()
	results$peptideResults <- results$alphaResults <- results$W <- list()

	results$residualDimensions <- vector(mode = "integer", length = length(toCorrect))
	names(results$residualDimensions) <- toCorrect
	results$residualDimensions[] <- -1L

	#Load previous results set. 
	if(!is.null(filename) && file.exists(filename))
	{
		load(filename)
	}
	pb <- progress_bar$new(total = length(toCorrect), format = "[:bar] :percent :eta")

	symmetrised <- YWithoutNA %*% t(YWithoutNA)
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
				indices <- which(!is.na(Y[,peptide]))
				#If tere are no non-missing observations, then we can just mark this peptide as uncorrected / uncorrectable
				if(length(indices) > 0)
				{
					#Create the matrix that select rows of the input data matrix
					selectRows <- matrix(0, nrow = length(indices), ncol = nrow(YWithoutNA))
					selectRows[cbind(1:length(indices), indices)] <- 1
					#Submatrix of the entire data matrix
					submatrix <- YWithoutNA[indices, , drop = F]
					#Submatrix of the replicate matrix
					Msubset <- M[indices, , drop = F]
					#Some whole biological replicates are removed in the process. So remove those entire columns. This has consequences for the dimension tests later on. 
					Msubset <- Msubset[, apply(Msubset, 2, function(x) sum(x) > 0), drop=F]
					m <- nrow(submatrix)
					residualDimensions <- m - ncol(Msubset)
					results$residualDimensions[peptide] <- residualDimensions
					#You need at least two observations, across two different biological samples, in ordor to make any kind of correction
					#Also, if the dimensions don't work, we can't correct this variable
					if(ncol(Msubset) < 2 || nrow(Msubset) < 2 || min(residualDimensions, length(controls)) < k) 
					{
						#Mark this peptide as uncorrected / uncorrectable.
						newComputation <- TRUE
						results$alphaResults[peptide] <- list(c())
						results$peptideResults[[peptide]] <- rep(as.numeric(NA), nrow(Y))
						names(results$peptideResults[[peptide]]) <- rownames(YWithoutNA)
						results$W[peptide] <- list(c())
						next
					}
					#In the special case that we're trying to remove the maximum possible number of factors from this peptide, the eigen decomposition seems like it can be unstable. Fortunately, there is also the simpler GLS formulation, which may avoid those problems. 
					else if(residualDimensions == k)
					{
						#Use only the columns of Msubset that have one or more replicates to get out alphaHat
						orthogonalProjection <- diag(1, m) - Msubset %*% solve(t(Msubset) %*% Msubset) %*% t(Msubset)
						productDesign <- Msubset %*% t(Msubset)
						Uk <- svd(orthogonalProjection)$u[, 1:k]
						
						#These are used in the GLS formula
						productControls <- submatrix[, controls] %*% t(submatrix[, controls])
						invMatrix <- solve(productControls)

						adjusted <- vector(mode = "numeric", length = nrow(YWithoutNA))
						adjusted[] <- NA
						adjusted[indices] <- Msubset %*% solve(t(Msubset) %*% invMatrix %*% Msubset) %*% t(Msubset) %*% invMatrix %*% submatrix[, peptide]
						names(adjusted) <- rownames(YWithoutNA)
						results$peptideResults[[peptide]] <- adjusted
						if(withAlpha) results$alphaResults[[peptide]] <- t(Uk) %*% submatrix
						if(withW)
						{
							results$W[[peptide]] <- matrix(nrow = nrow(YWithoutNA), ncol = k)
							results$W[[peptide]][indices, ] <- productControls %*% Uk %*% solve(t(Uk) %*% productControls %*% Uk)
						}
					}
					else
					{
						orthogonalProjection <- diag(1, m) - Msubset %*% solve(t(Msubset) %*% Msubset) %*% t(Msubset)
						try({
							#Now the RUV-III code. This may throw exceptions, possibly for numerical reasons, hence the try / catch. 
							eigenDecomp <- eigs_sym(orthogonalProjection %*% selectRows %*% symmetrised %*% t(selectRows) %*% t(orthogonalProjection), k = min(m - ncol(Msubset), length(controls)), which = "LM")
							fullalpha <- t(eigenDecomp$vectors) %*% submatrix
							alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
							ac <- alpha[,controls, drop = FALSE]
							currentPeptideW <- submatrix[, controls] %*% t(ac) %*% solve(ac %*% t(ac))
							currentResult <- fullalpha[1:k, peptide]

							#Now we extend the corrected results out to a full matrix. 
							adjusted <- vector(mode = "numeric", length = nrow(YWithoutNA))
							adjusted[] <- NA
							adjusted[indices] <- submatrix[,peptide, drop=F] - currentPeptideW %*% cbind(currentResult)
							names(adjusted) <- rownames(YWithoutNA)

							results$peptideResults[[peptide]] <- adjusted
							if(withAlpha) results$alphaResults[[peptide]] <- currentResult
							if(withW)
							{
								results$W[[peptide]] <- rep(as.numeric(NA), nrow(YWithoutNA))
								results$W[[peptide]][indices] <- currentPeptideW
							}
							newComputation <- TRUE
						}, silent = TRUE)
					}
				} else
				{
					newComputation <- TRUE
					results$alphaResults[peptide] <- list(c())
					results$peptideResults[[peptide]] <- rep(as.numeric(NA), nrow(Y))
					names(results$peptideResults[[peptide]]) <- rownames(YWithoutNA)
					results$W[peptide] <- list(c())
				}
			}
			#This case only triggers if there is an exception in the try() above. 
			if(!(peptide %in% names(results$peptideResults)))
			{
				newComputation <- TRUE
				results$alphaResults[peptide] <- list(c())
				results$peptideResults[[peptide]] <- rep(as.numeric(NA), nrow(Y))
				names(results$peptideResults[[peptide]]) <- rownames(YWithoutNA)
				results$W[peptide] <- list(c())
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
	if(withExtra)
	{
		finalResult <- list(newY = do.call(cbind, results$peptideResults), residualDimensions = results$residualDimensions)
		if(withW) finalResult$W <- results$W
		if(withAlpha) finalResult$alpha <- results$alphaResults
		return(finalResult)
	}
	else
	{
		return(do.call(cbind, results$peptideResults))
	}
}
