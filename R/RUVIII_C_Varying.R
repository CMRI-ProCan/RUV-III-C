#' RUV-III-C, varying controls
#'
#' Apply RUV-III-C, a variation of RUV-III that only uses non-missing values
#' 
#' See the documentation of \link{RUVIII_C} for more information about the RUV-III-C method. This function is identical, except in this case the set of negative control variables actually used varies depending on the target variable to be normalized. Instead of putting in a list of negative control variables, the user specifies a list of potential negatve control variables. 
#'
#' When normalizing variable X, the algorithm begins by selecting the rows of the data matrix for which X is non-missing. Out of the potential negative control peptides, it selects those that are always non-missing across the selected subset. The standard version of RUV-III is then applied, similar to \link{RUVIII_C}. 
#' 
#' There are two implementations of this function, one in C++ and one in R. Select which version using the \code{version} argument, which must be either "CPP" or "R"
#' 
#' @param k The number of factors of unwanted variation to remove
#' @param ruvInputData The input data matrix. Must be a matrix, not a data.frame. It should contain missing (NA) values, rather than zeros. 
#' @param M The design matrix containing information about technical replicates. It should not contain an intercept term!
#' @param toCorrect The names of the variables to correct using RUV-III-C
#' @param filename The intermediate file in which to save the results. 
#' @param potentialControls The names of the control variables which are known to be constant across the observations
#' @param withW Should we keep the estimate of the unwanted variation factors W, for each variable which is corrected?
#' @param batchSize How often should we write to the intermediate file? The default of 1000 implies that results are written to file every 1000 variables. 
#' @param version The version of the underlying code to use. Must be either "CPP" or "R"
#'
#' @return If withW = FALSE, returns a matrix. If withW = TRUE, returns a list with entries named \code{newY} and \code{W}.
#'
#' @export
RUVIII_C_Varying <- function(k, ruvInputData, M, toCorrect, filename, potentialControls, withW = FALSE, batchSize = 1000, version = "CPP")
{
	if(version == "CPP")
	{
		return(RUVIII_C_Varying_CPP(k = k, input = ruvInputData, M = M, toCorrect = toCorrect, potentialControls = potentialControls, withW = withW))
	}
	else if(version == "R")
	{
		return(RUVIII_C_Varying_R(k = k, ruvInputData = ruvInputData, M = M, toCorrect = toCorrect, filename = filename, potentialControls = potentialControls, withW = withW, batchSize = batchSize))
	}
	else
	{
		stop("version must be either 'CPP' or 'R'")
	}
}
RUVIII_C_Varying_R <- function(k, ruvInputData, M, toCorrect, filename, potentialControls, withW = FALSE, batchSize = 1000)
{
	if(missing(filename))
	{
		stop("Function RUVIII_C_Varying requires a filename for intermediate results, or NULL if an intermediate file should not be used")
	}
	if(k <= 0)
	{
		stop("Input k, the number of factors of unwanted variation, must be positive")
	}
	if(is.null(rownames(ruvInputData)))
	{
		stop("Function RUVIII_NM_Varying requires row-names for input ruvInputData")
	}
	if(k >= nrow(ruvInputData))
	{
		stop("Input k cannot be larger than or equal to the number of rows in the input matrix")
	}
	if(k > length(potentialControls))
	{
		stop("Input k cannot be larger than the number of negative controls")
	}
	#Replace NAs with 0
	ruvInputDataWithoutNA <- ruvInputData
	ruvInputDataWithoutNA[is.na(ruvInputDataWithoutNA)] <- 0

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

	symmetrised <- ruvInputDataWithoutNA %*% t(ruvInputDataWithoutNA)
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
				#Find the controls which are non-missing, for those rows. 
				controlsThisPeptide <- names(which(apply(ruvInputData[indices, potentialControls, drop = F], 2, function(x) sum(is.na(x)) == 0)))
				#If there are no non-missing observations, then we can just mark this peptide as uncorrected / uncorrectable
				if(length(indices) > 0)
				{
					#Create the matrix that select rows of the input data matrix
					selectRows <- matrix(0, nrow = length(indices), ncol = nrow(ruvInputDataWithoutNA))
					selectRows[cbind(1:length(indices), indices)] <- 1
					#Submatrix of the entire data matrix
					submatrix <- ruvInputDataWithoutNA[indices, , drop = F]
					#Submatrix of the replicate matrix
					Msubset <- M[indices, , drop = F]
					#Some whole biological replicates are removed in the process. So remove those entire columns. This has consequences for the dimension tests later on. 
					Msubset <- Msubset[, apply(Msubset, 2, function(x) sum(x) > 0), drop=F]
					#Number of observations
					m <- nrow(submatrix)
					residualDimensions <- m - ncol(Msubset)
					results$residualDimensions[peptide] <- residualDimensions
					#You need at least two observations, across two different biological samples, in ordor to make any kind of correction
					if(ncol(Msubset) < 2 || nrow(Msubset) < 2) 
					{
						#Mark this peptide as uncorrected / uncorrectable.
						newComputation <- TRUE
						results$alphaResults[peptide] <- list(c())
						results$peptideResults[[peptide]] <- rep(as.numeric(NA), nrow(ruvInputData))
						names(results$peptideResults[[peptide]]) <- rownames(ruvInputDataWithoutNA)
						if(withW) results$W[peptide] <- list(c())
						next
					}
					#If the dimensions don't work, we can't correct this variable
					if(min(residualDimensions, length(controlsThisPeptide)) < k)
					{
						newComputation <- TRUE
						results$alphaResults[peptide] <- list(c())
						results$peptideResults[[peptide]] <- rep(as.numeric(NA), nrow(ruvInputData))
						names(results$peptideResults[[peptide]]) <- rownames(ruvInputDataWithoutNA)
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
						productControls <- submatrix[, controlsThisPeptide] %*% t(submatrix[, controlsThisPeptide])
						invMatrix <- solve(productControls)

						adjusted <- vector(mode = "numeric", length = nrow(ruvInputDataWithoutNA))
						adjusted[] <- NA
						adjusted[indices] <- Msubset %*% solve(t(Msubset) %*% invMatrix %*% Msubset) %*% t(Msubset) %*% invMatrix %*% submatrix[, peptide]
						names(adjusted) <- rownames(ruvInputDataWithoutNA)
						results$peptideResults[[peptide]] <- adjusted
						results$alphaResults[[peptide]] <- t(Uk) %*% submatrix
						if(withW)
						{
							results$W[[peptide]] <- matrix(nrow = nrow(ruvInputDataWithoutNA), ncol = k)
							results$W[[peptide]][indices, ] <- productControls %*% Uk %*% solve(t(Uk) %*% productControls %*% Uk)
						}
					}
					else
					{
						orthogonalProjection <- diag(1, m) - Msubset %*% solve(t(Msubset) %*% Msubset) %*% t(Msubset)
						try({
							#Now the RUV-III code. This may throw exceptions, possibly for numerical reasons, hence the try / catch. 
							eigenDecomp <- eigs_sym(orthogonalProjection %*% selectRows %*% symmetrised %*% t(selectRows) %*% t(orthogonalProjection), k = min(m - ncol(Msubset), length(controlsThisPeptide)), which = "LM")
							fullalpha <- t(eigenDecomp$vectors) %*% submatrix
							alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
							ac <- alpha[,controlsThisPeptide, drop = FALSE]
							currentPeptideW <- submatrix[, controlsThisPeptide] %*% t(ac) %*% solve(ac %*% t(ac))
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
								results$W[[peptide]] <- rep(as.numeric(NA), nrow(ruvInputDataWithoutNA))
								results$W[[peptide]][indices] <- currentPeptideW
							}
							newComputation <- TRUE
						}, silent = TRUE)
					}
				} else
				{
					newComputation <- TRUE
					results$alphaResults[peptide] <- list(c())
					results$peptideResults[[peptide]] <- rep(as.numeric(NA), nrow(ruvInputData))
					names(results$peptideResults[[peptide]]) <- rownames(ruvInputDataWithoutNA)
					results$W[peptide] <- list(c())
				}
			}
			#This case only triggers if there is an exception in the try() above. 
			if(!(peptide %in% names(results$peptideResults)))
			{
				newComputation <- TRUE
				results$alphaResults[peptide] <- list(c())
				results$peptideResults[[peptide]] <- rep(as.numeric(NA), nrow(ruvInputData))
				names(results$peptideResults[[peptide]]) <- rownames(ruvInputDataWithoutNA)
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
	if(withW)
	{
		return(list(newY = do.call(cbind, results$peptideResults), W = results$W, alpha = results$alpha, residualDimensions = results$residualDimensions))
	}
	else
	{
		return(do.call(cbind, results$peptideResults))
	}
}
