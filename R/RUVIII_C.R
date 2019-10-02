RUVIII_C <- function(k, ruvInputData, M, toCorrect, filename, controls, withW = FALSE, batchSize = 1000)
{
	if(missing(filename))
	{
		stop("Function RUVIII_NM_Varying requires a filename for intermediate results")
	}
	if(k <= 0)
	{
		stop("Input k, the number of factors of unwanted variation, must be positive")
	}
	if(is.null(rownames(ruvInputData)))
	{
		stop("Function RUVIII_NM_Varying requires row-names for input ruvInputData")
	}
	#Replace NAs with 0
	ruvInputDataWithoutNA <- ruvInputData
	ruvInputDataWithoutNA[is.na(ruvInputDataWithoutNA)] <- 0

	results <- list()
	results$peptideResults <- results$alphaResults <- results$W <- list()

	#Load previous results set. 
	if(file.exists(filename))
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
			if(peptideIndex %% batchSize == 0 && newComputation)
			{
				save(results, file = filename)
			}
			pb$tick()
		}
		#Exit the progress bar.
		pb$terminate()
		#If we've made any new computations write out to file. 
		if(newComputation) save(results, file = filename)
	}
	return(results)
}
