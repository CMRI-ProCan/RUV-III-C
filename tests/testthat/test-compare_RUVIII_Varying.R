test_that("Test NAs in only the corrected column",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)
		for(targetColumn in 1:2)
		{
			for(k in 1:2)
			{
				for(naIndex in seq(1, 10, by = 2))
				{
					copiedData <- data
					copiedData[naIndex, targetColumn] <- NA
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(18:20), withExtra = FALSE, filename = NULL, version = "CPP", progress = FALSE)
					result_standardRUVIII <- ruv::RUVIII(Y = copiedData[-naIndex, ], M = M[-naIndex, ], k = k, ctl = 18:20)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
					
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(18:20), withExtra = FALSE, filename = NULL, version = "R", progress = FALSE)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
				}
			}
		}
	})

test_that("Test NAs in the corrected column, and also an NA in a different row in a different non-control column",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)
		for(targetColumn in 1:2)
		{
			for(k in 1:2)
			{
				for(naIndex in seq(1, 10, by = 2))
				{
					copiedData <- data
					copiedData[naIndex, targetColumn] <- NA
					copiedData[naIndex+1, targetColumn+1] <- NA
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(18:20), withExtra = FALSE, filename = NULL, version = "CPP", progress = FALSE)

					RUVIIIdata <- copiedData
					RUVIIIdata[naIndex+1, targetColumn+1] <- 0
					result_standardRUVIII <- ruv::RUVIII(Y = RUVIIIdata[-naIndex, ], M = M[-naIndex, ], k = k, ctl = 18:20)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
					
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(18:20), withExtra = FALSE, filename = NULL, version = "R", progress = FALSE)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
				}
			}
		}
	})

test_that("Test NAs in the corrected column, and also in a control for the same observation",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)
		for(targetColumn in 1:2)
		{
			for(k in 1:2)
			{
				for(naIndex in seq(1, 10, by = 2))
				{
					copiedData <- data
					copiedData[naIndex, targetColumn] <- NA
					#Because the NA is in a control on the same line as the NA in the variable to correct, we actually use all the controls here
					copiedData[naIndex, "16"] <- NA
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withExtra = FALSE, filename = NULL, version = "CPP", progress = FALSE)
					result_standardRUVIII <- ruv::RUVIII(Y = copiedData[-naIndex, ], M = M[-naIndex, ], k = k, ctl = 16:20)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
					
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withExtra = FALSE, filename = NULL, version = "R", progress = FALSE)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
				}
			}
		}
	})
test_that("Test NAs in the corrected column, and also in a control for a different observation",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)
		for(targetColumn in 1:2)
		{
			for(k in 1:2)
			{
				for(naIndex in seq(1, 10, by = 2))
				{
					copiedData <- data
					copiedData[naIndex, targetColumn] <- NA
					#Because the NA is in a control on a different observation, this shrinks the set of controls used.
					copiedData[naIndex+1, "16"] <- NA
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withExtra = FALSE, filename = NULL, version = "CPP", progress = FALSE)
					
					RUVIIIdata <- copiedData
					RUVIIIdata[naIndex+1, "16"] <- 0
					result_standardRUVIII <- ruv::RUVIII(Y = RUVIIIdata[-naIndex, ], M = M[-naIndex, ], k = k, ctl = 17:20)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
					
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withExtra = FALSE, filename = NULL, version = "R", progress = FALSE)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
				}
			}
		}
	})

test_that("Test NAs in the corrected column, and also in a control for a different observation, and also in a variable that is not corrected and not a control",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)
		for(targetColumn in 1:2)
		{
			for(k in 1:2)
			{
				for(naIndex in seq(1, 10, by = 2))
				{
					copiedData <- data
					copiedData[naIndex, targetColumn] <- NA
					copiedData[naIndex+1, targetColumn+1] <- NA
					#Because the NA is in a control on the same line as the NA in the variable to correct, we actually use all the controls here
					copiedData[naIndex, "16"] <- NA
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withExtra = FALSE, filename = NULL, version = "CPP", progress = FALSE)

					RUVIIIdata <- copiedData
					RUVIIIdata[naIndex+1, targetColumn+1] <- 0
					result_standardRUVIII <- ruv::RUVIII(Y = RUVIIIdata[-naIndex, ], M = M[-naIndex, ], k = k, ctl = 16:20)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
					
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withExtra = FALSE, filename = NULL, version = "R", progress = FALSE)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
				}
			}
		}
	})
test_that("Test NAs in the corrected column, other noncontrol columns, and also in the controls",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)
		for(targetColumn in 1:2)
		{
			for(k in 1:2)
			{
				for(naIndex in seq(1, 10, by = 2))
				{
					copiedData <- data
					copiedData[naIndex, targetColumn] <- NA
					copiedData[naIndex+1, targetColumn+1] <- NA
					#Because the NA is in a control on the same line as the NA in the variable to correct, we actually use all the controls here
					copiedData[naIndex+1, "16"] <- NA
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withExtra = FALSE, filename = NULL, version = "CPP", progress = FALSE)

					RUVIIIdata <- copiedData
					RUVIIIdata[naIndex+1, targetColumn+1] <- 0
					RUVIIIdata[naIndex+1, "16"] <- 0
					result_standardRUVIII <- ruv::RUVIII(Y = RUVIIIdata[-naIndex, ], M = M[-naIndex, ], k = k, ctl = 17:20)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
					
					result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withExtra = FALSE, filename = NULL, version = "R", progress = FALSE)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
				}
			}
		}
	})
