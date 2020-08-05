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
					result <- RUVIII_C(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), controls = as.character(18:20), withExtra = FALSE, filename = NULL, version = "CPP", progress = FALSE)
					result_standardRUVIII <- ruv::RUVIII(Y = copiedData[-naIndex, ], M = M[-naIndex, ], k = k, ctl = 18:20)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
					
					result <- RUVIII_C(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), controls = as.character(18:20), withExtra = FALSE, filename = NULL, version = "R", progress = FALSE)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
				}
			}
		}
	})

test_that("Test NAs in other columns too",
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
					result <- RUVIII_C(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), controls = as.character(18:20), withExtra = FALSE, filename = NULL, version = "CPP", progress = FALSE)

					RUVIIIdata <- copiedData
					RUVIIIdata[naIndex+1, targetColumn+1] <- 0
					result_standardRUVIII <- ruv::RUVIII(Y = RUVIIIdata[-naIndex, ], M = M[-naIndex, ], k = k, ctl = 18:20)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
					
					result <- RUVIII_C(Y = copiedData, k = k, M = M, toCorrect = as.character(targetColumn), controls = as.character(18:20), withExtra = FALSE, filename = NULL, version = "R", progress = FALSE)
					expect_equal(result[-naIndex, , drop=F], result_standardRUVIII[, targetColumn, drop=F])
				}
			}
		}
	})
