test_that("Test residual dimensions", 
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
		for(k in 2:3)
		{
			for(version in c("R", "CPP"))
			{
				copiedData <- data
				copiedData[1, 1] <- NA
				copiedData[1:2, 2] <- NA
				copiedData[1:3, 3] <- NA
				copiedData[1:19, 4] <- NA
				result <- RUVIII_C(Y = copiedData, k = k, M = M, toCorrect = as.character(1:5), controls = as.character(15:20), withExtra = TRUE, filename = NULL, version = version, progress = FALSE)
				expect_identical(unname(result$residualDimensions["1"]), 9L)
				expect_identical(unname(result$residualDimensions["2"]), 9L)
				expect_identical(unname(result$residualDimensions["3"]), 8L)
				expect_identical(unname(result$residualDimensions["4"]), 0L)
				expect_identical(unname(result$residualDimensions["5"]), 10L)

				result <- RUVIII_C_Varying(Y = copiedData, k = k, M = M, toCorrect = as.character(1:5), potentialControls = as.character(15:20), withExtra = TRUE, filename = NULL, version = version, progress = FALSE)
				expect_identical(unname(result$residualDimensions["1"]), 9L)
				expect_identical(unname(result$residualDimensions["2"]), 9L)
				expect_identical(unname(result$residualDimensions["3"]), 8L)
				expect_identical(unname(result$residualDimensions["4"]), 0L)
				expect_identical(unname(result$residualDimensions["5"]), 10L)
			}
		}
	})
