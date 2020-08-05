test_that("Test with 50 controls, R", 
	{
		for(counter in 1:4)
		{
			simulatedData <- rnorm(231 * 1000)
			dim(simulatedData) <- c(231, 1000)
			colnames(simulatedData) <- paste0("Var", 1:1000)
			rownames(simulatedData) <- paste0("Obs", 1:231)
			simulatedData[, "Var1"] <- NA

			#Non NA indices
			indices <- c(25, 36, 37, 42, 58, 60, 66, 69, 75, 135, 137, 146, 147, 154, 155, 162, 186, 187, 188)
			simulatedData[indices, "Var1"] <- rnorm(length(indices))

			factor <- as.character(1:231)
			factor[c(36, 37)] <- "Pair1"
			factor[c(146, 147)] <- "Pair2"
			factor[c(154, 155)] <- "Pair3"
			factor[c(186, 187)] <- "Pair4"
			factor <- as.factor(factor)
			M <- model.matrix( ~ M - 1, data = data.frame(M = factor))
			result <- RUVIII_C(k = 4, Y = simulatedData, M = M, toCorrect = "Var1", filename = NULL, controls = paste0("Var", 950:1000), withExtra = TRUE, withW = TRUE, withAlpha = TRUE, version = "R", progress = FALSE)
			#Check that we made the correction that we expected to, here
			expect_equal(result$newY[indices, "Var1"], simulatedData[indices, "Var1"] - (result$W[["Var1"]][indices, ] %*% result$alpha[["Var1"]])[, "Var1"])

			vanillaM <- M[indices, ]
			vanillaM <- vanillaM[, apply(vanillaM, 2, sum) > 0]
			resultVanilla <- ruv::RUVIII(k = 4, Y = simulatedData[indices, ], M = vanillaM, ctl = 950:1000)
			expect_equal(resultVanilla[,1], result$newY[indices, "Var1"])
	
	
			result_varying <- RUVIII_C(k = 4, Y = simulatedData, M = M, toCorrect = "Var1", filename = NULL, controls = paste0("Var", 950:1000), withExtra = TRUE, withW = TRUE, withAlpha = TRUE, version = "R", progress = FALSE)
			expect_equal(result_varying$newY[indices, "Var1"], simulatedData[indices, "Var1"] - (result_varying$W[["Var1"]][indices, ] %*% result_varying$alpha[["Var1"]])[, "Var1"])
			expect_equal(result$newY[indices, "Var1"], result_varying$newY[indices, "Var1"])
		}
	})
test_that("Test with 50 controls, CPP", 
	{
		for(counter in 1:4)
		{
			simulatedData <- rnorm(231 * 1000)
			dim(simulatedData) <- c(231, 1000)
			colnames(simulatedData) <- paste0("Var", 1:1000)
			rownames(simulatedData) <- paste0("Obs", 1:231)
			simulatedData[, "Var1"] <- NA

			#Non NA indices
			indices <- c(25, 36, 37, 42, 58, 60, 66, 69, 75, 135, 137, 146, 147, 154, 155, 162, 186, 187, 188)
			simulatedData[indices, "Var1"] <- rnorm(length(indices))

			factor <- as.character(1:231)
			factor[c(36, 37)] <- "Pair1"
			factor[c(146, 147)] <- "Pair2"
			factor[c(154, 155)] <- "Pair3"
			factor[c(186, 187)] <- "Pair4"
			factor <- as.factor(factor)
			M <- model.matrix( ~ M - 1, data = data.frame(M = factor))
			result <- RUVIII_C(k = 4, Y = simulatedData, M = M, toCorrect = "Var1", filename = NULL, controls = paste0("Var", 950:1000), withExtra = TRUE, withW = TRUE, withAlpha = TRUE, version = "CPP", progress = FALSE)
			#Check that we made the correction that we expected to, here
			expect_equal(result$newY[indices, "Var1"], simulatedData[indices, "Var1"] - (result$W[["Var1"]][indices, ] %*% result$alpha[["Var1"]])[, "Var1"])

			vanillaM <- M[indices, ]
			vanillaM <- vanillaM[, apply(vanillaM, 2, sum) > 0]
			resultVanilla <- ruv::RUVIII(k = 4, Y = simulatedData[indices, ], M = vanillaM, ctl = 950:1000)
			expect_equal(resultVanilla[,1], result$newY[indices, "Var1"])
	
	
			result_varying <- RUVIII_C(k = 4, Y = simulatedData, M = M, toCorrect = "Var1", filename = NULL, controls = paste0("Var", 950:1000), withExtra = TRUE, withW = TRUE, withAlpha = TRUE, version = "CPP", progress = FALSE)
			expect_equal(result_varying$newY[indices, "Var1"], simulatedData[indices, "Var1"] - (result_varying$W[["Var1"]][indices, ] %*% result_varying$alpha[["Var1"]])[, "Var1"])
			expect_equal(result$newY[indices, "Var1"], result_varying$newY[indices, "Var1"])
		}
	})
test_that("Test with 25 controls, R", 
	{
		for(counter in 1:4)
		{
			simulatedData <- rnorm(231 * 1000)
			dim(simulatedData) <- c(231, 1000)
			colnames(simulatedData) <- paste0("Var", 1:1000)
			rownames(simulatedData) <- paste0("Obs", 1:231)
			simulatedData[, "Var1"] <- NA

			#Non NA indices
			indices <- c(25, 36, 37, 42, 58, 60, 66, 69, 75, 135, 137, 146, 147, 154, 155, 162, 186, 187, 188)
			simulatedData[indices, "Var1"] <- rnorm(length(indices))
			simulatedData[25, 950:974] <- NA

			factor <- as.character(1:231)
			factor[c(36, 37)] <- "Pair1"
			factor[c(146, 147)] <- "Pair2"
			factor[c(154, 155)] <- "Pair3"
			factor[c(186, 187)] <- "Pair4"
			factor <- as.factor(factor)
			M <- model.matrix( ~ M - 1, data = data.frame(M = factor))
			result_varying <- RUVIII_C_Varying(k = 4, Y = simulatedData, M = M, toCorrect = "Var1", filename = NULL, potentialControls = paste0("Var", 950:1000), withExtra = TRUE, withW = TRUE, withAlpha = TRUE, version = "R", progress = FALSE)
			#Check that we made the correction that we expected to, here
			expect_equal(result_varying$newY[indices, "Var1"], simulatedData[indices, "Var1"] - (result_varying$W[["Var1"]][indices, ] %*% result_varying$alpha[["Var1"]])[, "Var1"])

			vanillaM <- M[indices, ]
			vanillaM <- vanillaM[, apply(vanillaM, 2, sum) > 0]
			vanillaData <- simulatedData[indices, ]
			vanillaData[is.na(vanillaData)] <- 0
			resultVanilla <- ruv::RUVIII(k = 4, Y = vanillaData, M = vanillaM, ctl = 975:1000)
			expect_equal(resultVanilla[,1], result_varying$newY[indices, "Var1"])
		}
	})
test_that("Test with 25 controls, CPP", 
	{
		for(counter in 1:4)
		{
			simulatedData <- rnorm(231 * 1000)
			dim(simulatedData) <- c(231, 1000)
			colnames(simulatedData) <- paste0("Var", 1:1000)
			rownames(simulatedData) <- paste0("Obs", 1:231)
			simulatedData[, "Var1"] <- NA

			#Non NA indices
			indices <- c(25, 36, 37, 42, 58, 60, 66, 69, 75, 135, 137, 146, 147, 154, 155, 162, 186, 187, 188)
			simulatedData[indices, "Var1"] <- rnorm(length(indices))
			simulatedData[25, 950:974] <- NA

			factor <- as.character(1:231)
			factor[c(36, 37)] <- "Pair1"
			factor[c(146, 147)] <- "Pair2"
			factor[c(154, 155)] <- "Pair3"
			factor[c(186, 187)] <- "Pair4"
			factor <- as.factor(factor)
			M <- model.matrix( ~ M - 1, data = data.frame(M = factor))
			result_varying <- RUVIII_C_Varying(k = 4, Y = simulatedData, M = M, toCorrect = "Var1", filename = NULL, potentialControls = paste0("Var", 950:1000), withExtra = TRUE, withW = TRUE, withAlpha = TRUE, version = "CPP", progress = FALSE)
			#Check that we made the correction that we expected to, here
			expect_equal(result_varying$newY[indices, "Var1"], simulatedData[indices, "Var1"] - (result_varying$W[["Var1"]][indices, ] %*% result_varying$alpha[["Var1"]])[, "Var1"])

			vanillaM <- M[indices, ]
			vanillaM <- vanillaM[, apply(vanillaM, 2, sum) > 0]
			vanillaData <- simulatedData[indices, ]
			vanillaData[is.na(vanillaData)] <- 0
			resultVanilla <- ruv::RUVIII(k = 4, Y = vanillaData, M = vanillaM, ctl = 975:1000)
			expect_equal(resultVanilla[,1], result_varying$newY[indices, "Var1"])
		}
	})
