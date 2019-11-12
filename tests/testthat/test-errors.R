test_that("Test NAs in negative control columns",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20
		data[1,1] <- NA

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)

		expect_error(RUVIIIC:::RUVIII_C_CPP(input = data, k = 1, M = M, toCorrect = as.character(15:20), controls = as.character(1), withW = FALSE), "The negative control variables should never be NA")
                expect_error(RUVIIIC:::RUVIII_C_R(ruvInputData = data, k = 1, M = M, toCorrect = as.character(15:20), controls = as.character(1), withW = FALSE, filename = NULL), "The negative control variables should never be NA")
	})

test_that("Test that there's an error when k is too large",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20
		data[1,1] <- NA

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)

		#Larger than possible, considering the number of rows in the initial data matrix
		expect_error(RUVIIIC:::RUVIII_C_CPP(input = data, k = 30, M = M, controls = as.character(15:20), toCorrect = as.character(1), withW = FALSE), "Input k cannot be larger than or equal to the number of rows in the input matrix")
                expect_error(RUVIIIC:::RUVIII_C_R(ruvInputData = data, k = 30, M = M, controls = as.character(15:20), toCorrect = as.character(1), withW = FALSE, filename = NULL), "Input k cannot be larger than or equal to the number of rows in the input matrix")
	})

test_that("Test that there's an error when k is too large",
	{
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20
		data[1,1] <- NA

                M <- data.frame(rep = factor(rep(1:10, each = 2)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)

		#Larger than possible, considering the number of controls
		expect_error(RUVIIIC:::RUVIII_C_CPP(input = data, k = 10, M = M, controls = as.character(15:20), toCorrect = as.character(1), withW = FALSE), "Input k cannot be larger than the number of negative controls")
                expect_error(RUVIIIC:::RUVIII_C_R(ruvInputData = data, k = 10, M = M, controls = as.character(15:20), toCorrect = as.character(1), withW = FALSE, filename = NULL), "Input k cannot be larger than the number of negative controls")
	})

test_that("Test that there's no bug in input to Spectra::SymEigsSolver", 
	  {
                set.seed(1)
                #Random data, 20 x 20
                data <- rnorm(n = 20*20)
                dim(data) <- c(20, 20)

                colnames(data) <- 1:20
                rownames(data) <- 1:20

                M <- data.frame(rep = factor(rep(1:2, each = 10)))
                M <- model.matrix(~ rep - 1, data = M)
                M <- data.matrix(M)

                copiedData <- data
                naIndex <- 1
                targetColumn <- 1
                copiedData[naIndex, targetColumn] <- NA
                result <- RUVIII_C(ruvInputData = copiedData, k = 3, M = M, toCorrect = as.character(targetColumn), controls = as.character(16:20), withW = FALSE, filename = NULL, version = "CPP")
                result <- RUVIII_C(ruvInputData = copiedData, k = 3, M = M, toCorrect = as.character(targetColumn), controls = as.character(16:20), withW = FALSE, filename = NULL, version = "R")
                result <- RUVIII_C_Varying(ruvInputData = copiedData, k = 3, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withW = FALSE, filename = NULL, version = "CPP")
		result <- RUVIII_C_Varying(ruvInputData = copiedData, k = 3, M = M, toCorrect = as.character(targetColumn), potentialControls = as.character(16:20), withW = FALSE, filename = NULL, version = "R")
		expect_identical(TRUE, TRUE)
	  })
