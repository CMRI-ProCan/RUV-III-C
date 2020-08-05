test_that("Test simple scenario 1", 
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
		#Set 50 values to NA
		data[sample(1:length(data), 50)] <- NA
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 1, M = M, potentialControls = as.character(11:20), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 1, M = M, potentialControls = as.character(11:20), toCorrect = as.character(toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 2, M = M, potentialControls = as.character(11:20), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 2, M = M, potentialControls = as.character(11:20), toCorrect = as.character(toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
	  })
#Same as test 1, but reverse column order
test_that("Test simple scenario 2", 
	  {
		set.seed(2)
	  	#Random data, 20 x 20
		data <- rnorm(n = 20*20)
		dim(data) <- c(20, 20)

		colnames(data) <- 1:20
		rownames(data) <- 1:20

		M <- data.frame(rep = factor(rep(1:10, each = 2)))
		M <- model.matrix(~ rep - 1, data = M)
		M <- data.matrix(M)
		#Set 50 values to NA
		data[sample(1:length(data), 50)] <- NA
		#Reverse column order, so NAs now at the end of the matrix. But still in columns with the same names
		data <- data[, 20:1]
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 1, M = M, potentialControls = as.character(11:20), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 1, M = M, potentialControls = as.character(11:20), toCorrect = as.character(toCorrect:(toCorrect + 1)), filename = NULL, withExtra = FALSE, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 2, M = M, potentialControls = as.character(11:20), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 2, M = M, potentialControls = as.character(11:20), toCorrect = as.character(toCorrect:(toCorrect + 1)), filename = NULL, withExtra = FALSE, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
	  })
#Same as test 2, but with the column names still 1 -> 20
test_that("Test simple scenario 3", 
	  {
		set.seed(3)
	  	#Random data, 20 x 20
		data <- rnorm(n = 20*20)
		dim(data) <- c(20, 20)

		M <- data.frame(rep = factor(rep(1:10, each = 2)))
		M <- model.matrix(~ rep - 1, data = M)
		M <- data.matrix(M)
		#Set 50 values to NA
		data[sample(1:length(data), 50)] <- NA
		#Reverse column order. 
		data <- data[, 20:1]

		#NAs are now in columns named 6:20
		colnames(data) <- 1:20
		rownames(data) <- 1:20
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 1, M = M, potentialControls = as.character(1:9), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 1, M = M, potentialControls = as.character(1:9), toCorrect = as.character(toCorrect:(toCorrect + 1)), filename = NULL, withExtra = FALSE, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 2, M = M, potentialControls = as.character(1:9), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 2, M = M, potentialControls = as.character(1:9), toCorrect = as.character(toCorrect:(toCorrect + 1)), filename = NULL, withExtra = FALSE, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
	  })
test_that("Test that there was nothing special about numeric column names", 
	  {
		set.seed(3)
	  	#Random data, 20 x 20
		data <- rnorm(n = 20*20)
		dim(data) <- c(20, 20)

		M <- data.frame(rep = factor(rep(1:10, each = 2)))
		M <- model.matrix(~ rep - 1, data = M)
		M <- data.matrix(M)
		#Set 50 values to NA
		data[sample(1:length(data), 50)] <- NA
		#Reverse column order. 
		data <- data[, 20:1]

		#NAs are now in columns named 6:20
		colnames(data) <- paste0("C", 1:20)
		rownames(data) <- paste0("C", 1:20)
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 1, M = M, potentialControls = paste0("C", 1:9), toCorrect = paste0("C", toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 1, M = M, potentialControls = paste0("C", 1:9), toCorrect = paste0("C", toCorrect:(toCorrect + 1)), filename = NULL, withExtra = FALSE, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 2, M = M, potentialControls = paste0("C", 1:9), toCorrect = paste0("C", toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 2, M = M, potentialControls = paste0("C", 1:9), toCorrect = paste0("C", toCorrect:(toCorrect + 1)), filename = NULL, withExtra = FALSE, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
	  })
test_that("Test when some peptides can be corrected, and some can't", 
	  {
		set.seed(3)
	  	#Random data, 20 x 20
		data <- rnorm(n = 20*20)
		dim(data) <- c(20, 20)

		M <- data.frame(rep = factor(rep(1:10, each = 2)))
		M <- model.matrix(~ rep - 1, data = M)
		M <- data.matrix(M)
	
		#Almost all the values for the first 4 peptides are missing
		data[1:17, 1:4] <- NA
		data[sample(1:length(data), 50)] <- NA

		#NAs are now in columns named 6:20
		colnames(data) <- paste0("C", 1:20)
		rownames(data) <- paste0("C", 1:20)
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_Varying_CPP(Y = data, k = 5, M = M, potentialControls = paste0("C", 11:20), toCorrect = paste0("C", toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_Varying_R(Y = data, k = 5, M = M, potentialControls = paste0("C", 11:20), toCorrect = paste0("C", toCorrect:(toCorrect + 1)), filename = NULL, withExtra = FALSE, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
	  })
