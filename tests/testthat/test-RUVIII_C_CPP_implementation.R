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
		#Set 30 values in the first 15 columns to NA
		data[sample(1:(15*20), 30)] <- NA
		for(toCorrect in 1:14)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 1, M = M, controls = as.character(16:20), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 1, M = M, controls = as.character(16:20), toCorrect = as.character(toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
		for(toCorrect in 1:14)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 2, M = M, controls = as.character(16:20), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 2, M = M, controls = as.character(16:20), toCorrect = as.character(toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
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
		#Set 30 values in the first 15 columns to NA
		data[sample(1:(15*20), 30)] <- NA
		#Reverse column order, so NAs now at the end of the matrix. But still in columns with the same names
		data <- data[, 20:1]
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 1, M = M, controls = as.character(16:20), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 1, M = M, controls = as.character(16:20), toCorrect = as.character(toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 2, M = M, controls = as.character(16:20), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 2, M = M, controls = as.character(16:20), toCorrect = as.character(toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
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
		#Set 30 values in the first 15 columns to NA
		data[sample(1:(15*20), 30)] <- NA
		#Reverse column order. 
		data <- data[, 20:1]

		#NAs are now in columns named 6:20
		colnames(data) <- 1:20
		rownames(data) <- 1:20
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 1, M = M, controls = as.character(1:5), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 1, M = M, controls = as.character(1:5), toCorrect = as.character(toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 2, M = M, controls = as.character(1:5), toCorrect = as.character(toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 2, M = M, controls = as.character(1:5), toCorrect = as.character(toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
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
		#Set 30 values in the first 15 columns to NA
		data[sample(1:(15*20), 30)] <- NA
		#Reverse column order. 
		data <- data[, 20:1]

		#NAs are now in columns named 6:20
		colnames(data) <- paste0("C", 1:20)
		rownames(data) <- paste0("C", 1:20)
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 1, M = M, controls = paste0("C", 1:5), toCorrect = paste0("C", toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 1, M = M, controls = paste0("C", 1:5), toCorrect = paste0("C", toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 2, M = M, controls = paste0("C", 1:5), toCorrect = paste0("C", toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 2, M = M, controls = paste0("C", 1:5), toCorrect = paste0("C", toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
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
		data[1:18, 1:4] <- NA

		#NAs are now in columns named 6:20
		colnames(data) <- paste0("C", 1:20)
		rownames(data) <- paste0("C", 1:20)
		for(toCorrect in 1:19)
		{
			cppImplementation <- RUVIIIC:::RUVIIIC_CPP(Y = data, k = 3, M = M, controls = paste0("C", 17:20), toCorrect = paste0("C", toCorrect:(toCorrect+1)), withExtra = FALSE, withW = FALSE, withAlpha = FALSE, progress = FALSE)
			rImplementation <- RUVIIIC:::RUVIIIC_R(Y = data, k = 3, M = M, controls = paste0("C", 17:20), toCorrect = paste0("C", toCorrect:(toCorrect + 1)), withExtra = FALSE, filename = NULL, withW = FALSE, withAlpha = FALSE)
			expect_equal(cppImplementation, rImplementation)
			expect_identical(is.na(cppImplementation), is.na(rImplementation))
			if(!all(is.na(cppImplementation))) expect_lt(max(abs(cppImplementation - rImplementation), na.rm=TRUE), 1e-8)
		}
	  })
