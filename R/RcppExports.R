RUVIII_C_CPP <- function(Y , k, M, controls, toCorrect, withExtra, withW = FALSE, withAlpha = FALSE) {
        if(any(is.na(Y[, controls])))
        {
                stop("The negative control variables should never be NA")
        }
	if(k >= nrow(Y))
	{
		stop("Input k cannot be larger than or equal to the number of rows in the Y matrix")
	}
	if(k > length(controls))
	{
		stop("Input k cannot be larger than the number of negative controls")
	}
    	.Call('_RUVIIIC_RUVIIIC', PACKAGE = 'RUVIIIC', Y, k, M, controls, toCorrect, withExtra, withW, withAlpha)
}
RUVIII_C_Varying_CPP <- function(Y, k, M, potentialControls, toCorrect, withExtra, withW = FALSE, withAlpha = FALSE) {
	if(k >= nrow(Y))
	{
		stop("Input k cannot be larger than or equal to the number of rows in the Y matrix")
	}
	if(k > length(potentialControls))
	{
		stop("Input k cannot be larger than the number of negative controls")
	}
    	.Call('_RUVIIIC_RUVIIIC_Varying', PACKAGE = 'RUVIIIC', Y, k, M, potentialControls, toCorrect, withExtra, withW, withAlpha)
}

