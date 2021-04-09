library(testthat)
library(RUVIIIC)

try(RUVIIIC::omp_set_num_threads(2L), silent = TRUE)
test_check("RUVIIIC")
