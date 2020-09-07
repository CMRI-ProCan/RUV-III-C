#' @title Set or get the number of OpenMP threads
#' @description Set or get the number of threads used by OpenMP for parallel regions. 
#' 
#' @details These functions mostly exist because of the CRAN check, which only allows 1 - 2 threads.
#' As there is a per-thread memory requirement for \code{RUVIII_C} and \code{RUVIII_C_Varying}, 
#' these functions can also be useful for controlling memory usage. However in that case, setting
#' environment variable OMP_NUM_THREADS also works. 
#' 
#' @param num The desired number of threads
#' @return In the case of \code{omp_get_num_threads}, the current number of threads. 
#' @rdname openmp
#' @export
omp_set_num_threads <- function(num) {
	RUVIIIC_omp_set_num_threads(num)
}

#' @rdname openmp
#' @export
omp_get_num_threads <- function() {
	RUVIIIC_omp_get_num_threads()
}

