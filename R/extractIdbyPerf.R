#' Extract ID by Performance: extracts ID value from loocv object based on the percentile of bias value of individuals
#' 
#' `extractIdbyPerf` is a utility function that returns one single id number based on the 
#' `perc` (numeric) argument provided. This is useful to determine the n-th percentile person
#' in terms of the overall performance values 
#' (could be a combination of bias, coverage, and precision) obtained from the LOOCV process.
#' 
#' @param test_proc  object from \code{\link{preproc}}. 
#' @param perc       Floating point value indicating N/100 percentile ID to look for.
#' 
#' @return           A numeric value representing patient id number 
#' 
#' @export

extractIdbyPerf <- function(test_proc, perc){
    val <- test_proc$test_o %>%
        arrange(.data$pred) %>%
        .[,2] %>% unlist %>% as.vector %>%
        quantile(., perc) %>% as.vector

    id <- test_proc$test_o %>%
        filter(.data$pred == val) %>%
        dplyr::select(id) %>% unlist %>% as.vector
    return(id)
}
