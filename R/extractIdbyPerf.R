#' Extract ID by Bias: extracts ID value from loocv and percentile based on bias
#' 
#' @param test_proc  object from \code{\link{preproc}}. 
#' @param perc       Floating point value indicating N/100 percentile ID to look for.
#' @return           A numeric value of patient id number 
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
