#' Extract ID by Bias: extracts ID value from loocv object based on the percentile of bias value of individuals
#' 
#' `extractIdbyBias` is a utility function that returns one single id number based on the 
#' `perc` (numeric) argument provided. This is useful to determine the n-th percentile person
#' in terms of the bias values obtained from the LOOCV process.
#' 
#' @param loocvObj   object from \code{\link{loocv_function}}. 
#' @param perc       Floating point value indicating N/100 percentile ID to look for.
#' 
#' @return           A numeric value representing patient id number 
#' 
#' @export

extractIdbyBias <- function(loocvObj, perc){
    rbindlist(loocvObj$pred_res$pred_test) %>%
        rename(test_id = 1)  %>%
        group_by(.data$test_id) %>%
        summarise(bias = mean(abs(.data$c50 - .data$tug)))  %>%
        arrange(.data$bias) %>%
        .[,1]  %>% unlist %>% as.vector %>%
        quantile(., perc) %>% as.vector
}
