#' Extract ID by Bias: extracts ID value from loocv and percentile based on bias
#' 
#' @param loocvObj   object from \code{\link{loocv_function}}. 
#' @param perc       Floating point value indicating N/100 percentile ID to look for.
#' @return           A numeric value of patient id number 
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
