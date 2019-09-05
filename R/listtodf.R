#' Converts Performance Values From a List To Dataframe
#' 
#' Utility function to convert results from the LOOCV for all 
#' nearest neighbors from a list of bias, precision, and coverage
#' values into an aggregated summarised data frame for generating
#' tables and pltos 
#' 
#' @param listfile - The leave-one-out-cross-validation result list 
#' component created using the \code{loocv_function}. 
#' If \code{x <- loocv_function()}, then \code{x$loocv_res} would be 
#' the \code{listfile} to be inputted into the function.
#' 
#' @return A data frame that contains performance values (bias, coverage probabilities, zscore, rmse, and dropped cases) by the number of matches (per row)
listtodf <- function(listfile){
    temp <- lapply(listfile, function(x) {
                       c(mean(x$bias, na.rm=TRUE),
                         mean(x$iqrcoverage, na.rm=TRUE),
                         mean(x$coverage95c, na.rm=TRUE),
                         mean(rbindlist(x$zscore)$zsc[!is.na(rbindlist(x$zscore)$zsc) & !is.infinite(rbindlist(x$zscore)$zsc)]),
                         mean(x$iqr, na.rm=TRUE),
                         mean(x$rmse, na.rm=TRUE),
                         x$dropped_cases)
                      })
    temp <- as.data.frame(temp)
    colnames(temp) <- regmatches(names(listfile), regexpr("\\d+",names(listfile)))
    rownames(temp) <- paste0(c('bias','iqrcoverage','coverage95c','zscore','iqrdif','rmse','dropped_cases'))
    temp$measure  <-  c('bias','iqrcoverage','coverage95c','zscore','iqrdif',"rmse","dropped_cases")
    temp <- melt(temp, id.vars=c("measure"))
    temp <- temp %>%
        mutate(variable = as.numeric(as.character(.data$variable))) %>%
        rename(nearest_n = .data$variable)
    return(temp)
}
