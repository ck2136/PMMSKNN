#' Convert Performance metrics into dataframe RMSE and COVERAGE
#' 
#' Utility function to aggregate all nearest_n results 
#' for all people.
#' @param listfile DOCUMENT
#' @return DOCUMENT
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
