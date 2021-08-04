#' Performance Metric Calculation Function: After LOOCV
#' 
#' @param loocv_res   - The \code{loocv_res} components of an object 
#' produced by \code{\link{loocv_function}}
#' @param outcome    - Name of the outcomes variable (type=string)
#' @param nearest_n - Numeric vector indicating the nearest number of matches 
#' to select from. The \code{loocv_function} will iterate through 
#' the number within the vector and select the number of matches that has
#' optimal bias, coverage, and precision.
#' The number of nearest_n should be within the range of number of individuals
#' That are in the data. For example if there are 500 individuals in the data,
#' one could do \code{nearest_n <- 10:100}
#' @param opt_cov - Float value indicating optimal coverage value used for `perfrank`.
#' @param perf_round_by - Integer value to indicate what decimal point will the performance values should be rounded by. Default is `perf_round_by = 4`, set to smaller value to be less coarse about ranking `nearest_n` values.
#' @param train - Logical indicating whether or not aggregating from training data or testing data. Deafult `TRUE`
#' 
#' @return Returns a data frame containing the performance measures (\code{bias}, \code{rmse}, \code{zscore}, \code{coverage}, \code{precision}, dropped cases due to failure in predicting values (\code{dcase}), number of matches (\code{nearest_n}), and the normalized scores for each of the performances measures (\code{bsc},\code{rsc},\code{covsc},\code{presc},\code{zsc}, and finally the sum of the scores: \code{totscore})) based on the number of nearest neighbors mathced (specified as \code{nearest_n}). 
#' 
#' @export

loocv_perf <- function(loocv_res, 
                       outcome,
                      nearest_n=nearest_n,
                      opt_cov = 0.5,
                      perf_round_by=perf_round_by,
                      train=TRUE
                      ){
    
    if(train){
        # For Training Data aggregation ---------------------
        
        perflist <- lapply(loocv_res, function(x) {
            
            #IF error object is returned then return NA columns
            if(is.null(x$pred_train)){
               return(
                      data.frame(
                                 rmse = NA,
                                 cov = NA,
                                 prec = NA
                      )
               )
            }
            resdf <- rbindlist(x$pred_train)
            prec <- c75 <- c25 <- cov <- NULL
            resdf[, cov := ifelse(get(outcome) >= c25 & get(outcome) <= c75, 1, 0)]
            resdf[, prec := c75 - c25]
            
            # get rmse
            # sqrt(sum((rbindlist(x$pred_train)[,outcome] - rbindlist(x$pred_train)[,"c50"])^2))
            data.frame(
                rmse = sqrt(sum((resdf[,outcome, with = FALSE] - resdf[,"c50"])^2)/nrow(resdf)),
                cov = mean(resdf$cov),
                prec = mean(resdf$prec)
            )
        })

        ## For cases where there are no data due to errors in fitting 
        #if(length(perflist) == 0){
           #return(
                  #data.frame(
                             #rmse = rep(NA,  length(nearest_n)),
                             #cov = rep(NA,  length(nearest_n)),
                             #prec = rep(NA,  length(nearest_n)),
                             #nearest_n = nearest_n,
                             #prec = rep(NA,  length(nearest_n))
                  #)
           #)
        #}
        
        perfdf <- rbindlist(perflist) %>%
            bind_cols(
                nearest_n = nearest_n
            ) %>%
            mutate(
                covdiff = round(abs(.data$cov - opt_cov), perf_round_by),
                rsc = (.data$rmse - min(.data$rmse, na.rm=TRUE)) / (max(.data$rmse, na.rm=TRUE) - min(.data$rmse, na.rm=TRUE)),
                covsc = if((max(.data$covdiff, na.rm=TRUE) - min(.data$covdiff, na.rm=TRUE)) == 0) { return(0)} else { (.data$covdiff - min(.data$covdiff, na.rm =TRUE)) / (max(.data$covdiff, na.rm=TRUE) - min(.data$covdiff, na.rm=TRUE))},
                presc = (.data$prec - min(.data$prec, na.rm=TRUE)) / (max(.data$prec, na.rm=TRUE) - min(.data$prec, na.rm=TRUE)),
            ) %>%
            bind_cols(
                dcase = sapply(loocv_res, function(x) {
                    c(x$dropped_cases)
                })
            ) %>%
            mutate(totscore = .data$rsc + .data$covsc + .data$presc + .data$dcase)  
        
        
    } else {

        # For Testing Data aggregation ---------------------
       if(is.null(loocv_res$pred_test)){
          resdf <- rbindlist(loocv_res[[1]]$pred_test)
       } else {
          resdf <- rbindlist(loocv_res$pred_test)
       }
        prec <- c75 <- c25 <- cov <- NULL
        resdf[, cov := ifelse(get(outcome) >= c25 & get(outcome) <= c75, 1, 0)]
        resdf[, prec := c75 - c25]
            
        perfdf <- data.frame(
            rmse = sqrt(sum((resdf[,outcome, with = FALSE] - resdf[,"c50"])^2)/nrow(resdf)),
            cov = mean(resdf$cov),
            prec = mean(resdf$prec)
        )
        
        perfdf <- perfdf %>%
            bind_cols(
                nearest_n = nearest_n
            ) %>%
            mutate(
                covdiff = round(abs(.data$cov - opt_cov), perf_round_by),
            ) 
        
    }
    
    

    return(perfdf)
}
