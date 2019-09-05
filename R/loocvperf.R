#' Performance Metric Calculation Function: After LOOCV
#' 
#' @param loocv_res   - The \code{loocv_res} components of an object 
#' produced by \code{\link{loocv_function}}
#' @param train_o - Data frame. Specifically, training data with patient_id ordered based on fitted distal outcome value using predicted mean matching.
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$train_o} would be used for this parameter.
#' @param bias - Column indicating which bias score to use. Default
#' is \code{'raw'}. Options: \code{'raw','rmse','zsc'}.
#' @param nearest_n - Numeric vector indicating the nearest number of matches 
#' to select from. The \code{loocv_function} will iterate through 
#' the number within the vector and select the number of matches that has
#' optimal bias, coverage, and precision.
#' The number of nearest_n should be within the range of number of individuals
#' That are in the data. For example if there are 500 individuals in the data,
#' one could do \code{nearest_n <- 10:100}
#' 
#' @return Returns a data frame containing the performance measures (\code{bias}, \code{rmse}, \code{zscore}, \code{coverage}, \code{precision}, dropped cases due to failure in predicting values (\code{dcase}), number of matches (\code{nearest_n}), and the normalized scores for each of the performances measures (\code{bsc},\code{rsc},\code{covsc},\code{presc},\code{zsc}, and finally the sum of the scores: \code{totscore})) based on the number of nearest neighbors mathced (specified as \code{nearest_n}). 
#' 
#' @export

loocvperf <- function(loocv_res, 
                      train_o,
                      bias="raw",
                      nearest_n=nearest_n
                      ){
    # --
    # Merging zscore, avg abs bias, avg rmse into dataframe
    # -- 
    perfdf <- data.frame(
                         # BIAS = MAE
                         bias = sapply(loocv_res, function(x) { 
                                           mean(abs(x$bias))
}),
                         # RMSE
                         rmse = sapply(loocv_res, function(x) { 
                                           mean(x$rmse)
}),
                         # ZSC
                         zscore = sapply(loocv_res, function(x) { 
                                             rbindlist(x$zscore, idcol="id") %>%
                                                 left_join(
                                                           train_o %>% dplyr::select(.data$id) %>%
                                                               rename(patient_id = .data$id) %>%
                                                               cbind(., id=seq(1,nrow(train_o))),
                                                           by = "id"
                                                           )%>% 
                                             filter(!is.na(.data$zsc) & !is.infinite(.data$zsc)) %>%
                                             summarise(zscore = mean(.data$zsc))
                                         #unlist(x$zscore$zsc)[!is.na(unlist(x$zscore)) & !is.infinite(unlist(x$zscore))]
}) %>% unlist %>% as.vector,
                         # COV
                         coverage = sapply(loocv_res, function(x) { 
                                               mean(x$iqrcoverage)
}),
                         # PREC
                         precision = sapply(loocv_res, function(x) { 
                                                rbindlist(x$precisionvec, idcol = "id") %>%
                                                    left_join(
                                                              train_o %>% dplyr::select(.data$id) %>%
                                                                  rename(patient_id = .data$id) %>%
                                                                  cbind(., id=seq(1,nrow(train_o))),
                                                              by = "id"
                                                              )%>% 
                                                summarise(precision = mean(.data$prec))
                                            #mean(sapply(x$precisionvec, function(y) {
                                            #mean(y$prec)
                                            #}), na.rm=TRUE)
}) %>% unlist %>% as.vector,
                         # DROPPED CASES,
                         dcase = sapply(loocv_res, function(x) { 
                                            mean(x$dropped_cases)
}),
                         # NNCOL
                         nearest_n = nearest_n
                         ) %>% 
    # --
    # Standardizing 
    # -- 
    mutate(
           # Nearest N to Values
           #nearest_n = as.numeric(gsub("nearest_", "", .data$nearest_n)),
           # MINMAX
           bsc = (.data$bias - min(.data$bias, na.rm=TRUE)) / (max(.data$bias, na.rm=TRUE) - min(.data$bias, na.rm=TRUE)),
           rsc = (.data$rmse - min(.data$rmse, na.rm=TRUE)) / (max(.data$rmse, na.rm=TRUE) - min(.data$rmse, na.rm=TRUE)),
           covsc = ifelse((max(abs(.data$coverage - 0.50), na.rm=TRUE) - min(abs(.data$coverage - 0.50), na.rm=TRUE)) == 0,0, (abs(.data$coverage - 0.50) - min(abs(.data$coverage - 0.50), na.rm =TRUE)) / (max(abs(.data$coverage - 0.50), na.rm=TRUE) - min(abs(.data$coverage - 0.50), na.rm=TRUE))),
           presc = (.data$precision - min(.data$precision, na.rm=TRUE)) / (max(.data$precision, na.rm=TRUE) - min(.data$precision, na.rm=TRUE)),
           zsc = (.data$zscore - min(.data$zscore, na.rm=TRUE)) / (max(.data$zscore, na.rm=TRUE) - min(.data$zscore, na.rm=TRUE))
    ) 

    # --
    # Choice of Total Score
    # -- 
    if(bias == "raw"){
        perfdf <- perfdf %>%
            mutate(totscore = .data$bsc + .data$covsc + .data$presc + .data$dcase)  
    } else if (bias == "zsc"){
        perfdf <- perfdf %>%
            mutate(totscore = .data$zsc + .data$covsc + .data$presc + .data$dcase)  
    } else { # RMSE as bias
        perfdf <- perfdf %>%
            mutate(totscore = .data$rsc + .data$covsc + .data$presc + .data$dcase)  
    }
    #mutate(totscore = .data$zsc + .data$covsc + .data$presc + .data$dcase)  %>%
    #.[complete.cases(.),]

    return(perfdf)
}
