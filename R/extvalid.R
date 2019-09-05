#' Calculate external validation measures of bias, coverage and precision
#' 
#' @param loocvobj An object produced by \code{\link{loocv_function}}
#' @param test_proc An object produced by \code{\link{preproc}}
#' 
#' @return A data frame that contains columns of performance measures. The performance measures include difference in prediction (bias), root mean square error (rmse), difference in z-score (zscore), 50% coverage probability (Probability that predictions are included in the inter-quartile range (i.e. 25% to 75%))
#' 
#' @export
extvalid <- function(loocvobj,
                     test_proc){
  
  # for matching train and test patients via min value difference in Fitted of train (i.e. pmm)
  temp <- data.frame(train_id = sapply(1:length(test_proc$test_o$pred), function(x) {
    #which.min(abs(test_proc$test_o$pred[x] - test_proc$train_o$Fitted))
    test_proc$train_o[which.min(abs(test_proc$test_o$pred[x] - test_proc$train_o$Fitted)),"id"]
  }), 
  test_id = test_proc$test_o$id)
  precisionvec = test_proc$test_post %>% 
    dplyr::select(.data$patient_id, .data$time) %>% 
    rename(test_id = .data$patient_id) %>% 
    left_join(
      temp ,
      by = "test_id"
    ) %>% 
    left_join(
      # merge training data precision by time and training id
      rbindlist(lapply(loocvobj$pred_res$precisionvec, data.frame), idcol="id")  %>%
        left_join(
          # create a table that has ordered training id and sequence from 1 to 398 because the precisionvector doesn't have id's recorded
          test_proc$train_o %>% dplyr::select(.data$id) %>%
            rename(train_id = .data$id) %>%
            cbind(., id=seq(1,nrow(test_proc$train_o))) ,
          by = "id"
        ) ,
      by = c("time","train_id")
    ) %>%
    dplyr::select(.data$prec) %>% unlist %>% as.vector
  
  res <- data.frame(
                    bias = mean(abs(loocvobj$pred_res$bias), na.rm=TRUE),
                    rmse = mean(loocvobj$pred_res$rmse, na.rm=TRUE),
                    zscore = mean(rbindlist(loocvobj$pred_res$zscore)$zsc[!is.na(rbindlist(loocvobj$pred_res$zscore)$zsc) & !is.infinite(rbindlist(loocvobj$pred_res$zscore)$zsc)]),
                    coverage = mean(loocvobj$pred_res$iqrcoverage, na.rm=TRUE),
                    precision = test_proc$test_post %>% 
                      dplyr::select(.data$patient_id, .data$time) %>% 
                      rename(test_id = .data$patient_id) %>% 
                      left_join(
                        temp ,
                        by = "test_id"
                      ) %>% 
                      left_join(
                        # merge training data precision by time and training id
                        rbindlist(lapply(loocvobj$pred_res$precisionvec, data.frame), idcol="id")  %>%
                          left_join(
                            # create a table that has ordered training id and sequence from 1 to 398 because the precisionvector doesn't have id's recorded
                            test_proc$train_o %>% dplyr::select(.data$id) %>%
                              rename(train_id = .data$id) %>%
                              cbind(., id=seq(1,nrow(test_proc$train_o))) ,
                            by = "id"
                          ) ,
                        by = c("time","train_id")
                      )  %>% 
                      # here we use na.rm but really need to think about if there are any dropped predictions.. probably not a good idea to use the model at all
                      summarise(precision = mean(.data$prec, na.rm=TRUE))
  )
  if(any(is.na(precisionvec))){
    warning("NA values in precision")
  }
  if(any(is.na(rbindlist(loocvobj$pred_res$zscore)$zsc))){
    warning("NA values in zscore")
  }
  if(any(is.na(loocvobj$pred_res$iqrcoverage))){
    warning("NA values in coverage")
  }
  res
}
