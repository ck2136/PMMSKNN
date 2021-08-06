#' Generate table of training and testing (i.e., external) validation measures of bias, coverage and precision
#' 
#' @param loocvobj An object produced by \code{\link{loocv_function}}
#' 
#' @return A data frame that contains columns of performance measures. The performance measures include training and testing set measures of bias, coverage and precision as well as difference in  root mean square error (rmse), 50% coverage probability (Probability that predictions are included in the inter-quartile range (i.e. 25% to 75%)), and precision. The last column 'dominant' displays whether the prediction algorithm performed better in the training set or in the test set.
#' 
#' @export
extvalid <- function(loocvobj){
  
  res <- loocvobj$loocv_score %>%
    select(1:5) %>%
    mutate(train_test = "train") %>%
    filter(nearest_n == loocvobj$nearest_n) %>%
    bind_rows(
      loocvobj$test_score %>%
        mutate(train_test = "test")
    )  %>%
    tidyr::pivot_longer(rmse:prec) %>%
    select(3:5) %>%
    tidyr::pivot_wider(names_from =  train_test) %>%
    mutate(
      diff = train - test,
      pdiff = (train - test)/train * 100,
      dominant = case_when(
        name == "rmse" & train > test ~ "test",
        name == "rmse" & train < test ~ "train",
        name == "rmse" & train == test ~ "tie",
        
        name == "cov" & abs(train-0.5) > abs(test-0.5) ~ "test",
        name == "cov" & abs(train-0.5) < abs(test-0.5) ~ "train",
        name == "cov" & abs(train-0.5) == abs(test-0.5) ~ "tie",
        
        name == "prec" & train > test ~ "test",
        name == "prec" & train < test ~ "train",
        name == "prec" & train == test ~ "tie"
      )
    )
  

    
  res
}
