#' Match ID Function: create vector of (Training set) ID values for Test set cases
#' 
#' @param test_proc   Object processed after using \code{\link{preproc}}. 
#' For example, \code{x <- preproc()}, then \code{test_proc = x}
#' @param mtype - Integer value indicating matching type. Default is set to 1 which follows the
#' matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. \code{mtype} values are 
#' from \code{0} to \code{4}
#' @param i          (Test patient) id indicator.
#' @param n          Number of matches
#' 
#' @return           A vector of (patient) id numbers
#' 
#' @export

matchIdExtractTest <- function(
                           test_proc,
                           mtype,  
                           i,
                           n
                           ) {
    if(mtype == 0){
        temp <- data.frame(train_id = sapply(1:length(test_proc$test_o$pred), 
                                             function(x) {
                                                 test_proc$train_o[which.min(abs(test_proc$test_o$pred[x] - test_proc$train_o$Fitted)), "id"]
                                             }), 
                           test_id = test_proc$test_o$id)

        matches <- test_proc$train_o %>%
            bind_cols(diff = abs(test_proc$train_o[test_proc$train_o$id == temp[temp$test_id == i, "train_id"], "Fitted"] - test_proc$train_o$Fitted)) %>%
            arrange(diff) %>%
            dplyr::select(.data$id)  %>%
            head(n = n) %>% unlist %>% as.vector

        return(matches)
    } else if(mtype == 1) {
        temp <- data.frame(train_id = sapply(1:length(test_proc$test_o$pred), 
                                             function(x) {
                                                 test_proc$train_o[which.min(abs(test_proc$test_o$ydotavg[x] - test_proc$train_o$Fitted)), "id"]
                                             }), 
                           test_id = test_proc$test_o$id)

        matches <- test_proc$train_o %>%
            bind_cols(diff = abs(test_proc$train_o[test_proc$train_o$id == temp[temp$test_id == i, "train_id"], "Fitted"] - test_proc$train_o$Fitted)) %>%
            arrange(diff) %>%
            dplyr::select(.data$id)  %>%
            head(n = n) %>% unlist %>% as.vector

        return(matches)
    } else if(mtype == 2) {
        temp <- data.frame(train_id = sapply(1:length(test_proc$test_o$pred), 
                                             function(x) {
                                                 test_proc$train_o[which.min(abs(test_proc$test_o$ydotavg[x] - test_proc$train_o$ydotavg)), "id"]
                                             }), 
                           test_id = test_proc$test_o$id)

        matches <- test_proc$train_o %>%
            bind_cols(diff = abs(test_proc$train_o[test_proc$train_o$id == temp[temp$test_id == i, "train_id"], "Fitted"] - test_proc$train_o$Fitted)) %>%
            arrange(diff) %>%
            dplyr::select(.data$id)  %>%
            head(n = n) %>% unlist %>% as.vector
        return(matches)
    } else if(mtype == 4) {
        temp <- data.frame(train_id = sapply(1:length(test_proc$test_o$pred), 
                                             function(x) {
                                                 test_proc$train_o[which.min(abs(test_proc$test_o$ydotavg[x] - test_proc$train_o$ydotavg)), "id"]
                                             }), 
                           test_id = test_proc$test_o$id)

        matches <- test_proc$train_o %>%
            bind_cols(diff = abs(test_proc$train_o[test_proc$train_o$id == temp[temp$test_id == i, "train_id"], "Fitted"] - test_proc$train_o$Fitted)) %>%
            arrange(diff) %>%
            dplyr::select(.data$id)  %>%
            head(n = n) %>% unlist %>% as.vector
        return(matches)
    } else {
        stop("Matching Number not between 1-3 or 5")
    }
}
