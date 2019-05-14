#' Match Testing Data Generating Function: create matched trainig and testing pair data
#' 
#' @param test_o        dataset. Testing data with patient_id ordered from \code{\link{preproc}}. 
#' @param ord_data      dataset. Training data with patient_id ordered from \code{\link{preproc}}. 
#' @param mtype         Integer value indicating matching type. Default is set to 1 which follows the
#'  matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. 
#' @return              A vector of patient id numbers
#' @export
matchTestDataGen <- function(
                              test_o=test_o, 
                              ord_data=ord_data, 
                              mtype=1
                              ){

    if(mtype == 0){
        matcheddf <- data.frame(train_id = sapply(1:length(test_o$pred), function(x) {
                                                 ord_data[which.min(abs(test_o$pred[x] - ord_data$Fitted)), "id"]
}), 
                           test_id = test_o$id)
    } else if(mtype == 1){
        matcheddf <- data.frame(train_id = sapply(1:length(test_o$ydotavg), function(x) {
                                                 ord_data[which.min(abs(test_o$ydotavg[x] - ord_data$Fitted)), "id"]
}), 
                           test_id = test_o$id)
    } else if(mtype == 2){
        matcheddf <- data.frame(train_id = sapply(1:length(test_o$ydotavg), function(x) {
                                                 ord_data[which.min(abs(test_o$ydotavg[x] - ord_data$ydotavg)), "id"]
}), 
                           test_id = test_o$id)
    } else if(mtype == 4){
        matcheddf <- data.frame(train_id = sapply(1:length(test_o$ydotavg), function(x) {
                                                 ord_data[which.min(abs(test_o$ydotavg[x] - ord_data$ydotavg)), "id"]
}), 
                           test_id = test_o$id)
    }

    return(matcheddf)
}
