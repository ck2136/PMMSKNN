#' Match ID Function: create vector of ID values
#' 
#' @param ord_data Data frame. Specifically, training data with patient_id ordered based on fitted distal outcome value using predicted mean matching.
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$train_o} would be used for this parameter.
#' @param mtype - Integer value indicating matching type. Default is set to 1 which follows the
#' matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. \code{mtype} values are 
#' from \code{0} to \code{4}
#' @param loocv - Logical (\code{TRUE/FALSE}) that specifies whether 
#' or not to perform leave-one-out cross validation or just output 
#' predictions without hyperparameter tuning. If \code{loocv=FALSE}, then
#' users need to specify the value of the nearest_n 
#' @param n          Number of matches
#' @param m - For \code{mtype = 4}, which is type 4 matching from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}, the Number of repititions for obtaining \eqn{\dot{y}} in terms of the predictive mean matching process.
#' @param i          (Patient) Id indicator.
#' 
#' @return           A vector of patient id numbers
#' 
#' @export

matchIdExtract <- function(ord_data, 
                           mtype,  
                           loocv,
                           n, 
                           m,
                           i) {
    # FOR LOOCV
    if(loocv){
        if(mtype == 1){
            matches <- ord_data  %>% 
                bind_cols(diff = abs(ord_data$ydotavg[i] - ord_data$Fitted)) %>%
                .[-i,] %>%
                arrange(diff) %>%
                dplyr::select(.data$id) %>%
                head(n = n) %>% unlist %>% as.vector
            return(matches)
        } else if(mtype == 0) {
            matches <- ord_data  %>% 
                bind_cols(diff = abs(ord_data$Fitted[i] - ord_data$Fitted)) %>%
                .[-i,] %>%
                arrange(diff) %>%
                dplyr::select(.data$id) %>%
                head(n = n) %>% unlist %>% as.vector
            return(matches)
        } else if(mtype == 2) {
            matches <- ord_data  %>% 
                bind_cols(diff = abs(ord_data$ydotavg[i] - ord_data$ydotavg)) %>%
                .[-i,] %>%
                arrange(diff) %>%
                dplyr::select(.data$id) %>%
                head(n = n) %>% unlist %>% as.vector
            return(matches)
        } else if(mtype == 4) {
            matches <- sapply(1:m, function(x) {
                ord_data %>% 
                    select(matches(paste0("ydot","\\d")),.data$id,contains("Fitted"))  %>%
                    data.frame %>%
                    bind_cols(diff = abs(.[i,x] - .[,"Fitted"])) %>%
                    .[-i, ] %>%
                    arrange(diff) %>%
                    dplyr::select(.data$id) %>%
                    head(n=n)
            }) %>% unlist %>% as.vector
            return(matches)
        } else {
            stop("Matching Number not between 0-2 or 4")
        }
        
    # FOR TESTING
    } else {
        if(mtype == 1){
            matches <- ord_data  %>% 
                bind_cols(diff = abs(ord_data$ydotavg[i] - ord_data$Fitted)) %>%
                # .[-i,] %>%
                arrange(diff) %>%
                dplyr::select(.data$id) %>%
                head(n = n) %>% unlist %>% as.vector
            return(matches)
        } else if(mtype == 0) {
            matches <- ord_data  %>% 
                bind_cols(diff = abs(ord_data$Fitted[i] - ord_data$Fitted)) %>%
                # .[-i,] %>%
                arrange(diff) %>%
                dplyr::select(.data$id) %>%
                head(n = n) %>% unlist %>% as.vector
            return(matches)
        } else if(mtype == 2) {
            matches <- ord_data  %>% 
                bind_cols(diff = abs(ord_data$ydotavg[i] - ord_data$ydotavg)) %>%
                # .[-i,] %>%
                arrange(diff) %>%
                dplyr::select(.data$id) %>%
                head(n = n) %>% unlist %>% as.vector
            return(matches)
        } else if(mtype == 4) {
            matches <- sapply(1:m, function(x) {
                ord_data %>% 
                    select(matches(paste0("ydot","\\d")),.data$id,contains("Fitted"))  %>%
                    data.frame %>%
                    bind_cols(diff = abs(.[i,x] - .[,"Fitted"])) %>%
                    # .[-i, ] %>%
                    arrange(diff) %>%
                    dplyr::select(.data$id) %>%
                    head(n=n)
            }) %>% unlist %>% as.vector
            return(matches)
        } else {
            stop("Matching Number not between 0-2 or 4")
        }
        
    }
}
