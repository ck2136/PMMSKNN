#' Match ID Function: create vector of ID values
#' 
#' @param ord_data   dataset. Training data from \code{\link{preproc}}. 
#' @param mtype       Integer value indicating matching type. Default is set to 1 which follows the
#'  matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. 
#' @param n          Number of matches
#' @param m          Number of repititions of obtaining \eqn{\dot{y}}
#' @param i          Patient id indicator.
#' @return           A vector of patient id numbers
#' @export

matchIdExtract <- function(ord_data, 
                           mtype,  
                           n, 
                           m,
                           i) {
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
}
