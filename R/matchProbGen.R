#' Match Probability Generating Function: create data frame of ID's with probability weights
#' 
#' @param ord_data   dataset. Training data from \code{\link{preproc}}. 
#' @param mtype       Integer value indicating matching type. Default is set to 1 which follows the
#'  matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. 
#' @param n          Number of matches
#' @param m          Number of repititions of obtaining \eqn{\dot{y}}
#' @param i          Patient id indicator.
#' @return           A vector of patient id numbers
#' @export
matchProbGen <- function(ord_data, 
                           mtype,  
                           n, 
                           m,
                           i) {

    if(mtype == 0){
        matchprob <- ord_data %>% 
            bind_cols(diff = abs(ord_data$Fitted[i] - ord_data$Fitted)) %>%
            .[-i,] %>%
            arrange(diff) %>%
            dplyr::select(.data$id, .data$diff) %>%
            head(n = n) %>% 
            mutate(weight = diff/sum(diff)) %>%
            dplyr::select(.data$id, .data$weight) %>%
            as.data.frame %>%
            rename(patient_id = 1)

    } else if(mtype == 1){
        matchprob <- ord_data %>% 
            bind_cols(diff = abs(ord_data$ydotavg[i] - ord_data$Fitted)) %>%
            .[-i,] %>%
            arrange(diff) %>%
            dplyr::select(.data$id, .data$diff) %>%
            head(n = n) %>% 
            mutate(weight = diff/sum(diff)) %>%
            dplyr::select(.data$id, .data$weight) %>%
            as.data.frame %>%
            rename(patient_id = 1)

    } else if(mtype == 2){
        matchprob <- ord_data %>% 
            bind_cols(diff = abs(ord_data$ydotavg[i] - ord_data$ydotavg)) %>%
            .[-i,] %>%
            arrange(diff) %>%
            dplyr::select(.data$id, .data$diff) %>%
            head(n = n) %>% 
            mutate(weight = diff/sum(diff)) %>%
            dplyr::select(.data$id, .data$weight) %>%
            as.data.frame %>%
            rename(patient_id = 1)

    } else if(mtype == 4){
        matchprob <- lapply(1:m, function(x) {
                                ord_data %>% 
                                    select(matches(paste0("ydot","\\d")),.data$id,contains("Fitted"))  %>%
                                    data.frame %>%
                                    bind_cols(diff = abs(.[i,x] - ord_data$Fitted))  %>% 
                                    .[-i, ] %>% 
                                    arrange(diff) %>%
                                    dplyr::select(.data$id, .data$diff) %>%
                                    head(n=n) %>%
                                    mutate(weight = diff/sum(diff)) %>%
                                    dplyr::select(.data$id, .data$weight) %>%
                                    as.data.frame %>%
                                    rename(patient_id = 1)
                           }) %>% rbindlist %>%
                                group_by(.data$patient_id) %>% 
                                summarise(weight = sum(.data$weight)) %>%
                                mutate(weight = .data$weight/sum(.data$weight))

                            return(matchprob)
    }
    return(matchprob)

}
