#' Match Training Data Generating Function: create Data based on pmm
#' 
#' @param train_post    dataset. Post operative training data from \code{\link{preproc}}. 
#' @param ord_data      dataset. Training data with patient_id ordered from \code{\link{preproc}}. 
#' @param mtype         Integer value indicating matching type. Default is set to 1 which follows the
#'  matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. 
#' @param n             Number of matches
#' @param m             Number of repititions of obtaining \eqn{\dot{y}}
#' @param i             Patient id indicator.
#' @param time_elapsed  Name of the time variable. (type=string)
#' @param outcome       Name of the outcomes variable
#' @param seed          Seed for randomly selecting matches based on differences in predicted y
#' @param matchprobweight Logical that specifies whether to utilize probability sampling
#' @return              A vector of patient id numbers
#' @export
matchTrainDataGen <- function(
                              train_post=train_post, 
                              ord_data, mtype=1, n, m = 5,i,
                              time_elapsed=time_elapsed,
                              outcome = outcome,
                              seed = seed,
                              matchprobweight){
    # If using new matching method where matches from all m subsets are extracted

    if(matchprobweight){
        matchprob <- matchProbGen(ord_data, mtype, n,m, i)
        matches <- matchIdExtract(ord_data, mtype, n,m, i) 
        # create dataset with weights
        matchindexdf <- data.frame(
                                   patient_id = matches
                                   ) %>% 
                                left_join(
                                          matchprob,
                                          by = "patient_id"
                                ) 
        #matchmodel <- train_post[train_post$patient_id %in% matches, ] %>%
            #left_join(
                      #matchprob,
                      #by = "patient_id"
                      #) %>%
        #dplyr::select_("patient_id", time_elapsed, outcome, "weight") 

        set.seed(seed)
        matchindex <- sample(x = seq_len(nrow(matchindexdf)), 
                             size=nrow(matchindexdf), 
                             prob=matchindexdf$weight, 
                             replace=TRUE)

        # create final weighted matchdataset
        matches <- matchindexdf[unique(matchindex),"patient_id"]
        matchmodel <- train_post[train_post$patient_id %in% matches, ]
        return(matchmodel)
    } else {
        matches <- matchIdExtract(ord_data, mtype, n,m, i) 
        matchmodel <- train_post[train_post$patient_id %in% matches, ]
        return(matchmodel)
    }
}
