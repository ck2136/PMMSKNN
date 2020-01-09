#' Match Training Data Generating Function: create Data based on pmm
#' 
#' @param train_post - Data frame that contains the post-baseline observations from the training dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param ord_data Data frame. Specifically, training data with patient_id ordered based on fitted distal outcome value using predicted mean matching.
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$train_o} would be used for this parameter.
#' @param mtype - Integer value indicating matching type. Default is set to 1 which follows the
#' matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. \code{mtype} values are 
#' from \code{0} to \code{4}
#' @param n             Number of matches
#' @param m - For \code{mtype = 4}, which is type 4 matching from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}, the Number of repititions for obtaining \eqn{\dot{y}} in terms of the predictive mean matching process.
#' @param i             (Patient) Id indicator.
#' @param time_elapsed  Name of the time variable. (type=string)
#' @param outcome       Name of the outcomes variable.
#' @param seed          Seed for randomly selecting matches based on differences in predicted y
#' @param loocv - Logical (\code{TRUE/FALSE}) that specifies whether 
#' or not to perform leave-one-out cross validation or just output 
#' predictions without hyperparameter tuning. If \code{loocv=FALSE}, then
#' users need to specify the value of the nearest_n 
#' @param matchprobweight Logical (\code{TRUE/FALSE}) that specifies whether to utilize probability sampling for selecting patients according to difference in the distance between the neighbor.
#' 
#' @return              A vector of (patient) id numbers
#' 
#' @export
matchTrainDataGen <- function(
                              train_post=train_post, 
                              ord_data, mtype=1, n, m = 5,i,
                              time_elapsed=time_elapsed,
                              outcome = outcome,
                              seed = seed, loocv=loocv,
                              matchprobweight){
    # If using new matching method where matches from all m subsets are extracted

    if(matchprobweight){
        matchprob <- matchProbGen(ord_data, mtype, n,m, i)
        matches <- matchIdExtract(ord_data, mtype, loocv=loocv, n,m, i) 
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
        matches <- matchIdExtract(ord_data, mtype, loocv=loocv, n,m, i) 
        matchmodel <- train_post[train_post$patient_id %in% matches, ]
        return(matchmodel)
    }
}
