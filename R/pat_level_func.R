#' Calculate bias, coverage and CI width for neighbors-based prediction
#' 
#' The function calculate three frequency-based parameters to 
#' demonstate statistical quality of the neighbors-based prediction.  #' The parameters are bias, coverage and the 50% prediction interval 
#' width.
#' 
#' @param ref   - Fitted gamlss object from \code{\link{fitrefgamlss}}
#' @param nearest - Numeric vector indicating the nearest number of matches 
#' to select from. The \code{loocv_function} will iterate through 
#' the number within the vector and select the number of matches that has
#' optimal bias, coverage, and precision.
#' The number of nearest_n should be within the range of number of individuals
#' That are in the data. For example if there are 500 individuals in the data,
#' one could do \code{nearest_n <- 10:100}
#' @param dist_fam  - gamlss distribution specification using the \code{\link{gamlss.dist}} package. The specification for a normal distribution would be \code{gamlss.dist::NO}. For other distributions see \code{\link{gamlss.dist}}.
#' @param train_post - Data frame that contains the post-baseline observations from the training dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param ord_data Data frame. Specifically, training data with patient_id ordered based on fitted distal outcome value using predicted mean matching.
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$train_o} would be used for this parameter.
#' @param test_post Data frame that contains the post-baseline observations from the testing dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param test_o  Data frame. Specifically, testing data with patient_id ordered based on fitted distal outcome value using predicted mean matching.
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$test_o} would be used for this parameter.
#' @param traintestmatchdf  - Matched train test dataframe based on \code{\link{matchTestDataGen}}
#' @param outcome    - Name of the outcomes variable (type=string)
#' @param time_elapsed - Name of the time variable. (type=string)
#' @param plot - Logical (\code{TRUE/FALSE}) that specifies whether to output individual precision plots
#' @param matchprobweight - Logical (\code{TRUE/FALSE}) that specifies whether to utilize probability sampling
#'  when doing the mean matching. If TRUE, matches nearest n weighted on differnce in 
#'  predicted outcome.
#' @param time_window - vector of numbers for `centiles.pred()`, `xvalues` argument. For example, specify such as \code{c(10:30)}
#' @param interval - Int value specifying the interval of individuals to skip 
#' @param cs  - Logical that specifies whether to use cubic spline. 
#' The default is set to \code{cs = FALSE}. 
#' @param dfspec - Logical (\code{TRUE/FALSE}) that specifies whether to 
#' specify degrees of freedoms for the location, scale, and shape parameters
#' for the distribution specified with \code{dist_fam}.
#' Default value is \code{NULL}.
#' @param d_f_m - Numeric value that specifies the degrees of freedom for the cubic spline specified for the mean parameter of the distribution specified according to \code{dist_fam}
#' @param ptr_m - Numeric value that specifies the power transformation of time variable. Default value is 1.
#' @param d_f_s - Numeric value that specifies the degrees of freedom for the cubic spline specified for the scale parameter of the distribution specified according to \code{dist_fam}
#' @param d_f_n - Numeric value that specifies the degrees of freedom for the cubic spline specified for the shape parameter, specifically the \eqn{\nu} parameter, of the distribution specified according to \code{dist_fam}
#' @param d_f_t - Numeric value that specifies the degrees of freedom for the cubic spline specified for the shape parameter, specifically the \eqn{\tau} parameter, of the distribution specified according to \code{dist_fam}
#' @param thresh_val - Numeric value indicating value of bias to ignore (not include in output) in terms of the leave-one-out cross validation process. The default is set to \code{thresh_val = 10000}
#' @param printtrace - Logical (\code{TRUE/FALSE}) that specifies printing of gamlss parameter estimation process within the \code{\link{gamlss}} function
#' @param userchoose - Int value indicating the choice that the user wants to use for the number of nearest matches 
#' @param seed - Seed for probability sampling of the nearest matches
#' @param parallel - Number of cores used for the leave-one-out cross validation process. Default = 1
#' @param loocv - Logical (\code{TRUE/FALSE}) that specifies whether 
#' or not to perform leave-one-out cross validation or just output 
#' predictions without hyperparameter tuning. If \code{loocv=FALSE}, then
#' users need to specify the value of the nearest_n 
#' @param mtype - Integer value indicating matching type. Default is set to 1 which follows the
#' matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. \code{mtype} values are 
#' from \code{0} to \code{4}
#' @param biasm - Column indicating which bias measure to use for 
#' choosing the optimal nearest number of neighbors. 
#' Default is \code{'raw'}. Options: \code{'raw','rmse','zsc'}.
#' @param m - For \code{mtype = 4}, which is type 4 matching from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}, the Number of repititions for obtaining \eqn{\dot{y}} in terms of the predictive mean matching process.
#' 
#' @return Returns a list of 3 lists and a value. 1) \code{pred_res} contains a list of predicted values for the training data (\code{pred_train}) and test data (\code{pred_test}), the performance (\code{bias},\code{rmse},\code{zscore},\code{iqrcoverage},\code{precisionvec}, and number of dropped cases in fitting the gamlss model(\code{dropped_cases})); 2) \code{loocv_res} contains the same lists described above for each models fitted using different values of number of nearest neighbors; 3) \code{loocv_score} contains the summarised performance measures as a data frame; 4) \code{nearest_n} contains the optimal number of matches based on the aggregate performance metric.
#' 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# LOOCV Function ----
# - - - - - - - - - - - - - - - - - - - -#
pat_level_func <- function(
                           ref=ref,nearest, 
                           dist_fam = dist_fam, # for gamlss distribution
                           train_post=train_post, 
                           ord_data=ord_data, 
                           test_post=test_post, 
                           test_o=test_o,  # datasets
                           traintestmatchdf=traintestmatchdf,
                           outcome=outcome, time_elapsed=time_elapsed, 
                           plot = plot,
                           matchprobweight=matchprobweight,
                           time_window=time_window,
                           interval=interval,
                           cs=cs,
                           dfspec=dfspec,
                           d_f_m=d_f_m, ptr_m=ptr_m,
                           d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                           thresh_val = thresh_val,
                           printtrace=printtrace,
                           userchoose=userchoose,
                           seed=seed,
                           parallel=parallel,
                           loocv=loocv,
                           mtype=mtype,
                           biasm=biasm,
                           m=m
                           ) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # EXTRACT REF OBJECT MATERIAL
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    mint <- ref$mint; maxt <- ref$maxt; gamlss_dist <- ref$gamlss_dist; 
    spl <- ref$spl; spls <- ref$spls; spln <- ref$spln; splt <- ref$splt; ref <- ref$ref;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # PARALLEL/NON-PARALLEL IF THEN 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if(!is.null(parallel)){
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Parallel solution
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        registerDoParallel(cores = parallel)
        #registerDoSNOW(makeCluster(parallel, type = "SOCK"))
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Loop based on the nearest_n specified as user input 
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        loocvres <- foreach(
                            n=nearest,
                            .export=c("ord_data", "train_post", "test_post","spl","ptr_m", "d_f_m", "d_f_s",
                                      "d_f_n", "d_f_t", "dist_fam", "spls","spln","splt",
                                      "outcome","time_elapsed", "plot", "matchprobweight","time_window",
                                      "interval","thresh_val","printtrace","userchoose","seed", "ref",
                                      "loocv","mint","maxt", "traintestmatchdf", "m"),                                                                                                              .packages=c("dplyr","gamlss"),
                            .combine=list,
                            .multicombine = TRUE
                            ) %dopar% {
            #for(n in nearest){
            #cnt = 0
            misses = 0

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Instantiate array of values
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            perfout <- list(
                            dfList = vector("list", length = length(ord_data$id)), 
                            dfList_test = vector("list", length = length(ord_data$id)), 
                            centilepred = vector("list", length = length(ord_data$id)), 
                            biasvec = vector("numeric", length = length(ord_data$id)), 
                            rmsevec = vector("numeric", length = length(ord_data$id)), 
                            coveragevec = vector("numeric", length = length(ord_data$id)), 
                            coveragevec95 = vector("numeric", length = length(ord_data$id)), 
                            iqrvec = vector("numeric", length = length(ord_data$id)), 
                            precisionvec = vector("list", length = length(ord_data$id)), 
                            crazymatch = vector("list", length = length(ord_data$id)), 
                            zscore = vector("list", length = length(ord_data$id))
            )

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # iterate through all patients and match patients according to order. ord_data is the training data
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if(is.null(interval)){
                patlistint <- seq(1,length(ord_data$id))
            } else {
                patlistint <- seq(1,length(ord_data$id), by=interval)
            }
            for (i in patlistint) {

                # - - - - - - - - - - - - - - - - - - - - - - #
                # Matching cases using absolute value of difference between Fitted values
                # - - - - - - - - - - - - - - - - - - - - - - #
                matchmodel <- matchTrainDataGen(train_post=train_post, 
                                                ord_data=ord_data, 
                                                mtype=mtype, n=n,m=m, i=i,
                                                time_elapsed=time_elapsed,
                                                outcome = outcome,
                                                seed = seed,
                                                matchprobweight = matchprobweight)


                # - - - - - - - - - - - - - - - - - - - - - - #
                # Fit GAMLSS model to nearest n matches
                # - - - - - - - - - - - - - - - - - - - - - - #
                plmr <- update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE))

                # - - - - - - - - - - - - - - - - - - - - - - #
                # Error Check if plmr object is updated properly
                # - - - - - - - - - - - - - - - - - - - - - - #
                # if something wrong with iteration... then just don't do prediction it eats up all the memory anyways
                if(typeof(plmr) != 'list') {
                    misses = misses + 1 # increment number of misses
                    message('Something wrong with plm model. Prediction not included')
                } else {
                    #message(paste0("Getting predictions: "))
                    
                    invisible(perfout <- plmout(
                           plmr=plmr, # updated patient level model
                           i=i, time_window=time_window, mint=mint,maxt=maxt,
                           perfout = perfout, 
                           loocv=loocv,matchmodel=matchmodel,
                           traintestmatchdf=traintestmatchdf,
                           outcome=outcome, time_elapsed=time_elapsed, 
                           ord_data=ord_data, train_post=train_post, 
                           test_post=test_post, thresh_val = thresh_val
                    ))


                }
                #message(paste0("Current count is: ",i))
                #message(paste0("and Current n: ",n))
            }

            # - - - - - - - - - - - - - - - - - - - - - - - - -#  
            # Include missing points and rename and return
            # - - - - - - - - - - - - - - - - - - - - - - - - -#  
            perfout <- c(perfout, dropped_cases=misses)
            names(perfout) <- c("pred_train","pred_test",
                                "centilerange","bias",
                                "rmse","iqrcoverage",
                                "coverage95c","iqr",
                                "precisionvec","crazymatch", 
                                "zscore","dropped_cases")
            perfout

        }

        # return either LOOCV array or final prediction array
        return(loocvres)

    } else {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Non parallel solution
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Empty list that will be populated with performance measures and predictions
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nn_arr = vector("list", length=length(nearest))

        for(n in nearest){
            #cnt = 0
            misses = 0

            perfout <- list(
                            dfList = vector("list", length = length(ord_data$id)), 
                            dfList_test = vector("list", length = length(ord_data$id)), 
                            centilepred = vector("list", length = length(ord_data$id)), 
                            biasvec = vector("numeric", length = length(ord_data$id)), 
                            rmsevec = vector("numeric", length = length(ord_data$id)), 
                            coveragevec = vector("numeric", length = length(ord_data$id)), 
                            coveragevec95 = vector("numeric", length = length(ord_data$id)), 
                            iqrvec = vector("numeric", length = length(ord_data$id)), 
                            precisionvec = vector("list", length = length(ord_data$id)), 
                            crazymatch = vector("list", length = length(ord_data$id)), 
                            zscore = vector("list", length = length(ord_data$id))
            )

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # iterate through all patients and match patients according to order. ord_data is the training data
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if(is.null(interval)){
                patlistint <- seq(1,length(ord_data$id))
            } else {
                patlistint <- seq(1,length(ord_data$id), by=interval)
            }

            for (i in patlistint) {

                # - - - - - - - - - - - - - - - - - - - - - - #
                # Matching cases using absolute value of difference between Fitted values
                # Matching by probability weighting also possible
                # - - - - - - - - - - - - - - - - - - - - - - #

                matchmodel <- matchTrainDataGen(train_post=train_post, 
                                                ord_data=ord_data, 
                                                mtype=mtype, n=n,m=m, i=i,
                                                time_elapsed=time_elapsed,
                                                outcome = outcome,
                                                seed = seed,
                                                matchprobweight = matchprobweight)

                # - - - - - - - - - - - - - - - - - - - - - - #
                # Fit GAMLSS model to nearest n matches
                # - - - - - - - - - - - - - - - - - - - - - - #

                plmr <- tryCatch(
                                 {
                                     update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE))
                                 },
                                 error=function(cond){
                                     message(cond)
                                     return(1)
                                 },
                                 warning = function(cond){
                                     message(cond)
                                     update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE))
                                 },
                                 finally = {
                                 }
                )
                #plmr <- update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE))

                # if something wrong with iteration... then just don't do prediction it eats up all the memory anyways
                if(typeof(plmr) != 'list') {
                    misses = misses + 1 # increment number of misses
                    message('Something wrong with plm model. Prediction not included')
                } else {
                    #message(paste0("Getting predictions: "))

                    invisible(perfout <- plmout(
                           plmr=plmr, # updated patient level model
                           i=i, time_window=time_window, mint=mint,maxt=maxt,
                           perfout = perfout, 
                           loocv=loocv,matchmodel=matchmodel,
                           traintestmatchdf=traintestmatchdf,
                           outcome=outcome, time_elapsed=time_elapsed, 
                           ord_data=ord_data, train_post=train_post, 
                           test_post=test_post, thresh_val = thresh_val
                    ))

                }
                #message(paste0("Current count is: ",i))
                #message(paste0("and Current n: ",n))

            }

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # NON LOOCVPrediction Vector Only
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            # - - - - - - - - - - - - - - - - - - - - - - - - -#  
            # Include missing points and rename and return
            # - - - - - - - - - - - - - - - - - - - - - - - - -#  
            perfout <- c(perfout, dropped_cases=misses)
            if(!loocv){
                names(perfout) <- c("pred_train","pred_test",
                                    "centilerange","biasvec",
                                    "rmse","iqrcoverage",
                                    "coverage95c","iqr",
                                    "precisionvec","crazymatch", 
                                    "zscore","dropped_cases")
                nn_arr <- perfout
            } else {
                names(perfout) <- c("pred_train","pred_test",
                                    "centilerange","bias",
                                    "rmse","iqrcoverage",
                                    "coverage95c","iqr",
                                    "precisionvec","crazymatch", 
                                    "zscore","dropped_cases")
                nn_arr[[which(nearest == n)]] <- perfout
            }
        }

        # return either LOOCV array or final prediction array
        return(nn_arr)

    }

}
