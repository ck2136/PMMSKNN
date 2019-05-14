#' Calculate bias, coverage and CI width for neighbors-based prediction
#' 
#' The function calculate three frequency-based parameters to 
#' demonstate statistical quality of the neighbors-based prediction.  #' The parameters are bias, coverage and the 50% prediction interval 
#' width.
#' @param ref   Fitted gamlss object from \code{\link{fitrefgamlss}}
#' @param nearest Numeric vector with number of matches per scenario
#' @param dist_fam  gamlss distribution specification
#' @param train_post - datasets, typically the \code{train_post} list component 
#'  of the object produced by \code{\link{preproc}}.
#' @param ord_data Idem, component \code{train_o}
#' @param test_post Idem, component \code{test_post}
#' @param test_o  Idem, component \code{test_o}
#' @param traintestmatchdf  Matched train test dataframe based on \code{\link{matchTestDataGen}}
#' @param outcome    Name of the outcomes variable
#' @param time_elapsed - Name of the time variable. (type=string)
#' @param plot - Logical that specifies whether to output individual precision plots
#' @param matchprobweight - Logical that specifies whether to utilize probability sampling
#'  when doing the mean matching. If TRUE, matches nearest n weighted on differnce in 
#'  predicted outcome.
#' @param time_window - vector of numbers for `centiles.pred()`, `xvalues` argument 
#' @param interval - Int value specifying the interval of individuals to skip 
#' @param cs Logical that specifies whether to use cubic spline. 
#'  The default \code{cs = FALSE} uses ...
#' @param dfspec Logical that specifies whether to the user sets  
#' degrees of freedom (...not clear to me what it does, and how it 
#' interacts with the next set of arguments)
#' @param d_f_m ... explain: arguments probably mean different things for different distributions. Might be preferable to package it with \code{dist_fam}
#' @param ptr_m -
#' @param d_f_s -
#' @param d_f_n - 
#' @param d_f_t -
#' @param thresh_val -
#' @param printtrace - Logical that specifies printing of gamlss parameter estimation
#' @param userchoose - Int value indicating the choice that the user wants to use for the number of nearest matches 
#' @param seed - Seed for probability sampling of the nearest matches
#' @param parallel - Number of cores used for parallel computing. Default = 1
#' @param loocv - Whether or not to perform leave one out cross validation or just go straight to prediction. Should
#'  have the userchoose value specified if `loocv=FALSE`
#' @param mtype - Integer value indicating matching type. Default is set to 1 which follows the
#'  matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. 
#' @param biasm - Column indicating which bias score to use for choosing optimal n. Default
#' is \code{'raw'}. Options: \code{'raw','rmse','zsc'}.
#' @param m             Number of repititions of obtaining \eqn{\dot{y}}
#' @return There are many possible return values 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# LOOCV Function ----
# - - - - - - - - - - - - - - - - - - - -#
# now using that we need to come up with the right number of matches that will give us the best bias and coverage.
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
                                      "loocv","mint","maxt", "traintestmatchdf"),                                                                                                              .packages=c("dplyr","gamlss"),
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
                message(paste0("Current count is: ",i))
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
                                     return(NA)
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
                    message(paste0("Getting predictions: "))

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
                message(paste0("Current count is: ",i))
                message(paste0("and Current n: ",n))

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
