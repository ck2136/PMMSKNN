#' Calculate bias, coverage and CI width for neighbors-based prediction
#' 
#' The function calculate three frequency-based parameters to 
#' demonstate statistical quality of the neighbors-based prediction.
#' The parameters are bias, coverage and the 50% prediction interval 
#' width.
#' @param nearest_n Numeric vector with number of matches per scenario
#' @param dist_fam  gamlss distribution specification
#' @param train_post - datasets, typically the \code{train_post} list component 
#'  of the object produced by \code{\link{preproc}}.
#' @param ord_data Idem, component \code{train_o}
#' @param test_post Idem, component \code{test_post}
#' @param test_o  Idem, component \code{test_o}
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
#' @param ptr_m - Numeric value for any power transformation of time variable.
#' For example a \eqn{time^2} would be done as \code{ptr_m = 2}.
#' @param d_f_s - Numeric value to impose a Degree of freedom on the spline term for
#' the \eqn{\sigma} parameter.
#' @param d_f_n - Numeric value to impose a Degree of freedom on the spline term for
#' the \eqn{\nu} parameter.
#' @param d_f_t - Numeric value to impose a Degree of freedom on the spline term for
#' the \eqn{\tau} parameter.
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
#' @param \dots Passed down to \code{gamlss}
#' @return There are many possible return values 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# LOOCV Function ----
# - - - - - - - - - - - - - - - - - - - -#
# now using that we need to come up with the right number of matches that will give us the best bias and coverage.
loocv_function <- function(nearest_n = seq(20,150,by=10), # number to play with 
                           dist_fam = NULL, # for gamlss distribution
                           train_post, 
                           ord_data, 
                           test_post, 
                           test_o,  # datasets
                           outcome, time_elapsed, plot = FALSE,
                           matchprobweight=FALSE,
                           time_window=NULL,
                           interval=NULL,
                           cs=FALSE,
                           dfspec=NULL,
                           d_f_m=1.8, ptr_m=1,
                           d_f_s=1.2,
                           d_f_n=1,
                           d_f_t=1,
                           thresh_val = 10000,
                           printtrace=FALSE,
                           userchoose=NULL,
                           seed=1234,
                           parallel=NULL,
                           loocv=TRUE,
                           mtype=1,
                           biasm="raw",
                           m,
                           ...) {

    
    # - - - - - - - - - - - - - - - - - - - - - # 
    # FIT REFERENCE GAMLSS MODEL
    # - - - - - - - - - - - - - - - - - - - - - # 

    ref <- fitrefgamlss(
                        dist_fam = dist_fam, # for gamlss distribution
                        train_post=train_post, 
                        test_post=test_post, 
                        outcome=outcome, time_elapsed=time_elapsed, 
                        time_window=time_window,
                        cs=cs,
                        dfspec=dfspec,
                        d_f_m=d_f_m, ptr_m=ptr_m,
                        d_f_s=d_f_s,
                        d_f_n=d_f_n,
                        d_f_t=d_f_t,
                        ...) 

    # - - - - - - - - - - - - - - - - - - - - - # 
    # NEAREST NEIGHBOR MATCHING
    # - - - - - - - - - - - - - - - - - - - - - # 

    # Outerloop will iterate through a list of nearest_n numbers
    # Innerloop will iterate through all patients and store the result of the whole iteration into an array
    # nearest neighbor should be an even number for at least those that aren't edge cases

    
    # iterate through the nearest_n list calculate bias coverage for each iteration
    # for e.g. if we go form 10:100 by 10's we will have a list with 10 results 

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Match training and testing patients based on minimum difference in the predicted outcome at out_time 
    # Need this to get test data predictions. Default matching is type 2 matching
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    traintestmatchdf <- matchTestDataGen(
                                         test_o=test_o, 
                                         ord_data=ord_data, 
                                         mtype=mtype
    )


    # - - - - - - - - - - - - - - - - - - - - - - #
    # Calculation of Weighted Z-score, Coverage, and Bias to select Optimal N #
    # - - - - - - - - - - - - - - - - - - - - - - #

    if(length(nearest_n) != 1){

        if(loocv){
            loocv_test_result <- pat_level_func(
                                                ref=ref,nearest=nearest_n, # number to play with 
                                                dist_fam = dist_fam, # for gamlss distribution
                                                train_post=train_post, ord_data=ord_data, 
                                                test_post=test_post, test_o=test_o,  # datasets
                                                traintestmatchdf=traintestmatchdf,
                                                outcome=outcome, time_elapsed=time_elapsed, 
                                                plot = plot,matchprobweight=matchprobweight,
                                                time_window=time_window, interval=interval,
                                                cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                                d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                                thresh_val=thresh_val, printtrace=printtrace,
                                                userchoose=userchoose, seed=seed,
                                                parallel=parallel, loocv=loocv, mtype=mtype,
                                                biasm=biasm,m=m
            )
            # - - - - - - - - - - - - - - - - - - - - - - #
            # Allow user to choose optimal number of nearest neighbor
            # - - - - - - - - - - - - - - - - - - - - - - #
            if(!is.null(userchoose)){
                if(!any(nearest_n %in% userchoose)){
                    stop("Optimal n not in the range of nearest_n specified!")
                }
                opt_n_index <- which(nearest_n == userchoose)
            } else {
                # - - - - - - - - - - - - - - - - - - - - - - #
                # Use a weighted scoring fucntion of bias, coverage, and precision to choose optimal m
                # - - - - - - - - - - - - - - - - - - - - - - #
                opt_n_index <- loocvperf(loocv_test_result, ord_data, bias=biasm, nearest_n) %>%
                    dplyr::select(.data$totscore) %>% unlist(.) %>% which.min(.) %>% as.vector(.)
            }
        }
        names(loocv_test_result) <- paste0("nearest_",nearest_n)

        print(paste0(mean(loocv_test_result[[opt_n_index]]$rmse, na.rm=TRUE), " from ", names(loocv_test_result))[[opt_n_index]])
        print(paste0("Number of misses is: ",loocv_test_result[[opt_n_index]]$dropped_cases,' cases'))
        print(paste0("Distribution chosen for matched GAMLSS: ", ref$gamlss_dist))
        print(paste0("Optimal Number of Matches is: ", nearest_n[[opt_n_index]]))
        # - - - - - - - - - - - - - - - - - - - - - - #
        # Run the result on test set
        # - - - - - - - - - - - - - - - - - - - - - - #
        predict_test_result <- pat_level_func(
                                              ref=ref,nearest=nearest_n[opt_n_index], # number to play with 
                                              dist_fam = dist_fam, # for gamlss distribution
                                              train_post=train_post, ord_data=ord_data, 
                                              test_post=test_post, test_o=test_o,  # datasets
                                              traintestmatchdf=traintestmatchdf,
                                              outcome=outcome, time_elapsed=time_elapsed, 
                                              plot = plot,matchprobweight=matchprobweight,
                                              time_window=time_window, interval=interval,
                                              cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                              d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                              thresh_val = thresh_val, printtrace=printtrace,
                                              userchoose=userchoose, seed=seed,
                                              parallel=parallel, loocv=FALSE, mtype=mtype,
                                              biasm=biasm,m=m
        ) 
        # extract prediction results from the optimum nearest n 
        return(list(pred_res = predict_test_result,
                    loocv_res =  loocv_test_result,
                    loocv_score = loocvperf(loocv_test_result, ord_data, bias=biasm, nearest_n),
                    nearest_n=nearest_n[opt_n_index]))
    } else {
        if(loocv){

            loocv_test_result <- pat_level_func(
                                                ref=ref,nearest=nearest_n, # number to play with 
                                                dist_fam = dist_fam, # for gamlss distribution
                                                train_post=train_post, ord_data=ord_data, 
                                                test_post=test_post, test_o=test_o,  # datasets
                                                traintestmatchdf=traintestmatchdf,
                                                outcome=outcome, time_elapsed=time_elapsed, 
                                                plot = plot,matchprobweight=matchprobweight,
                                                time_window=time_window, interval=interval,
                                                cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                                d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                                thresh_val = thresh_val, printtrace=printtrace,
                                                userchoose=userchoose, seed=seed,
                                                parallel=parallel, loocv=loocv, mtype=mtype,
                                                biasm=biasm,m=m
            ) 
            predict_test_result <- pat_level_func(
                                                  ref=ref,nearest=nearest_n[opt_n_index], # number to play with 
                                                  dist_fam = dist_fam, # for gamlss distribution
                                                  train_post=train_post, ord_data=ord_data, 
                                                  test_post=test_post, test_o=test_o,  # datasets
                                                  traintestmatchdf=traintestmatchdf,
                                                  outcome=outcome, time_elapsed=time_elapsed, 
                                                  plot = plot,matchprobweight=matchprobweight,
                                                  time_window=time_window, interval=interval,
                                                  cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                                  d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                                  thresh_val = thresh_val, printtrace=printtrace,
                                                  userchoose=userchoose, seed=seed,
                                                  parallel=parallel, loocv=FALSE, mtype=mtype,
                                                  biasm=biasm,m=m
            ) 
            retlist <- list(pred_res = predict_test_result,
                            loocv_res =  list(loocv_test_result),
                            nearest_n=nearest_n)
            names(retlist$loocv_res) <- c(paste0("nearest_",nearest_n))
            retlist$loocv_score <- loocvperf(retlist$loocv_res, ord_data, bias=biasm, nearest_n)
            return(retlist)
        } else {
            predict_test_result <- pat_level_func(
                                                  ref=ref,nearest=nearest_n, # number to play with 
                                                  dist_fam = dist_fam, # for gamlss distribution
                                                  train_post=train_post, ord_data=ord_data, 
                                                  test_post=test_post, test_o=test_o,  # datasets
                                                  traintestmatchdf=traintestmatchdf,
                                                  outcome=outcome, time_elapsed=time_elapsed, 
                                                  plot = plot,matchprobweight=matchprobweight,
                                                  time_window=time_window, interval=interval,
                                                  cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                                  d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                                  thresh_val = thresh_val, printtrace=printtrace,
                                                  userchoose=userchoose, seed=seed,
                                                  parallel=parallel, loocv=FALSE, mtype=mtype,
                                                  biasm=biasm,m=m
            ) 
            return(list(pred_res = predict_test_result,
                        nearest_n=nearest_n))
        }
    }


}
