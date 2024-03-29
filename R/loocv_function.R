#' Leave One Out Cross Validation function that Calculates
#' Bias, Coverage and Precision for neighbors-based prediction
#' 
#' The function calculates three frequency-based parameters to 
#' demonstate statistical quality of the neighbors-based prediction.
#' The parameters (i.e. performances measures) are 
#' bias, coverage and the precision (50% prediction interval) 
#' 
#' @param nearest_n - Numeric vector indicating the nearest number of matches 
#' to select from. The \code{loocv_function} will iterate through 
#' the number within the vector and select the number of matches that has
#' optimal bias, coverage, and precision.
#' The number of nearest_n should be within the range of number of individuals
#' That are in the data. For example if there are 500 individuals in the data,
#' one could do \code{nearest_n <- 10:100}
#' @param dist_fam  Gamlss distribution specification using the \code{\link{gamlss.dist}} package. The specification for a normal distribution would be \code{gamlss.dist::NO}. For other distributions see \code{\link{gamlss.dist}}.
#' @param preproc  Preprocessed object as an ouput from \code{preproc()} function. The object contains preprocessed training and testing data.frame formatted as a list component that was generated from the \code{\link{preproc}} function
#' @param plot  Logical (\code{TRUE/FALSE}) that specifies whether to output individual precision plots
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
#' @param m - For \code{mtype = 4}, which is type 4 matching from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}, the Number of repititions for obtaining \eqn{\dot{y}} in terms of the predictive mean matching process (default = 5).
#' @param perfrank - String indicating how to rank the performance of the LOOCV. Default is `perfrank == "cov"`, which prioritizes LOOCV based on prefering coverage values that are close to `opt_cov`. Then the lowest `rmse` value then `prec` value is prefered,
#' @param opt_cov - Float value indicating optimal coverage value used for `perfrank`. Defaults to 0.5
#' @param perf_round_by - Integer value to indicate what decimal point will the performance values should be rounded by. Default is `perf_round_by = 4`, set to smaller value to be less coarse about ranking `nearest_n` values.
#' @param \dots Passed down to \code{gamlss}
#' 
#' @return Returns a list of 3 lists and a value. 1) \code{pred_res} contains a list of predicted values for the training data (\code{pred_train}) and test data (\code{pred_test}), the performance (\code{bias},\code{rmse},\code{zscore},\code{iqrcoverage},\code{precisionvec}, and number of dropped cases in fitting the gamlss model(\code{dropped_cases})); 2) \code{loocv_res} contains the same lists described above for each models fitted using different values of number of nearest neighbors; 3) \code{loocv_score} contains the summarised performance measures as a data frame; 4) \code{nearest_n} contains the optimal number of matches based on the aggregate performance metric.
#' 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# LOOCV Function ----
# - - - - - - - - - - - - - - - - - - - -#
# now using that we need to come up with the right number of matches that will give us the best bias and coverage.
loocv_function <- function(nearest_n = seq(20,150,by=10), # number to play with 
                           dist_fam = NULL, # for gamlss distribution
                           preproc,
                           plot = FALSE,
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
                           m=5,
                           perfrank="cov",
                           opt_cov = 0.5,
                           perf_round_by=4,
                           ...) {

    
    # - - - - - - - - - - - - - - - - - - - - - # 
    # FIT REFERENCE GAMLSS MODEL
    # - - - - - - - - - - - - - - - - - - - - - # 

    # Even if a time_window is specified, the entier training data will be used to estimate a model. THe time_window only specifies the prediction window that the user wants to capture. All other prediction time points will not be accounted for in the output. 
    
    ref <- invisible(
        fitrefgamlss(
            dist_fam = dist_fam, # for gamlss distribution
            train_post=preproc$train_post, 
            test_post=preproc$test_post, 
            outcome=preproc$varname[1], time_elapsed=preproc$varname[2], 
            time_window=time_window,
            cs=cs,
            dfspec=dfspec,
            d_f_m=d_f_m, ptr_m=ptr_m,
            d_f_s=d_f_s,
            d_f_n=d_f_n,
            d_f_t=d_f_t,
            ...) 
    )

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
                                         test_o=preproc$test_o, 
                                         ord_data=preproc$train_o, 
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
                                                train_post=preproc$train_post, ord_data=preproc$train_o, 
                                                test_post=preproc$test_post, test_o=preproc$test_o,  # datasets
                                                traintestmatchdf=traintestmatchdf,
                                                outcome=preproc$varname[1], time_elapsed=preproc$varname[2], idname=preproc$varname[3],
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
                perfdf <- loocv_perf(
                    loocv_test_result,
                    outcome=preproc$varname[1],
                    nearest_n=nearest_n,
                    opt_cov = opt_cov,
                    perf_round_by=perf_round_by
                )
            } else {
                
                # PERFORMANCE CALCULATION --------------------------------
                
                if(perfrank=="totscore"){
                    
                    # - - - - - - - - - - - - - - - - - - - - - - #
                    # Use a weighted scoring fucntion of bias, coverage, and precision to choose optimal m
                    # - - - - - - - - - - - - - - - - - - - - - - #
                    # opt_n_index <- loocvperf(loocv_test_result, ord_data, bias=biasm, nearest_n) %>%
                    #     dplyr::select(.data$totscore) %>% unlist(.) %>% which.min(.) %>% as.vector(.)
                    
                    perfdf <- loocv_perf(
                        loocv_test_result,
                        outcome=preproc$varname[1],
                        nearest_n=nearest_n,
                        opt_cov = opt_cov,
                        perf_round_by=perf_round_by
                    )
                    
                    opt_n <- perfdf %>%
                        arrange(perfdf$totscore)  %>%
                        head(1) %>%
                        .[,"nearest_n"] %>% .[[1]]
                    
                    
                } else if(perfrank=="cov"){
                    
                    
                    perfdf <- loocv_perf(
                        loocv_test_result,
                        outcome=preproc$varname[1],
                        nearest_n=nearest_n,
                        opt_cov = opt_cov,
                        perf_round_by=perf_round_by
                    )
                    
                    
                    opt_n <- perfdf %>%
                        arrange(.data$covdiff, .data$rmse, .data$prec)  %>%
                        head(1) %>%
                        .[,"nearest_n"]  %>% .[[1]]
                    
                    
                } else if(perfrank=="bias"){
                    
                    perfdf <- loocv_perf(
                        loocv_test_result,
                        outcome=preproc$varname[1],
                        nearest_n=nearest_n,
                        opt_cov = opt_cov,
                        perf_round_by=perf_round_by
                    )
                    
                    opt_n <- perfdf %>%
                        arrange(.data$rmse, .data$covdiff, .data$prec)  %>%
                        head(1) %>%
                        .[,"nearest_n"]  %>% .[[1]]
                    
                }
                
                # select first n if they are all the same
                # if(length(opt_n) > 1){
                #     opt_n <- opt_n[1] # select first one
                # } 
                
                opt_n_index <- which(perfdf[,"nearest_n"] == opt_n)
                
            }
            names(loocv_test_result) <- paste0("nearest_",nearest_n)
            
            print(paste0("bias: ", perfdf$rmse, "cov: ", perfdf$cov, "prec: ", perfdf$prec, " from ", names(loocv_test_result))[[opt_n_index]])
            print(paste0("Number of misses is: ",loocv_test_result[[opt_n_index]]$dropped_cases,' cases'))
            print(paste0("Distribution chosen for matched GAMLSS: ", ref$gamlss_dist))
            print(paste0("Optimal Number of Matches is: ", nearest_n[[opt_n_index]]))
        } else {
            # if loocv = FALSE
            if(is.null(userchoose)){
                stop("User must choose the nearest number of matches via specifying userchoose = N!")
            }
            nearest_n <-  userchoose
            opt_n_index <- which(nearest_n == userchoose)
        }
        # - - - - - - - - - - - - - - - - - - - - - - #
        # Run the result on test set
        # - - - - - - - - - - - - - - - - - - - - - - #
        predict_test_result <- pat_level_func(
                                              ref=ref,nearest=nearest_n[opt_n_index], # number to play with 
                                              dist_fam = dist_fam, # for gamlss distribution
                                              train_post=preproc$train_post, ord_data=preproc$train_o, 
                                              test_post=preproc$test_post, test_o=preproc$test_o,  # datasets
                                              traintestmatchdf=traintestmatchdf,
                                              outcome=preproc$varname[1], time_elapsed=preproc$varname[2], idname=preproc$varname[3],
                                              plot = plot,matchprobweight=matchprobweight,
                                              time_window=time_window, interval=interval,
                                              cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                              d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                              thresh_val = thresh_val, printtrace=printtrace,
                                              userchoose=userchoose, seed=seed,
                                              parallel=parallel, loocv=FALSE, mtype=mtype,
                                              biasm=biasm,m=m
        ) 
        # calculate performance for test data
        
        perfdf_test <- loocv_perf(
            predict_test_result,
            outcome=preproc$varname[1],
            nearest_n=nearest_n[opt_n_index],
            opt_cov = opt_cov,
            perf_round_by=perf_round_by,
            train=FALSE
        )
        
        # extract prediction results from the optimum nearest n 
        return(list(pred_res = predict_test_result,
                    test_score = perfdf_test,
                    loocv_res =  loocv_test_result,
                    loocv_score = perfdf,
                    # loocv_score = loocvperf(loocv_test_result, ord_data, bias=biasm, nearest_n),
                    nearest_n=nearest_n[opt_n_index]))
    } else {
        if(loocv){

            loocv_test_result <- pat_level_func(
                                                ref=ref,nearest=nearest_n, # number to play with 
                                                dist_fam = dist_fam, # for gamlss distribution
                                                train_post=preproc$train_post, ord_data=preproc$train_o, 
                                                test_post=preproc$test_post, test_o=preproc$test_o,  # datasets
                                                traintestmatchdf=traintestmatchdf,
                                                outcome=preproc$varname[1], time_elapsed=preproc$varname[2], idname=preproc$varname[3],
                                                plot = plot,matchprobweight=matchprobweight,
                                                time_window=time_window, interval=interval,
                                                cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                                d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                                thresh_val = thresh_val, printtrace=printtrace,
                                                userchoose=userchoose, seed=seed,
                                                parallel=parallel, loocv=loocv, mtype=mtype,
                                                biasm=biasm,m=m
            ) 
            
            perfdf <- loocv_perf(
                loocv_test_result,
                outcome=preproc$varname[1],
                nearest_n=nearest_n,
                opt_cov = opt_cov,
                perf_round_by=perf_round_by
            )
            
            predict_test_result <- pat_level_func(
                                                  ref=ref,nearest=nearest_n, # number to play with 
                                                  dist_fam = dist_fam, # for gamlss distribution
                                                  train_post=preproc$train_post, ord_data=preproc$train_o, 
                                                  test_post=preproc$test_post, test_o=preproc$test_o,  # datasets
                                                  traintestmatchdf=traintestmatchdf,
                                                  outcome=preproc$varname[1], time_elapsed=preproc$varname[2], idname=preproc$varname[3],
                                                  plot = plot,matchprobweight=matchprobweight,
                                                  time_window=time_window, interval=interval,
                                                  cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                                  d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                                  thresh_val = thresh_val, printtrace=printtrace,
                                                  userchoose=userchoose, seed=seed,
                                                  parallel=parallel, loocv=FALSE, mtype=mtype,
                                                  biasm=biasm,m=m
            ) 
            
            perfdf_test <- loocv_perf(
                predict_test_result,
                outcome=preproc$varname[1],
                nearest_n=nearest_n,
                opt_cov = opt_cov,
                perf_round_by=perf_round_by,
                train=FALSE
            )
            
            retlist <- list(pred_res = predict_test_result,
                            loocv_res = loocv_test_result,
                            test_score = perfdf_test,
                            loocv_res = loocv_test_result,
                            loocv_score = perfdf,
                            nearest_n=nearest_n)
            names(retlist$loocv_res) <- c(paste0("nearest_",nearest_n))
            return(retlist)
        } else {
            predict_test_result <- pat_level_func(
                                                  ref=ref,nearest=nearest_n, # number to play with 
                                                  dist_fam = dist_fam, # for gamlss distribution
                                                  train_post=preproc$train_post, ord_data=preproc$train_o, 
                                                  test_post=preproc$test_post, test_o=preproc$test_o,  # datasets
                                                  traintestmatchdf=traintestmatchdf,
                                                  outcome=preproc$varname[1], time_elapsed=preproc$varname[2], idname=preproc$varname[3],
                                                  plot = plot,matchprobweight=matchprobweight,
                                                  time_window=time_window, interval=interval,
                                                  cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                                                  d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                                                  thresh_val = thresh_val, printtrace=printtrace,
                                                  userchoose=userchoose, seed=seed,
                                                  parallel=parallel, loocv=FALSE, mtype=mtype,
                                                  biasm=biasm,m=m
            ) 
            
            perfdf_test <- loocv_perf(
                predict_test_result,
                outcome=preproc$varname[1],
                nearest_n=nearest_n,
                opt_cov = opt_cov,
                perf_round_by=perf_round_by,
                train=FALSE
            )
            return(list(pred_res = predict_test_result,
                        test_score = perfdf_test,
                        nearest_n=nearest_n))
        }
    }


}
