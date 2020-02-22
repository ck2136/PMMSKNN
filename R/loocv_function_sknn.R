#' Calculate bias, coverage and CI width for neighbors-based prediction
#' 
#' The function calculate three frequency-based parameters to 
#' demonstate statistical quality of the neighbors-based prediction.
#' The parameters are bias, coverage and the 50% prediction interval 
#' width.
#' @param nearest_n Numeric vector with number of matches per scenario
#' @param dist_fam  gamlss distribution specification
#' @param fulldata - dataset, the full data typically the \code{train_post} list component 
#' @param train_test    Column name indicating whether individual belongs to
#' Training or Testing dataset. Train = 1, Test = 2 by default.
#' @param patid         Column name indicating patient id.
#' @param formula       Formula indicating the variables used for matching.
#' (e.g. \code{ ~ var1 + var2 + var3 }).
#' @param train_post - datasets, typically the \code{train_post} list component 
#'  of the object produced by \code{\link{preproc}}.
#' @param test_post Idem, component \code{test_post}
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
#' @param parallel - Number of cores used for parallel computing. Default = 1
#' @param loocv - Whether or not to perform leave one out cross validation or just go straight to prediction. Should
#'  have the userchoose value specified if `loocv=FALSE`
#' @param biasm - Column indicating which bias score to use for choosing optimal n. Default
#' is \code{'raw'}. Options: \code{'raw','rmse','zsc'}.
#' @param seed - Seed for probability sampling of the nearest matches
#' @param perfrank - String indicating how to rank the performance of the LOOCV. Default is `perfrank == "cov"`, which prioritizes LOOCV based on prefering coverage values that are close to 0.5. Then the lowest `rmse` value then `prec` value is prefered,
#' @param perf_round_by - Integer value to indicate what decimal point will the performance values should be rounded by. Default is `perf_round_by = 4`, set to smaller value to be less coarse about ranking `nearest_n` values.
#' @param \dots Passed down to \code{gamlss}
#' 
#' @return There are many possible return values 
#' 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# LOOCV Function ----
# - - - - - - - - - - - - - - - - - - - -#
# now using that we need to come up with the right number of matches that will give us the best bias and coverage.
loocv_function_sknn <- function(nearest_n = seq(20,150,by=10), # number to play with 
                           dist_fam = NULL, # for gamlss distribution
                           fulldata,
                           train_test = "train_test",
                           patid = "patient_id",
                           formula,
                           train_post, 
                           test_post, 
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
                           parallel=NULL,
                           loocv=TRUE,
                           seed = 1234,
                           biasm="raw",
                           perfrank="cov",
                           perf_round_by=4,
                           ...) {

    # - - - - - - - - - - - - - - - - - - - - - # 
    # DATASET MANIPULATION
    # - - - - - - - - - - - - - - - - - - - - - # 
    train <- fulldata %>%
        filter(train_test == 1) %>% 
        distinct_(.dots = patid, .keep_all=TRUE) %>%
        dplyr::select(patid, all.vars(formula))
    test <- fulldata %>%
        filter(train_test == 2) %>% 
        distinct_(.dots = patid, .keep_all=TRUE) %>%
        dplyr::select(patid, all.vars(formula))

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


    # - - - - - - - - - - - - - - - - - - - - - # 
    # GENERATE ARRAY OF MATCHED TRAINING AND TESTING ID
    # - - - - - - - - - - - - - - - - - - - - - # 

    traintestmatchdf <- matchIdExtractsknn(
                               data = fulldata,
                               train_test = train_test,
                               patid = patid,
                               formula = formula
    )
    
    # iterate through the nearest_n list calculate bias coverage for each iteration
    # for e.g. if we go form 10:100 by 10's we will have a list with 10 results 

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # pat_level_func() is the loocv function 
    # Patient level function this is the workhorse 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if(!is.null(parallel)){
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Parallel solution
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # - - - - - - - - - - - - - - - - - - - - - - #
        # Calculation of Weighted Z-score, Coverage, and Bias to select Optimal N #
        # - - - - - - - - - - - - - - - - - - - - - - #

        if(length(nearest_n) != 1){

            if(loocv){
                loocv_test_result <- pat_level_func_sknn(
                    ref=ref,nearest=nearest_n, # number to play with 
                    dist_fam = dist_fam, # for gamlss distribution
                    patid = patid,
                    train_post=train_post, 
                    test_post=test_post, 
                    traintestmatchdf=traintestmatchdf,
                    train = train,
                    test = test,
                    outcome=outcome, time_elapsed=time_elapsed, 
                    plot = plot,matchprobweight=matchprobweight,
                    time_window=time_window, interval=interval,
                    cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                    d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                    thresh_val=thresh_val, printtrace=printtrace,
                    userchoose=userchoose, seed=seed,
                    parallel=parallel, loocv=loocv, 
                    biasm=biasm
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
                        outcome=outcome,
                        nearest_n=nearest_n,
                        perf_round_by=perf_round_by
                    )
                } else {
                    
                    
                    if(perfrank=="totscore"){
                        
                        # - - - - - - - - - - - - - - - - - - - - - - #
                        # Use a weighted scoring fucntion of bias, coverage, and precision to choose optimal m
                        # - - - - - - - - - - - - - - - - - - - - - - #
                        opt_n_index <- loocvperf(loocv_test_result, ord_data, bias=biasm, nearest_n) %>%
                            dplyr::select(.data$totscore) %>% unlist(.) %>% which.min(.) %>% as.vector(.)
                        
                        
                    } else if(perfrank=="cov"){
                        
                        
                        perfdf <- loocv_perf(
                            loocv_test_result,
                            outcome=outcome,
                            nearest_n=nearest_n,
                            perf_round_by=perf_round_by
                        )
                        
                        opt_n <- perfdf %>%
                            arrange(.data$covdiff, .data$rmse, .data$prec)  %>%
                            head(1) %>%
                            .[,"nearest_n"] 
                        
                        
                    } else if(perfrank=="bias"){
                        
                        perfdf <- loocv_perf(
                            loocv_test_result,
                            outcome=outcome,
                            nearest_n=nearest_n,
                            perf_round_by=perf_round_by
                        )
                        
                        opt_n <- perfdf %>%
                            arrange(.data$rmse, .data$covdiff, .data$prec)  %>%
                            head(1) %>%
                            .[,"nearest_n"] 
                        
                        
                    } 
                    if(length(opt_n) > 1){
                        opt_n <- opt_n[1] # select first one
                    } 
                    
                    opt_n_index <- which(perfdf[,"nearest_n"] == opt_n)
                    
                }
                names(loocv_test_result) <- paste0("nearest_",nearest_n)
                
                print(paste0(mean(loocv_test_result[[opt_n_index]]$rmse, na.rm=TRUE), " from ", names(loocv_test_result))[[opt_n_index]])
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
            predict_test_result <- pat_level_func_sknn(
                nearest=nearest_n[opt_n_index], loocv=FALSE, ref=ref,
                dist_fam = dist_fam, # for gamlss distribution
                patid = patid,
                train_post=train_post, 
                test_post=test_post, 
                traintestmatchdf=traintestmatchdf,
                train = train,
                test = test,
                outcome=outcome, time_elapsed=time_elapsed, 
                plot = plot,matchprobweight=matchprobweight,
                time_window=time_window, interval=interval,
                cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                thresh_val=thresh_val, printtrace=printtrace,
                userchoose=userchoose, seed=seed,
                parallel=parallel,
                biasm=biasm
                )
            
            perfdf_test <- loocv_perf(
                predict_test_result,
                outcome=outcome,
                nearest_n=nearest_n[opt_n_index],
                perf_round_by=perf_round_by,
                train=FALSE
            )
            
            # extract prediction results from the optimum nearest n 
            
            # extract prediction results from the optimum nearest n 
            return(list(pred_res = predict_test_result,
                        test_score = perfdf_test,
                        loocv_res =  loocv_test_result,
                        loocv_score = perfdf,
                        nearest_n=nearest_n[opt_n_index]))
            
            # return(list(pred_res = predict_test_result,
            #             loocv_res =  loocv_test_result,
            #             loocv_score = loocvperf(loocv_test_result, train %>% dplyr::select(patid) %>% rename(id = 1), bias=biasm, nearest_n),
            #             nearest_n=nearest_n[opt_n_index]))
        } else {
            if(loocv){
                loocv_test_result <- pat_level_func_sknn(
                    ref=ref,nearest=nearest_n, # number to play with 
                    dist_fam = dist_fam, # for gamlss distribution
                    patid = patid,
                    train_post=train_post, 
                    test_post=test_post, 
                    traintestmatchdf=traintestmatchdf,
                    train = train,
                    test = test,
                    outcome=outcome, time_elapsed=time_elapsed, 
                    plot = plot,matchprobweight=matchprobweight,
                    time_window=time_window, interval=interval,
                    cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                    d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                    thresh_val=thresh_val, printtrace=printtrace,
                    userchoose=userchoose, seed=seed,
                    parallel=parallel, loocv=loocv, 
                    biasm=biasm
                    )
                perfdf <- loocv_perf(
                    loocv_test_result,
                    outcome=outcome,
                    nearest_n=nearest_n,
                    perf_round_by=perf_round_by
                )
                predict_test_result <- pat_level_func_sknn(
                    ref=ref,nearest=nearest_n, # number to play with 
                    dist_fam = dist_fam, # for gamlss distribution
                    patid = patid,
                    train_post=train_post, 
                    test_post=test_post, 
                    traintestmatchdf=traintestmatchdf,
                    train = train,
                    test = test,
                    outcome=outcome, time_elapsed=time_elapsed, 
                    plot = plot,matchprobweight=matchprobweight,
                    time_window=time_window, interval=interval,
                    cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                    d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                    thresh_val=thresh_val, printtrace=printtrace,
                    userchoose=userchoose, seed=seed,
                    parallel=parallel, loocv=FALSE, 
                    biasm=biasm
                    )
                perfdf_test <- loocv_perf(
                    predict_test_result,
                    outcome=outcome,
                    nearest_n=nearest_n,
                    perf_round_by=perf_round_by,
                    train=FALSE
                )
                retlist <- list(pred_res = predict_test_result,
                                loocv_res = loocv_test_result,
                                test_score = perfdf_test,
                                loocv_res = loocv_test_result,
                                loocv_score = perfdf,
                                nearest_n=nearest_n)
                # retlist <- list(pred_res = predict_test_result,
                #                 loocv_res =  list(loocv_test_result),
                #                 nearest_n=nearest_n)
                names(retlist$loocv_res) <- c(paste0("nearest_",nearest_n))
                # retlist$loocv_score <- loocvperf(retlist$loocv_res, train %>% dplyr::select(patid) %>% rename(id = 1), bias=biasm, nearest_n)
                return(retlist)
            } else {
                predict_test_result <- pat_level_func_sknn(
                    nearest=nearest_n, ref=ref,
                    dist_fam = dist_fam, # for gamlss distribution
                    patid = patid,
                    train_post=train_post, 
                    test_post=test_post, 
                    traintestmatchdf=traintestmatchdf,
                    train = train,
                    test = test,
                    outcome=outcome, time_elapsed=time_elapsed, 
                    plot = plot,matchprobweight=matchprobweight,
                    time_window=time_window, interval=interval,
                    cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                    d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                    thresh_val=thresh_val, printtrace=printtrace,
                    userchoose=userchoose, seed=seed,
                    parallel=parallel, loocv=FALSE, 
                    biasm=biasm
                    )
                perfdf_test <- loocv_perf(
                    predict_test_result,
                    outcome=outcome,
                    nearest_n=nearest_n,
                    perf_round_by=perf_round_by,
                    train=FALSE
                )
                return(list(pred_res = predict_test_result,
                            test_score = perfdf_test,
                            nearest_n=nearest_n))
                # return(list(pred_res = predict_test_result,
                #             nearest_n=nearest_n))
            }
        }

    } else {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Non Parallel solution
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # if we're doing this on training data
        loocv_test_result <- pat_level_func_sknn(
            nearest=nearest_n, ref=ref,
            dist_fam = dist_fam, # for gamlss distribution
            patid = patid,
            train_post=train_post, 
            test_post=test_post, 
            traintestmatchdf=traintestmatchdf,
            train = train,
            test = test,
            outcome=outcome, time_elapsed=time_elapsed, 
            plot = plot,matchprobweight=matchprobweight,
            time_window=time_window, interval=interval,
            cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
            d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
            thresh_val=thresh_val, printtrace=printtrace,
            userchoose=userchoose, seed=seed,
            parallel=parallel, loocv=TRUE, 
            biasm=biasm
            )

        # - - - - - - - - - - - - - - - - - - - - - - #
        # Calculation of Weighted Z-score, Coverage, and Bias to select Optimal N #
        # - - - - - - - - - - - - - - - - - - - - - - #

        #opt <- n <- index <- which.min(as.vector(sapply(loocv <- test <- result, function(x) { mean(x$rmse, na.rm=TRUE)})))
        if(length(nearest_n) != 1){

            if(!is.null(userchoose)){
                opt_n_index <- which(nearest_n == userchoose)
            } else {
                
                if(!is.null(userchoose)){
                    if(!any(nearest_n %in% userchoose)){
                        stop("Optimal n not in the range of nearest_n specified!")
                    }
                    opt_n_index <- which(nearest_n == userchoose)
                } else {
                    
                    # PERFORMANCE CALCULATION --------------------------------
                    
                    if(perfrank=="totscore"){
                        
                        # - - - - - - - - - - - - - - - - - - - - - - #
                        # Use a weighted scoring fucntion of bias, coverage, and precision to choose optimal m
                        # - - - - - - - - - - - - - - - - - - - - - - #
                        opt_n_index <- loocvperf(loocv_test_result, ord_data, bias=biasm, nearest_n) %>%
                            dplyr::select(.data$totscore) %>% unlist(.) %>% which.min(.) %>% as.vector(.)
                        
                    } else if(perfrank=="cov"){
                        perfdf <- loocv_perf(
                            loocv_test_result,
                            outcome=outcome,
                            nearest_n=nearest_n,
                            perf_round_by=perf_round_by
                        )
                        opt_n <- perfdf %>%
                            arrange(.data$covdiff, .data$rmse, .data$prec)  %>%
                            head(1) %>%
                            .[,"nearest_n"] 
                        
                    } else if(perfrank=="bias"){
                        
                        perfdf <- loocv_perf(
                            loocv_test_result,
                            outcome=outcome,
                            nearest_n=nearest_n,
                            perf_round_by=perf_round_by
                        )
                        
                        opt_n <- perfdf %>%
                            arrange(.data$rmse, .data$covdiff, .data$prec)  %>%
                            head(1) %>%
                            .[,"nearest_n"] 
                        
                    } 
                    
                    if(length(opt_n) > 1){
                        opt_n <- opt_n[1] # select first one
                    } 
                    
                    opt_n_index <- which(perfdf[,"nearest_n"] == opt_n)
                }
            }

            print(paste0(mean(loocv_test_result[[opt_n_index]]$rmse, na.rm=TRUE), " from ", names(loocv_test_result))[[opt_n_index]])
            print(paste0("Number of misses is: ",loocv_test_result[[opt_n_index]]$dropped_cases,' cases'))
            print(paste0("Distribution chosen for matched GAMLSS: ", ref$gamlss_dist))
            print(paste0("Optimal Number of Matches is: ", nearest_n[[opt_n_index]]))
            predict_test_result <- pat_level_func_sknn(
                nearest=nearest_n[opt_n_index], ref=ref,
                dist_fam = dist_fam, # for gamlss distribution
                patid = patid,
                train_post=train_post, 
                test_post=test_post, 
                traintestmatchdf=traintestmatchdf,
                train = train,
                test = test,
                outcome=outcome, time_elapsed=time_elapsed, 
                plot = plot,matchprobweight=matchprobweight,
                time_window=time_window, interval=interval,
                cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                thresh_val=thresh_val, printtrace=printtrace,
                userchoose=userchoose, seed=seed,
                parallel=parallel, loocv=FALSE, 
                biasm=biasm
                )
            
            perfdf_test <- loocv_perf(
                predict_test_result,
                outcome=outcome,
                nearest_n=nearest_n[opt_n_index],
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

            predict_test_result <- pat_level_func_sknn(
                nearest=nearest_n, ref=ref,
                dist_fam = dist_fam, # for gamlss distribution
                patid = patid,
                train_post=train_post, 
                test_post=test_post, 
                traintestmatchdf=traintestmatchdf,
                train = train,
                test = test,
                outcome=outcome, time_elapsed=time_elapsed, 
                plot = plot,matchprobweight=matchprobweight,
                time_window=time_window, interval=interval,
                cs=cs, dfspec=dfspec, d_f_m=d_f_m, ptr_m=ptr_m,
                d_f_s=d_f_s, d_f_n=d_f_n, d_f_t=d_f_t,
                thresh_val=thresh_val, printtrace=printtrace,
                userchoose=userchoose, seed=seed,
                parallel=parallel, loocv=FALSE, 
                biasm=biasm
                )
            
            perfdf_test <- loocv_perf(
                predict_test_result,
                outcome=outcome,
                nearest_n=nearest_n,
                perf_round_by=perf_round_by,
                train=FALSE
            )
            return(list(pred_res = predict_test_result,
                        test_score = perfdf_test,
                        nearest_n=nearest_n))
        
        }

    }

}
