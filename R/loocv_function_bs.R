#' Leave One Out Cross Validation function that Calculates
#' Bias, Coverage and Precision for neighbors-based prediction
#' Using a simplified BrokenStick prediction framework
#' 
#' The function calculates three frequency-based parameters to 
#' demonstate statistical quality of the neighbors-based prediction.
#' The parameters (i.e. performances measures) are 
#' bias, coverage and the precision (50% prediction interval) 
#' 
#' @param nearest_n  Numeric vector indicating the nearest number of matches 
#' to select from. The \code{loocv_function} will iterate through 
#' the number within the vector and select the number of matches that has
#' optimal bias, coverage, and precision.
#' The number of nearest_n should be within the range of number of individuals
#' That are in the data. For example if there are 500 individuals in the data,
#' one could do \code{nearest_n <- 10:100}
#' @param train_post  Data frame that contains the post-baseline observations from the training dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param ord_data Data frame. Specifically, training data with patient_id ordered based on fitted distal outcome value using predicted mean matching.
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$train_o} would be used for this parameter.
#' @param test_post Data frame that contains the post-baseline observations from the testing dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param test_o  Data frame. Specifically, testing data with patient_id ordered based on fitted distal outcome value using predicted mean matching.
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$test_o} would be used for this parameter.
#' @param outcome     Name of the outcomes variable (type=string)
#' @param plot  Logical (\code{TRUE/FALSE}) that specifies whether to output individual precision plots
#' @param matchprobweight  Logical (\code{TRUE/FALSE}) that specifies whether to utilize probability sampling
#'  when doing the mean matching. If TRUE, matches nearest n weighted on differnce in 
#'  predicted outcome.
#' @param time_window  vector of numbers for `centiles.pred()`, `xvalues` argument. For example, specify such as \code{c(10:30)}
#' @param interval  Int value specifying the interval of individuals to skip 
#' @param thresh_val  Numeric value indicating value of bias to ignore (not include in output) in terms of the leave-one-out cross validation process. The default is set to \code{thresh_val = 10000}
#' @param userchoose  Int value indicating the choice that the user wants to use for the number of nearest matches 
#' @param seed  Seed for probability sampling of the nearest matches
#' @param parallel  Number of cores used for the leave-one-out cross validation process. Default = 1
#' @param loocv  Logical (\code{TRUE/FALSE}) that specifies whether 
#' or not to perform leave-one-out cross validation or just output 
#' predictions without hyperparameter tuning. If \code{loocv=FALSE}, then
#' users need to specify the value of the nearest_n 
#' @param mtype  Integer value indicating matching type. Default is set to 1 which follows the
#' matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. \code{mtype} values are 
#' from \code{0} to \code{4}
#' @param biasm  Column indicating which bias measure to use for 
#' choosing the optimal nearest number of neighbors. 
#' Default is \code{'raw'}. Options: \code{'raw','rmse','zsc'}.
#' @param m - For \code{mtype = 4}, which is type 4 matching from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}, the Number of repititions for obtaining \eqn{\dot{y}} in terms of the predictive mean matching process.
#' @param perfrank  String indicating how to rank the performance of the LOOCV. Default is `perfrank == "cov"`, which prioritizes LOOCV based on prefering coverage values that are close to 0.5. Then the lowest `rmse` value then `prec` value is prefered,
#' @param opt_cov   Numeric value to indicate what the optimal coverage value is for calculating the performance. Default is `opt_cov = 0.5` (i.e. 50% coverage).
#' @param perf_round_by  Integer value to indicate what decimal point will the performance values should be rounded by. Default is `perf_round_by = 4`, set to smaller value to be less coarse about ranking `nearest_n` values.
#' @param \dots Passed down to \code{gamlss}
#' 
#' @return Returns a list of 3 lists and a value. 1) \code{pred_res} contains a list of predicted values for the training data (\code{pred_train}) and test data (\code{pred_test}), the performance (\code{bias},\code{rmse},\code{zscore},\code{iqrcoverage},\code{precisionvec}, and number of dropped cases in fitting the gamlss model(\code{dropped_cases})); 2) \code{loocv_res} contains the same lists described above for each models fitted using different values of number of nearest neighbors; 3) \code{loocv_score} contains the summarised performance measures as a data frame; 4) \code{nearest_n} contains the optimal number of matches based on the aggregate performance metric.
#' 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# LOOCV Function ----
# - - - - - - - - - - - - - - - - - - - -#
# now using that we need to come up with the right number of matches that will give us the best bias and coverage.
loocv_function_bs <- function(nearest_n = seq(20,150,by=10), # number to play with 
                           train_post, 
                           ord_data, 
                           test_post, 
                           test_o,  # datasets
                           bs_obj, # broken stick object
                           outcome, plot = FALSE,
                           matchprobweight=FALSE,
                           time_window=NULL,
                           interval=NULL,
                           thresh_val = 10000,
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

    # NEAREST NEIGHBOR MATCHING: Description ----------------------
    #   Outerloop will iterate through a list of nearest_n numbers
    #   Innerloop will iterate through all patients and store the result of the whole iteration into an array
    #   nearest neighbor should be an even number for at least those that aren't edge cases
    #   iterate through the nearest_n list calculate bias coverage for each iteration
    #   for e.g. if we go form 10:100 by 10's we will have a list with 10 results 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Match training and testing patients based on minimum difference in the predicted outcome at out_time 
    # Need this to get test data predictions. Default matching is type 2 matching
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # PMMSKNN Prediction Code -------------------------------------------------
    message("Initiating PMMSKNN Prediction Process")
    traintestmatchdf <- matchTestDataGen(
                                         test_o=test_o, 
                                         ord_data=ord_data, 
                                         mtype=mtype
    )
    # Specify if we'll be skipping on predicting some patients in the training set
    if(is.null(interval)){
        patlistint <- seq(1,length(ord_data$id))
    } else {
        patlistint <- seq(1,length(ord_data$id), by=interval)
    }
    
    if(!is.null(parallel)) {
        # Parallel -----------------------
        message("Initiating Parallel Process")
        plan(multiprocess, workers = parallel)
        
        # If Lenght of the nearest_n checking is more than 1 then go ahead and do the loocv probably
        if(length(nearest_n) != 1){
            
            # LOOCV on Training --------------------------
            # If we are going to do loocv only run the CV if loocv True
            if(loocv) {
                
                # generate for each of the nearest n combinations
                loocv_test_result <- lapply(
                    
                    nearest_n, function(n) {
                        
                        # generate bias values 
                        predreslist <- future_lapply(patlistint , function(nthind) {
                            
                            # identify current training id 
                            curid <- ord_data$id[nthind]
                            
                            # identify matches for the curid
                            matches <- matchIdExtract(ord_data, mtype, loocv, n, m, nthind) 
                            # extract the time  component of the curid participants
                            time_curid <- train_post[train_post$patient_id %in% curid, "time"] %>% unlist %>% as.vector
                            
                            # using the match ids get predicted values of matches
                            predvaldf <- predict(bs_obj, x=time_curid, ids=matches) %>% filter(x %in% time_curid) %>% group_by(x) %>% 
                                summarise(
                                    y_hat_avg = mean(.data$yhat), 
                                    y_sd = sd(.data$yhat)
                                ) 
                            
                            # get obs value for curid
                            # calculate confidence intervals using t dist
                            data.frame(
                                patient_id = train_post %>% filter(.data$patient_id %in% curid) %>% dplyr::select(patient_id) %>% unlist %>% as.vector,
                                time = train_post %>% filter(.data$patient_id %in% curid) %>% dplyr::select(time) %>% unlist %>% as.vector,
                                obsvals =train_post %>% filter(.data$patient_id %in% curid) %>% dplyr::select(outcome) %>% unlist %>% as.vector,
                                predvals = predvaldf$y_hat_avg,
                                lower95 = predvaldf$y_hat_avg - (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches)),
                                upper95 = predvaldf$y_hat_avg + (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches))
                            )
                            
                        })
                        message("Finished LOOCV: compiling results for ", n ," matches")
                        # rbind all
                        predresdf <- data.table::rbindlist(predreslist)
                        # generate rmse, cov, and prec 
                        predresdf <- predresdf %>%
                            mutate(
                                cov = ifelse(.data$obsvals > lower95 & .data$obsvals < upper95, 1, 0),
                                prec = upper95 - lower95,
                                nearest_n = n
                            )
                        
                        # generate summarized df
                        list(
                            perfdf = data.frame(
                                rmse = sqrt(sum((predresdf$obsvals - predresdf$predvals)^2)/length(predresdf$obsvals)),
                                cov = sum(predresdf$cov) / nrow(predresdf),
                                prec = mean(predresdf$prec),
                                nearest_n = predresdf$nearest_n[1]
                            ),
                            resdf = predresdf
                        )
                    }
                )
                
                message("Finished LOOCV: compiling results for ", max(nearest_n) ," matches")
                
                perfdf = data.table::rbindlist(
                    lapply(loocv_test_result, function(x) {
                        x[[1]] # extracting the 1st list which is the performance dataframe
                    })
                )
                
                perfdf <- perfdf %>%
                    mutate(
                        rsc = (.data$rmse - min(.data$rmse, na.rm=TRUE)) / (max(.data$rmse, na.rm=TRUE) - min(.data$rmse, na.rm=TRUE)),
                        presc = (.data$prec - min(.data$prec, na.rm=TRUE)) / (max(.data$prec, na.rm=TRUE) - min(.data$prec, na.rm=TRUE))
                    ) 
                
                # Coverage calculation
                ifelse((max(abs(.data$cov - opt_cov), na.rm=TRUE) - min(abs(.data$cov - opt_cov), na.rm=TRUE)) == 0, 
                       {.data$covsc = 0},
                       {perfdf <- perfdf %>%
                           mutate(
                               covsc = (abs(.data$cov - opt_cov) - min(abs(.data$cov - opt_cov), na.rm =TRUE)) / (max(abs(.data$cov - opt_cov), na.rm=TRUE) - min(abs(.data$cov - opt_cov), na.rm=TRUE))
                           )
                       })
                
                perfdf <- perfdf %>%
                    mutate(totscore = .data$rsc + .data$covsc + .data$presc)  
                
                if(perfrank=="totscore"){
                    
                    # Select optimal n 
                    if(nrow(perfdf %>% top_n(-1, totscore)) > 1){
                        opt_n <- perfdf %>%
                            arrange(totscore)  %>%
                            head(1) %>%
                            .[,"nearest_n"]
                    } else {
                        opt_n <- perfdf %>%
                            arrange(totscore)  %>%
                            head(1) %>%
                            .[,"nearest_n"]
                    }
                    
                } else if(perfrank=="cov"){
                    opt_n <- perfdf %>%
                        mutate(
                            covdiff = round(abs(cov - opt_cov), perf_round_by)
                        ) %>%
                        arrange(covdiff, rmse, prec)  %>%
                        head(1) %>%
                        .[,"nearest_n"]
                } else if(perfrank=="bias"){
                    opt_n <- perfdf %>%
                        mutate(
                            covdiff = round(abs(cov - opt_cov), perf_round_by),
                            rmse = round(rmse, perf_round_by)
                        ) %>%
                        arrange(rmse, covdiff, prec)  %>%
                        head(1) %>%
                        .[,"nearest_n"]
                }
                
                
            } else {
                # case when loocv is not going to run....
                # then we need to set userchoose the optimal number of matches
                if(is.null(userchoose)){
                    stop("Default number of matches not chosen")
                } else {
                    opt_n <- userchoose
                }
            } 
            
            # Run the Test Prediction Now ----------------
            
            message("Compiling results for test data")
            predict_test_result <- future_lapply(1:nrow(ord_data), function(nthind) {
                
                # identify current training id 
                curid <- ord_data$id[nthind]
                
                # extract the multiple potential test ids
                matched_test_ids <- traintestmatchdf %>% 
                    filter(train_id == curid) %>% 
                    dplyr::select(test_id) %>% 
                    unlist %>% as.vector 
                
                # identify training matches for the curid based on opt_n chosen above
                matches <- matchIdExtract(ord_data, mtype, loocv = FALSE, opt_n, m, nthind) 
                
                # time_curid <- train_post[train_post$patient_id %in% curid, "time"] %>% unlist %>% as.vector
                nest_pred <- test_post %>%
                    filter(.data$patient_id %in% matched_test_ids) %>%
                    tidyr::nest(-patient_id) %>%
                    mutate(
                        times = lapply(data, function(data) {
                            data %>% 
                                dplyr::select(time) %>% 
                                unlist %>% as.vector}),
                        predval = lapply(times, function(times) {
                            predict(bs_obj, x=times, ids=matches) %>% 
                                filter(x %in% times) %>% 
                                group_by(x) %>% 
                                summarise(
                                    y_hat_avg = mean(.data$yhat), 
                                    y_sd = sd(.data$yhat)
                                ) 
                        })
                    )
                
                predvaldf <- nest_pred %>%
                    tidyr::unnest(predval) 
                
                # calculate confidence intervals using t dist
                # get obs value for curid
                data.frame(
                    patient_id = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(patient_id) %>% unlist %>% as.vector,
                    time = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(time) %>% unlist %>% as.vector,
                    obsvals = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(outcome) %>% unlist %>% as.vector,
                    predvals = predvaldf$y_hat_avg,
                    lower95 = predvaldf$y_hat_avg - (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches)),
                    upper95 = predvaldf$y_hat_avg + (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches))
                )
                
            })
            message("Finished compiling results for test data")
            
            predict_test_result =predict_test_result[lapply(predict_test_result,nrow)>0] ## you can use sapply,rapply
            # rbind all
            predresdf <- data.table::rbindlist(predict_test_result)
            # generate rmse, cov, and prec 
            predresdf <- predresdf %>%
                mutate(
                    cov = ifelse(.data$obsvals > lower95 & .data$obsvals < upper95, 1, 0),
                    prec = upper95 - lower95
                )
            
            # Get test performance
            
            perfdf_test = data.frame(
                rmse = sqrt(sum((predresdf$obsvals - predresdf$predvals)^2)/length(predresdf$obsvals)),
                cov = sum(predresdf$cov) / nrow(predresdf),
                prec = mean(predresdf$prec)
            )
            
            message("Finished PMM Prediction Process")
            # extract prediction results from the optimum nearest n  --------------
            
            if(loocv){
                return(
                    list(
                        pred_res = predresdf,
                        test_score = perfdf_test,
                        loocv_res =  loocv_test_result,
                        loocv_score = perfdf,
                        nearest_n=opt_n)
                )
            } else {
                return(
                    list(
                        pred_res = predresdf,
                        test_score = perfdf_test,
                        nearest_n=opt_n)
                )
            }
            
        } else {
            
            # If LOOCV isn't necessary, then just go straight to testing data
            if(is.null(userchoose)){
                stop("Default number of matches not chosen")
            } else {
                opt_n <- userchoose
                
            }
            
            message("Compiling results for test data")
            predict_test_result <- lapply(1:nrow(ord_data), function(nthind) {
                
                # identify current training id 
                curid <- ord_data$id[nthind]
                
                # extract the multiple potential test ids
                matched_test_ids <- traintestmatchdf %>% filter(train_id == curid) %>% dplyr::select(test_id) %>% unlist %>% as.vector 
                
                # identify training matches for the curid based on opt_n chosen above
                matches <- matchIdExtract(ord_data, mtype, loocv=FALSE, opt_n, m, nthind) 
                
                # time_curid <- train_post[train_post$patient_id %in% curid, "time"] %>% unlist %>% as.vector
                nest_pred <- test_post %>%
                    filter(.data$patient_id %in% matched_test_ids) %>%
                    tidyr::nest(-patient_id) %>%
                    mutate(
                        time_vals = lapply(data, function(data) {
                            data %>% 
                                dplyr::select(time) %>% 
                                unlist %>% as.vector}),
                        predval = lapply(time_vals, function(time_vals) {
                            predict(bs_obj, x=time_vals, ids=matches) %>% 
                                filter(x %in% time_vals) %>% 
                                group_by(x) %>% 
                                summarise(
                                    y_hat_avg = mean(.data$yhat), 
                                    y_sd = sd(.data$yhat)
                                ) 
                        })
                    )
                
                predvaldf <- nest_pred %>%
                    tidyr::unnest(predval) 
                
                # calculate confidence intervals using t dist
                # get obs value for curid
                data.frame(
                    patient_id = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(patient_id) %>% unlist %>% as.vector,
                    time = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(time) %>% unlist %>% as.vector,
                    obsvals = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(outcome) %>% unlist %>% as.vector,
                    predvals = predvaldf$y_hat_avg,
                    lower95 = predvaldf$y_hat_avg - (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches)),
                    upper95 = predvaldf$y_hat_avg + (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches))
                )
                
            })
            message("Finished compiling results for test data")
            
            # rbind all
            predresdf <- data.table::rbindlist(predict_test_result)
            # generate rmse, cov, and prec 
            predresdf <- predresdf %>%
                mutate(
                    cov = ifelse(.data$obsvals > lower95 & .data$obsvals < upper95, 1, 0),
                    prec = upper95 - lower95,
                    nearest_n = opt_n
                    # nearest_n = n
                )
            
            perfdf_test = data.frame(
                rmse = sqrt(sum((predresdf$obsvals - predresdf$predvals)^2)/length(predresdf$obsvals)),
                cov = sum(predresdf$cov) / nrow(predresdf),
                prec = mean(predresdf$prec)
            )
            
            # extract prediction results from the optimum nearest n  --------------
            message("Finished PMM Prediction Process Only For Test Case (i.e. No LOOCV)")
            return(
                list(
                    pred_res = predict_test_result,
                    test_score = perfdf_test,
                    predresdf = predresdf,
                    nearest_n=opt_n)
            )
        }
        stopCluster(cl)
    } else {
        # Non Parallel -----------------------
        message("Initiating Non Parallel Process")
        # If Lenght of the nearest_n checking is more than 1 then go ahead and do the loocv probably
        if(length(nearest_n) != 1){
            
            # LOOCV on Training --------------------------
            # If we are going to do loocv only run the CV if loocv True
            if(loocv) {
                
                # generate for each of the nearest n combinations
                loocv_test_result <- lapply(
                    
                    nearest_n, function(n) {
                        
                        # generate bias values 
                        predreslist <- lapply(patlistint , function(nthind) {
                            
                            # identify current training id 
                            curid <- ord_data$id[nthind]
                            
                            # identify matches for the curid
                            matches <- matchIdExtract(ord_data, mtype, loocv, n, m, nthind) 
                            # extract the time  component of the curid participants
                            time_curid <- train_post[train_post$patient_id %in% curid, "time"] %>% unlist %>% as.vector
                            
                            # using the match ids get predicted values of matches
                            predvaldf <- predict(bs_obj, x=time_curid, ids=matches) %>% filter(x %in% time_curid) %>% group_by(x) %>% 
                                summarise(
                                    y_hat_avg = mean(.data$yhat), 
                                    y_sd = sd(.data$yhat)
                                ) 
                            
                            # get obs value for curid
                            # calculate confidence intervals using t dist
                            
                            data.frame(
                                patient_id = train_post %>% filter(.data$patient_id %in% curid) %>% dplyr::select(patient_id) %>% unlist %>% as.vector,
                                time = train_post %>% filter(.data$patient_id %in% curid) %>% dplyr::select(time) %>% unlist %>% as.vector,
                                obsvals = train_post %>% filter(.data$patient_id %in% curid) %>% dplyr::select_(outcome) %>% unlist %>% as.vector,
                                predvals = predvaldf$y_hat_avg,
                                lower95 = predvaldf$y_hat_avg - (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches)),
                                upper95 = predvaldf$y_hat_avg + (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches))
                            )
                            
                        })
                        message("Finished compiling results for ", n ," matches")
                        # rbind all
                        predresdf <- data.table::rbindlist(predreslist)
                        # generate rmse, cov, and prec 
                        predresdf <- predresdf %>%
                            mutate(
                                cov = ifelse(.data$obsvals > lower95 & .data$obsvals < upper95, 1, 0),
                                prec = upper95 - lower95,
                                nearest_n = n
                            )
                        
                        # generate summarized df
                        list(
                            perfdf = data.frame(
                                rmse = sqrt(sum((predresdf$obsvals - predresdf$predvals)^2)/length(predresdf$obsvals)),
                                cov = sum(predresdf$cov) / nrow(predresdf),
                                prec = mean(predresdf$prec),
                                nearest_n = predresdf$nearest_n[1]
                            ),
                            resdf = predresdf
                        )
                    }
                )
                
                message("Finished compiling results for ", max(nearest_n) ," matches")
                
                perfdf = data.table::rbindlist(
                    lapply(loocv_test_result, function(x) {
                        x[[1]] # extracting the 1st list which is the performance dataframe
                    })
                )
                
                perfdf <- perfdf %>%
                    mutate(
                        rsc = (.data$rmse - min(.data$rmse, na.rm=TRUE)) / (max(.data$rmse, na.rm=TRUE) - min(.data$rmse, na.rm=TRUE)),
                        presc = (.data$prec - min(.data$prec, na.rm=TRUE)) / (max(.data$prec, na.rm=TRUE) - min(.data$prec, na.rm=TRUE))
                        #zsc = (.data$zscore - min(.data$zscore, na.rm=TRUE)) / (max(.data$zscore, na.rm=TRUE) - min(.data$zscore, na.rm=TRUE))
                    ) 
                
                # Coverage calculation
                ifelse((max(abs(.data$cov - opt_cov), na.rm=TRUE) - min(abs(.data$cov - opt_cov), na.rm=TRUE)) == 0, 
                       {.data$covsc = 0},
                       {perfdf <- perfdf %>%
                           mutate(
                               covsc = (abs(.data$cov - opt_cov) - min(abs(.data$cov - opt_cov), na.rm =TRUE)) / (max(abs(.data$cov - opt_cov), na.rm=TRUE) - min(abs(.data$cov - opt_cov), na.rm=TRUE))
                           )
                       })
                
                perfdf <- perfdf %>%
                    mutate(totscore = .data$rsc + .data$covsc + .data$presc)  
                
                if(perfrank=="totscore"){
                    
                    # Select optimal n 
                    if(nrow(perfdf %>% top_n(-1, totscore)) > 1){
                        opt_n <- perfdf %>%
                            arrange(totscore)  %>%
                            head(1) %>%
                            .[,"nearest_n"]
                    } else {
                        opt_n <- perfdf %>%
                            arrange(totscore)  %>%
                            head(1) %>%
                            .[,"nearest_n"]
                    }
                    
                } else if(perfrank=="cov"){
                    opt_n <- perfdf %>%
                        mutate(
                            covdiff = round(abs(cov - opt_cov), perf_round_by)
                        ) %>%
                        arrange(covdiff, rmse, prec)  %>%
                        head(1) %>%
                        .[,"nearest_n"]
                } else {
                    # assume bias as default
                    opt_n <- perfdf %>%
                        mutate(
                            covdiff = round(abs(cov - opt_cov), perf_round_by)
                        ) %>%
                        arrange(rmse, covdiff, prec)  %>%
                        head(1) %>%
                        .[,"nearest_n"]
                }
            } else {
                # case when loocv is not going to run....
                # then we need to set userchoose the optimal number of matches
                if(is.null(userchoose)){
                    stop("Default number of matches not chosen")
                } else {
                    opt_n <- userchoose
                }
            } 
            
            # Run the Test Prediction Now ----------------
            
            message("Compiling results for test data")
            predict_test_result <- lapply(1:nrow(ord_data), function(nthind) {
                
                # identify current training id 
                curid <- ord_data$id[nthind]
                
                # extract the multiple potential test ids
                matched_test_ids <- traintestmatchdf %>% filter(train_id == curid) %>% dplyr::select(test_id) %>% unlist %>% as.vector 
                
                # identify training matches for the curid based on opt_n chosen above
                matches <- matchIdExtract(ord_data, mtype, loocv=FALSE, opt_n, m, nthind) 
                
                # time_curid <- train_post[train_post$patient_id %in% curid, "time"] %>% unlist %>% as.vector
                nest_pred <- test_post %>%
                    filter(.data$patient_id %in% matched_test_ids) %>%
                    tidyr::nest(-patient_id) %>%
                    mutate(
                        times = lapply(data, function(data) {
                            data %>% 
                                dplyr::select(time) %>% 
                                unlist %>% as.vector}),
                        predval = lapply(times, function(times) {
                            predict(bs_obj, x=times, ids=matches) %>% 
                                filter(x %in% times) %>% 
                                group_by(x) %>% 
                                summarise(
                                    y_hat_avg = mean(.data$yhat), 
                                    y_sd = sd(.data$yhat)
                                ) 
                        })
                    )
                
                predvaldf <- nest_pred %>%
                    tidyr::unnest(predval) 
                
                # calculate confidence intervals using t dist
                # get obs value for curid
                
                data.frame(
                    patient_id = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(patient_id) %>% unlist %>% as.vector,
                    time = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(time) %>% unlist %>% as.vector,
                    obsvals = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(outcome) %>% unlist %>% as.vector,
                    predvals = predvaldf$y_hat_avg,
                    lower95 = predvaldf$y_hat_avg - (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches)),
                    upper95 = predvaldf$y_hat_avg + (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches))
                )
                
            })
            message("Finished compiling results for test data")
            
            predict_test_result =predict_test_result[lapply(predict_test_result,nrow)>0] ## you can use sapply,rapply
            # rbind all
            predresdf <- data.table::rbindlist(predict_test_result)
            # generate rmse, cov, and prec 
            predresdf <- predresdf %>%
                mutate(
                    cov = ifelse(.data$obsvals > lower95 & .data$obsvals < upper95, 1, 0),
                    prec = upper95 - lower95
                )
            
            # Get test performance
            
            perfdf_test = data.frame(
                rmse = sqrt(sum((predresdf$obsvals - predresdf$predvals)^2)/length(predresdf$obsvals)),
                cov = sum(predresdf$cov) / nrow(predresdf),
                prec = mean(predresdf$prec)
            )
            
            message("Finished PMM Prediction Process")
            # extract prediction results from the optimum nearest n  --------------
            return(
                list(
                    pred_res = predresdf,
                    test_score = perfdf_test,
                    loocv_res =  loocv_test_result,
                    loocv_score = perfdf,
                    nearest_n=opt_n)
            )
            
        } else {
            
            # If LOOCV isn't necessary, then just go straight to testing data
            if(is.null(userchoose)){
                stop("Default number of matches not chosen")
            } else {
                opt_n <- userchoose
                
            }
            
            message("Compiling results for test data")
            predict_test_result <- lapply(1:nrow(ord_data), function(nthind) {
                
                # identify current training id 
                curid <- ord_data$id[nthind]
                
                # extract the multiple potential test ids
                matched_test_ids <- traintestmatchdf %>% filter(train_id == curid) %>% dplyr::select(test_id) %>% unlist %>% as.vector 
                
                # identify training matches for the curid based on opt_n chosen above
                matches <- matchIdExtract(ord_data, mtype, loocv = FALSE, opt_n, m, nthind) 
                
                # time_curid <- train_post[train_post$patient_id %in% curid, "time"] %>% unlist %>% as.vector
                nest_pred <- test_post %>%
                    filter(.data$patient_id %in% matched_test_ids) %>%
                    tidyr::nest(-patient_id) %>%
                    mutate(
                        time_vals = lapply(data, function(data) {
                            data %>% 
                                dplyr::select(time) %>% 
                                unlist %>% as.vector}),
                        predval = lapply(time_vals, function(time_vals) {
                            predict(bs_obj, x=time_vals, ids=matches) %>% 
                                filter(x %in% time_vals) %>% 
                                group_by(x) %>% 
                                summarise(
                                    y_hat_avg = mean(.data$yhat), 
                                    y_sd = sd(.data$yhat)
                                ) 
                        })
                    )
                
                predvaldf <- nest_pred %>%
                    tidyr::unnest(predval) 
                
                # calculate confidence intervals using t dist
                # get obs value for curid
                
                data.frame(
                    patiet_id  = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(patient_id) %>% unlist %>% as.vector,
                    time = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(time) %>% unlist %>% as.vector,
                    obsvals = test_post %>% filter(.data$patient_id %in% matched_test_ids) %>% dplyr::select(outcome) %>% unlist %>% as.vector,
                    predvals = predvaldf$y_hat_avg,
                    lower95 = predvaldf$y_hat_avg - (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches)),
                    upper95 = predvaldf$y_hat_avg + (qt(0.975, length(matches) - 1) * predvaldf$y_sd)/sqrt(length(matches))
                )
                
            })
            message("Finished compiling results for test data")
            
            # rbind all
            predresdf <- data.table::rbindlist(predict_test_result)
            # generate rmse, cov, and prec 
            predresdf <- predresdf %>%
                mutate(
                    cov = ifelse(.data$obsvals > lower95 & .data$obsvals < upper95, 1, 0),
                    prec = upper95 - lower95,
                    nearest_n = opt_n
                )
            
            # extract prediction results from the optimum nearest n  --------------
            message("Finished PMM Prediction Process")
            perfdf_test = data.frame(
                rmse = sqrt(sum((predresdf$obsvals - predresdf$predvals)^2)/length(predresdf$obsvals)),
                cov = sum(predresdf$cov) / nrow(predresdf),
                prec = mean(predresdf$prec)
            )
            return(
                list(
                    pred_res = predict_test_result,
                    test_score = perfdf_test,
                    predresdf = predresdf,
                    nearest_n=opt_n)
            )
        }
        
        
    }
}
