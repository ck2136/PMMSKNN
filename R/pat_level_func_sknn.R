#' Calculate bias, coverage and CI width for neighbors-based prediction for SKNN algorithm
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
#' The number of nearest should be within the range of number of individuals
#' That are in the data. For example if there are 500 individuals in the data,
#' one could do \code{nearest <- 10:100}
#' @param dist_fam  - gamlss distribution specification using the \code{\link{gamlss.dist}} package. The specification for a normal distribution would be \code{gamlss.dist::NO}. For other distributions see \code{\link{gamlss.dist}}.
#' @param train_post - Data frame that contains the post-baseline observations from the training dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param test_post Data frame that contains the post-baseline observations from the testing dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param test  Data frame. Specifically, full testing data 
#' @param train  Data frame. Specifically, full training data 
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$test} would be used for this parameter.
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
#' @param patid         Column name indicating patient id.
#' @param loocv - Logical (\code{TRUE/FALSE}) that specifies whether 
#' or not to perform leave-one-out cross validation or just output 
#' predictions without hyperparameter tuning. If \code{loocv=FALSE}, then
#' users need to specify the value of the nearest 
#' @param biasm - Column indicating which bias measure to use for 
#' choosing the optimal nearest number of neighbors. 
#' Default is \code{'raw'}. Options: \code{'raw','rmse','zsc'}.
#' 
#' @return Returns a list of 3 lists and a value. 1) \code{pred_res} contains a list of predicted values for the training data (\code{pred_train}) and test data (\code{pred_test}), the performance (\code{bias},\code{rmse},\code{zscore},\code{iqrcoverage},\code{precisionvec}, and number of dropped cases in fitting the gamlss model(\code{dropped_cases})); 2) \code{loocv_res} contains the same lists described above for each models fitted using different values of number of nearest neighbors; 3) \code{loocv_score} contains the summarised performance measures as a data frame; 4) \code{nearest} contains the optimal number of matches based on the aggregate performance metric.
#' 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# LOOCV Function ----
# - - - - - - - - - - - - - - - - - - - -#
pat_level_func_sknn <- function(
    nearest, 
    loocv=loocv, 
    parallel=parallel,
    train_post=train_post,
    test_post=test_post,
    train = train,
    ref=ref,
    patid = patid,
    dist_fam = dist_fam, # for gamlss distribution
    test=test,  # datasets
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
    biasm=biasm
){
    
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
        # Loop based on the nearest specified as user input 
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        loocvres <- foreach(n=nearest,
                            .export=c("train_post", "test_post","spl","ptr_m", "d_f_m", "d_f_s", "train","test", "gamlss_dist",
                                      "d_f_n", "d_f_t", "dist_fam", "spls","spln","splt",
                                      "outcome","time_elapsed", "plot", "matchprobweight","time_window",
                                      "interval","thresh_val","printtrace","userchoose", "ref",
                                      "loocv","mint","maxt", "traintestmatchdf"),                                                                                                              .packages=c("dplyr","gamlss"),
                            #.combine=list,
                            .errorhandling = "pass"
                            #.multicombine = TRUE
        ) %dopar% {
            #for(n in nearest){
            cnt = 0
            misses = 0
            
            dfList <- list()                #-- list to store training data's predicted C50 from gamlss model
            dfList_test <- list()           #-- list to store testing data's predicted C50 from gamlss model
            # centilepred <- list()           #-- store all centile for later test set merging
            # biasvec<-vector()               #-- store mean of bias
            # rmsevec <- vector()             #-- store rmse vector
            # coveragevec<-vector()           #-- store coverage vector
            # coveragevec95a<-vector()        #-- store mean of the n coverage in vector
            # #coveragevec95<-vector()
            # iqrvec<-vector()
            # ninefiveconf <- data.frame()
            # precisionvec <- list()
            crazymatch <- list()
            zsc_list <- list()
            
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # iterate through all patients and match patients according to order. ord_data is the training data
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if(is.null(interval)){
                patlistint <- seq(1,nrow(train))
            } else {
                patlistint <- seq(1,nrow(train), by=interval)
            }
            for (i in patlistint) {
                cnt = cnt + 1
                
                # - - - - - - - - - - - - - - - - - - - - - - #
                # Matching cases using absolute value of difference between Fitted values
                # - - - - - - - - - - - - - - - - - - - - - - #
                
                matchmodel <- train_post %>%
                    filter(.data$patient_id %in% (traintestmatchdf$nnarraytrain[,i] %>%
                                                      head(n=n)))
                
                # - - - - - - - - - - - - - - - - - - - - - - #
                # Fit GAMLSS model to nearest n matches
                # - - - - - - - - - - - - - - - - - - - - - - #
                
                plm <- function() {
                    
                    out <- tryCatch(
                        {
                            #message("TRY PART")
                            update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE))
                            
                        },
                        error = function(cond)
                        {
                            message(paste("ERROR in GAMLSS"))
                            #message("ORIGINAL ERROR:")
                            message(cond)
                            return(NA)
                        },
                        warning=function(cond)
                        {
                            #message(paste("GAMLSS WARNING"))
                            #message("ORIIGNAL WARNING:")
                            message(cond)
                            return(update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE)))
                            #return(NULL)
                        },
                        finally={
                            #message("PROCESSING GAMLSS")
                            #message("GAMLSS EOL SUCCESS")
                        }
                    )
                    return(out)
                }
                plmr <- plm()
                
                # if something wrong with iteration... then just don't do prediction it eats up all the memory anyways
                if(typeof(plmr) != 'list') {
                    
                    misses = misses + 1 # increment number of misses
                    message('Something wrong with plm model. Prediction not included')
                    
                } else {
                    
                    message(paste0("Getting predictions: "))
                    
                    # time window again
                    if(is.null(time_window)){
                        iqr<-centiles.pred(plmr, type="centiles",  xname = time_elapsed, xvalues=c(mint:maxt),
                                           cent=c(2.5, 25, 50, 75, 97.5),
                                           data=matchmodel,
                                           plot=FALSE)
                    } else {
                        iqr<-centiles.pred(plmr, type="centiles",  xname = time_elapsed, xvalues=c(time_window[1]:time_window[2]),
                                           cent=c(2.5, 25,50,75, 97.5),
                                           data=matchmodel,
                                           plot=FALSE)
                    }
                    
                    
                    # - - - - - - - - - - - - - - - - - - - - - - #
                    # If not LOOCV, meaning if we're doing just predictions obtain test data bias, coverage and precision
                    # - - - - - - - - - - - - - - - - - - - - - - #
                    if(loocv!=TRUE){
                        
                        #
                        # - - - Zscore on test set
                        #
                        
                        #-- first need to make the test set such that 
                        zscf <- function(){
                            out <- tryCatch(
                                {
                                    #message("ZSCORE PREDICTION TRY PART")
                                    #message(paste0("PREDICTING FOR TRAIN = ",i))
                                    testpred=data.frame()
                                    # - - - - - - - - - - - - - - - - - - - - - - #
                                    # Iterate through each testing obs(x) closest to the current train obs(i) 
                                    # get predicted values using centles.pred()
                                    # - - - - - - - - - - - - - - - - - - - - - - #
                                    for(x in (test_post %>%
                                              distinct_(.dots="patient_id") %>%
                                              .[which(traintestmatchdf$nnarraytest[1,]  == {
                                                  train_post %>%
                                                      distinct_(.dots="patient_id") %>%
                                                      .[i,] %>% unlist %>% as.vector}
                                              ),]  %>% unlist %>% as.vector)){
                                        #message(paste0("PREDICTING FOR TEST = ",x))
                                        testpred <- testpred %>% 
                                            bind_rows(
                                                data.frame(
                                                    zsc = centiles.pred(plmr, type="z-scores",
                                                                        xname="time",
                                                                        data=matchmodel,
                                                                        xvalues=test_post %>%
                                                                            filter(.data$patient_id %in% x) %>%
                                                                            dplyr::select(.data$time) %>%
                                                                            unlist %>% as.vector,
                                                                        yval=test_post %>%
                                                                            filter(.data$patient_id %in% x) %>%
                                                                            dplyr::select_(outcome) %>%
                                                                            unlist %>% as.vector
                                                    ),
                                                    test_id = rep(x, length(test_post %>%
                                                                                filter(.data$patient_id %in% x) %>%
                                                                                dplyr::select(.data$time) %>%
                                                                                unlist %>% as.vector
                                                    )),
                                                    time = test_post %>%
                                                        filter(.data$patient_id %in% x) %>%
                                                        dplyr::select(.data$time) %>%
                                                        unlist %>% as.vector
                                                ))
                                    }
                                    
                                    return(testpred)
                                },
                                error = function(cond)
                                {
                                    message(paste("ERROR in ZSCORE PREDICTION"))
                                    #message("ORIGINAL ERROR:")
                                    message(cond)
                                    return(NA)
                                },
                                warning=function(cond)
                                {
                                    #message(paste("GAMLSS WARNING"))
                                    #message("ORIIGNAL WARNING:")
                                    message(cond)
                                    testpred=data.frame()
                                    # - - - - - - - - - - - - - - - - - - - - - - #
                                    # Iterate through each testing obs(x) closest to the current train obs(i) 
                                    # get predicted values using centles.pred()
                                    # - - - - - - - - - - - - - - - - - - - - - - #
                                    for(x in (test_post %>%
                                              distinct_(.dots="patient_id") %>%
                                              .[which(traintestmatchdf$nnarraytest[1,]  == {
                                                  train_post %>%
                                                      distinct_(.dots="patient_id") %>%
                                                      .[i,] %>% unlist %>% as.vector}
                                              ),]  %>% unlist %>% as.vector)){
                                        #message(paste0("PREDICTING FOR TEST = ",x))
                                        testpred <- testpred %>% 
                                            bind_rows(
                                                data.frame(
                                                    zsc = centiles.pred(plmr, type="z-scores",
                                                                        xname="time",
                                                                        data=matchmodel,
                                                                        xvalues=test_post %>%
                                                                            filter(.data$patient_id %in% x) %>%
                                                                            dplyr::select(.data$time) %>%
                                                                            unlist %>% as.vector,
                                                                        yval=test_post %>%
                                                                            filter(.data$patient_id %in% x) %>%
                                                                            dplyr::select_(outcome) %>%
                                                                            unlist %>% as.vector
                                                    ),
                                                    test_id = rep(x, length(test_post %>%
                                                                                filter(.data$patient_id %in% x) %>%
                                                                                dplyr::select(.data$time) %>%
                                                                                unlist %>% as.vector
                                                    )),
                                                    time = test_post %>%
                                                        filter(.data$patient_id %in% x) %>%
                                                        dplyr::select(.data$time) %>%
                                                        unlist %>% as.vector
                                                ))
                                    }
                                    
                                    return(testpred)
                                    #return(NULL)
                                },
                                finally={
                                    #message("PROCESSING GAMLSS")
                                    #message("GAMLSS EOL SUCCESS")
                                }
                            )
                            return(out)
                        }
                        zsc <- zscf()
                        if(!is.data.frame(zsc)) {
                            print(zsc)
                            print(str(zsc))
                            message(zsc)
                            message(str(zsc))
                            
                            message('Something wrong with Z-SCORE prediction. Prediction not included')
                            zsc_list[[cnt]] <- data.frame(zsc = NA,
                                                          test_id = NA,
                                                          time = NA
                            )
                        } else {
                            zsc_list[[cnt]] <- zsc
                        }
                        
                        # - - - - --  - - - - - - -#
                        # Calculate and store Test-set bia, coverage, and precision
                        # - - - - --  - - - - - - -#
                        
                        # -- IQR for the test data set
                        
                        # iqr$iqr<-iqr$C75-iqr$C25
                        # iqrvec[cnt]<-mean(iqr$iqr)
                        
                        # select the test id's that corresopnd to current train patient
                        targetid<-test_post %>%
                            distinct_(.dots="patient_id") %>%
                            .[which(traintestmatchdf$nnarraytest[1,]  == {
                                train_post %>%
                                    distinct_(.dots="patient_id") %>%
                                    .[i,] %>% unlist %>% as.vector}
                            ),]  %>% unlist %>% as.vector
                        
                        # get the trainging set post op data
                        # targetrec<-test_post[which(test_post$patient_id %in% targetid), ] 
                        
                        # bias<-merge(iqr,targetrec, by=time_elapsed)
                        # bias$diff<-bias$C50-bias[,outcome]
                        
                        #store mean of bias and rmse in biasvec, rmsevec
                        # biasvec[cnt] <- mean(bias$diff)
                        # rmsevec[cnt] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))
                        
                        # coverage prob based on IQR 50%
                        # bias$cov1<-bias[,outcome]-bias$C25
                        # bias$cov2<-bias$C75-bias[,outcome]
                        # bias$cov3<-ifelse(bias$cov1>0,1,0)
                        # bias$cov4<-ifelse(bias$cov2>0,1,0)
                        # bias$cov5<-bias$cov3+bias$cov4
                        # bias$coverage<-ifelse(bias$cov5>1,1,0)
                        
                        # store mean of the n coverage in vector
                        # coveragevec[cnt]<-mean(bias$coverage)
                        
                        
                        if(any(iqr$C50 > thresh_val)){
                            crazymatch[[cnt]] <- matchmodel
                        } else {
                            crazymatch[[cnt]] <- NA
                        }
                        
                        #-- precision
                        # precisionvec[[cnt]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 
                        #-- store Testing C50
                        dfList_test[[cnt]] <- test_post[which(test_post$patient_id %in% targetid), c("patient_id",time_elapsed,outcome)] %>%
                            left_join(
                                data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                                by = "time"
                            )
                        
                        # #-- store Training C50
                        # dfList[[cnt]] <- train_post[which(train_post$patient_id %in% train[i, patid][[1]]), c("patient_id",time_elapsed,outcome)] %>%
                        #     left_join(
                        #         data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                        #         by = "time"
                        #     )
                        #-- store all centile for later test set merging
                        # centilepred[[cnt]] <- cbind(train[i, patid][[1]], iqr$time, iqr$C50)
                        
                    } else { # loocv = TRUE
                        
                        # - - - -
                        # z-scores vector
                        # - - - - 
                        zscf <- function(){
                            out <- tryCatch(
                                {
                                    #message("ZSCORE PREDICTION LOOCV TRY PART")
                                    trainzsc <- data.frame(
                                        zsc = centiles.pred(plmr, type="z-scores",
                                                            xname=time_elapsed,
                                                            data=matchmodel,
                                                            xvalues=train_post[train_post$patient_id %in% train[i, patid][[1]],time_elapsed][[1]],
                                                            yval=train_post[train_post$patient_id %in% train[i, patid][[1]],outcome][[1]]),
                                        train_id = rep(train[i, patid][[1]], length(train_post[train_post$patient_id %in% train[i, patid][[1]],time_elapsed][[1]]))
                                        ,
                                        time = train_post[train_post$patient_id %in% train[i, patid][[1]],time_elapsed][[1]]
                                    )
                                    return(trainzsc)
                                    
                                },
                                error = function(cond)
                                {
                                    message(paste("ERROR in ZSCORE PREDICTION"))
                                    #message("ORIGINAL ERROR:")
                                    message(cond)
                                    return(NA)
                                },
                                warning=function(cond)
                                {
                                    #message(paste("GAMLSS WARNING"))
                                    #message("ORIIGNAL WARNING:")
                                    message(cond)
                                    trainzsc <- data.frame(
                                        zsc = centiles.pred(plmr, type="z-scores",
                                                            xname=time_elapsed,
                                                            data=matchmodel,
                                                            xvalues=train_post[train_post$patient_id %in% train[i, patid][[1]],time_elapsed][[1]],
                                                            yval=train_post[train_post$patient_id %in% train[i, patid][[1]],outcome][[1]]),
                                        train_id = rep(train[i, patid][[1]], length(train_post[train_post$patient_id %in% train[i, patid][[1]],time_elapsed][[1]]))
                                        ,
                                        time = train_post[train_post$patient_id %in% train[i, patid][[1]],time_elapsed][[1]]
                                    )
                                    return(trainzsc)
                                    #return(NULL)
                                },
                                finally={
                                    #message("PROCESSING GAMLSS")
                                    #message("GAMLSS EOL SUCCESS")
                                }
                            )
                            return(out)
                        }
                        zsc <- zscf()
                        if(!is.data.frame(zsc)) {
                            
                            message('Something wrong with Z-SCORE prediction. Prediction not included')
                            #zsc_list[[cnt]] <- list(train = NA)
                            zsc_list[[cnt]] <- data.frame(zsc = NA,
                                                          test_id = NA,
                                                          time = NA
                            )
                            
                        } else {
                            zsc_list[[cnt]] <-  zsc
                            #zsc_list[[cnt]] <- list(train = zsc)
                            # below is when we want to also include the time point for the individual which is probably iportant... 
                        }
                        
                        message(paste0("Getting MLE  predictions: "))
                        
                        # - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - #
                        # code for storing either LOOCV result or extracting prediction #
                        # - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - #
                        
                        # iqr$iqr<-iqr$C75-iqr$C25
                        
                        # iqrvec[cnt]<-mean(iqr$iqr)
                        targetid<-train[n,patid][[1]]
                        # targetrec<-train_post[which(train_post$patient_id %in% targetid), ]
                        
                        # bias<-merge(iqr,targetrec, by=time_elapsed)
                        # bias$diff<-bias$C50-bias[,outcome]
                        
                        #store mean of bias and rmse in biasvec, rmsevec
                        # biasvec[cnt] <- mean(bias$diff)
                        # rmsevec[cnt] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))
                        # 
                        # # coverage prob based on IQR 50%
                        # bias$cov1<-bias[,outcome]-bias$C25
                        # bias$cov2<-bias$C75-bias[,outcome]
                        # bias$cov3<-ifelse(bias$cov1>0,1,0)
                        # bias$cov4<-ifelse(bias$cov2>0,1,0)
                        # bias$cov5<-bias$cov3+bias$cov4
                        # bias$coverage<-ifelse(bias$cov5>1,1,0)
                        # 
                        # # store mean of the n coverage in vector
                        # coveragevec[cnt]<-mean(bias$coverage)
                        # 
                        # # coverage prob based on IQR 95%
                        # bias$cov1<-bias[,outcome]-bias$C2.5
                        # bias$cov2<-bias$C97.5 -bias[,outcome]
                        # bias$cov3<-ifelse(bias$cov1>0,1,0)
                        # bias$cov4<-ifelse(bias$cov2>0,1,0)
                        # bias$cov5<-bias$cov3+bias$cov4
                        # bias$coverage<-ifelse(bias$cov5>1,1,0)
                        # 
                        # # store mean of the n coverage in vector
                        # coveragevec95a[cnt]<-mean(bias$coverage)
                        # #-- store Testing C50
                        # dfList_test[[cnt]] <- test_post[which(test_post$patient_id %in% targetid), c("patient_id",time_elapsed,outcome)] %>%
                        #     left_join(
                        #         data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                        #         by = "time"
                        #     )
                        
                        #-- store Training C50
                        dfList[[cnt]] <- train_post[which(train_post$patient_id %in% train[i, patid][[1]]), c("patient_id",time_elapsed,outcome)] %>%
                            left_join(
                                data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                                by = "time"
                            )
                        
                        #-- precision potentially remove later because this is 7.5mb per n so 7.5*14 ~ 100MB
                        # precisionvec[[cnt]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 
                        
                        
                    }
                    
                }
                message(paste0("Current count is: ",cnt))
                message(paste0("and Current n: ",i))
            }
            
            # - - - - - - - - - - - - - - - - - - - - - - - - -#  
            # Code to create dataframe of the time and C75-C25 #
            # - - - - - - - - - - - - - - - - - - - - - - - - -#  
            
            message(paste0("Merging measures: "))
            if(loocv!=TRUE){
                message(paste0("Constructing Prediction Data"))
                
                nn_arr  <- list(
                    pred_train = dfList, 
                    pred_test = dfList_test, 
                    # biasvec = Filter(Negate(is.na),Filter(Negate(is.na), biasvec)),
                    # coveragevec = Filter(Negate(is.na),Filter(Negate(is.na), coveragevec)),
                    # centilerange = centilepred, 
                    # precisionvec=precisionvec, 
                    zsc_list = zsc_list,
                    # rmse = rmsevec,
                    crazymatch = crazymatch)
                
                
            } else {
                
                # Get rid of NA or NULL values that are created
                nn_arr <- list(
                    #                          Filter(Negate(is.null),Filter(Negate(is.null), fin)),
                    dfList,
                    dfList_test,
                    # Filter(Negate(is.na),Filter(Negate(is.na), biasvec)),
                    # Filter(Negate(is.na),Filter(Negate(is.na), coveragevec)),
                    # Filter(Negate(is.na),Filter(Negate(is.na), coveragevec95a)),
                    #Filter(Negate(is.na),Filter(Negate(is.na), coveragevec95)),
                    # Filter(Negate(is.na),Filter(Negate(is.na), iqrvec)),
                    # rmsevec,
                    zsc_list,
                    # precisionvec = precisionvec,
                    misses 
                )
                
                # name the list objects
                message(paste0("Assigning names"))
                #names(All_list) <- c("bias","iqrcoverage","coverage95c","coverage95m","iqr","rmse", "dropped_cases")
                names(nn_arr) <- c(
                    "pred_train","pred_test",
                    # "dfList","dfList_test",
                    # "bias","iqrcoverage","coverage95c","iqr","rmse",
                    "zscore",
                    # "precisionvec", 
                    "dropped_cases")
            }
            nn_arr
            
        }
        
        # after loocv finishes for all nearest, print which n is the lowest
        
        
        stopImplicitCluster()
        # return either LOOCV array or final prediction array
        return(loocvres)
        
    } else {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Non parallel solution
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Empty list that will be populated with performance measures and predictions
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nn_arr = list()
        
        for(n in nearest){
            cnt = 0
            misses = 0
            
            dfList <- list()                #-- list to store training data's predicted C50 from gamlss model
            dfList_test <- list()           #-- list to store testing data's predicted C50 from gamlss model
            # centilepred <- list()           #-- store all centile for later test set merging
            # biasvec<-vector()               #-- store mean of bias
            # rmsevec <- vector()             #-- store rmse vector
            # coveragevec<-vector()           #-- store coverage vector
            # coveragevec95a<-vector()        #-- store mean of the n coverage in vector
            #coveragevec95<-vector()
            # iqrvec<-vector()
            # ninefiveconf <- data.frame()
            # precisionvec <- list()
            crazymatch <- list()
            zsc_list <- list()
            
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # iterate through all patients and match patients according to order. ord <- data is the training data
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if(is.null(interval)){
                patlistint <- seq(1,nrow(train))
            } else {
                patlistint <- seq(1,nrow(train), by=interval)
            }
            
            for (i in patlistint) {
                cnt = cnt + 1
                
                # - - - - - - - - - - - - - - - - - - - - - - #
                # Matching cases using absolute value of difference between Fitted values
                # Matching by probability weighting also possible
                # - - - - - - - - - - - - - - - - - - - - - - #
                
                if(loocv){
                    
                    matchmodel <- train_post %>%
                        filter(.data$patient_id %in% (traintestmatchdf$nnarraytrain[,i] %>%
                                                          head(n=n)))
                } else {
                    
                    matchmodel <- train_post %>%
                        filter(.data$patient_id %in% c(train[[patid]][i], traintestmatchdf$nnarraytrain[,i] %>%
                                                          head(n=(n-1))))
                    
                }
                
                # - - - - - - - - - - - - - - - - - - - - - - #
                # Fit GAMLSS model to nearest n matches
                # - - - - - - - - - - - - - - - - - - - - - - #
                
                plm <- function() {
                    
                    out <- tryCatch(
                        {
                            #message("TRY PART")
                            update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE))
                            
                        },
                        error = function(cond)
                        {
                            message(paste("ERROR in GAMLSS"))
                            #message("ORIGINAL ERROR:")
                            message(cond)
                            return(NA)
                        },
                        warning=function(cond)
                        {
                            #message(paste("GAMLSS WARNING"))
                            #message("ORIIGNAL WARNING:")
                            message(cond)
                            return(update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE)))
                            #return(NULL)
                        },
                        finally={
                            #message("PROCESSING GAMLSS")
                            #message("GAMLSS EOL SUCCESS")
                        }
                    )
                    return(out)
                }
                plmr <- plm()
                
                # if something wrong with iteration... then just don't do prediction it eats up all the memory anyways
                if(typeof(plmr) != 'list') {
                    
                    misses = misses + 1 # increment number of misses
                    message('Something wrong with plm model. Prediction not included')
                    
                } else {
                    
                    
                    message(paste0("Getting predictions: "))
                    
                    # time window again
                    if(is.null(time_window)){
                        iqr<-centiles.pred(plmr, type="centiles",  xname = time_elapsed, xvalues=c(mint:maxt),
                                           cent=c(2.5, 25, 50, 75, 97.5),
                                           data=matchmodel,
                                           plot=FALSE)
                    } else {
                        iqr<-centiles.pred(plmr, type="centiles",  xname = time_elapsed, xvalues=c(time_window[1]:time_window[2]),
                                           cent=c(2.5, 25,50,75, 97.5),
                                           data=matchmodel,
                                           plot=FALSE)
                    }
                    
                    # - - - - - - - - - - - - - - - - - - - - - - #
                    # If not LOOCV, meaning if we're doing just predictions obtain test data bias, coverage and precision
                    # - - - - - - - - - - - - - - - - - - - - - - #
                    if(loocv!=TRUE){
                        
                        
                        
                        #
                        # - - - Zscore on test set
                        #
                        
                        #-- first need to make the test set such that 
                        zscf <- function(){
                            
                            out <- tryCatch(
                                {
                                    
                                    #message("ZSCORE PREDICTION TRY PART")
                                    #message(paste0("PREDICTING FOR TRAIN = ",i))
                                    #message(paste0("PREDICTING FOR TRAIN = ",train[i,patid][[1]]))
                                    testpred=data.frame()
                                    # - - - - - - - - - - - - - - - - - - - - - - #
                                    # Iterate through each testing obs(x) closest to the current train obs(i) 
                                    # get predicted values using centles.pred()
                                    # - - - - - - - - - - - - - - - - - - - - - - #
                                    for(x in (test_post %>%
                                              distinct_(.dots="patient_id") %>%
                                              .[which(traintestmatchdf$nnarraytest[1,]  == {
                                                  train_post %>%
                                                      distinct_(.dots="patient_id") %>%
                                                      .[i,] %>% unlist %>% as.vector}
                                              ),]  %>% unlist %>% as.vector)){
                                        #message(paste0("PREDICTING FOR TEST = ",x))
                                        testpred <- testpred %>% 
                                            bind_rows(
                                                data.frame(
                                                    zsc = centiles.pred(plmr, type="z-scores",
                                                                        xname="time",
                                                                        data=matchmodel,
                                                                        xvalues=test_post %>%
                                                                            filter(.data$patient_id %in% x) %>%
                                                                            dplyr::select(.data$time) %>%
                                                                            unlist %>% as.vector,
                                                                        yval=test_post %>%
                                                                            filter(.data$patient_id %in% x) %>%
                                                                            dplyr::select_(outcome) %>%
                                                                            unlist %>% as.vector
                                                    ),
                                                    test_id = rep(x, length(test_post %>%
                                                                                filter(.data$patient_id %in% x) %>%
                                                                                dplyr::select(.data$time) %>%
                                                                                unlist %>% as.vector
                                                    )),
                                                    time = test_post %>%
                                                        filter(.data$patient_id %in% x) %>%
                                                        dplyr::select(.data$time) %>%
                                                        unlist %>% as.vector
                                                ))
                                    }
                                    
                                    return(testpred)
                                    
                                },
                                error = function(cond)
                                {
                                    message(paste("ERROR in ZSCORE PREDICTION"))
                                    #message("ORIGINAL ERROR:")
                                    message(cond)
                                    return(NA)
                                },
                                warning=function(cond)
                                {
                                    #message(paste("GAMLSS WARNING"))
                                    #message("ORIIGNAL WARNING:")
                                    message(cond)
                                    
                                    testpred=data.frame()
                                    # - - - - - - - - - - - - - - - - - - - - - - #
                                    # Iterate through each testing obs(x) closest to the current train obs(i) 
                                    # get predicted values using centles.pred()
                                    # - - - - - - - - - - - - - - - - - - - - - - #
                                    
                                    for(x in (test_post %>%
                                              distinct_(.dots="patient_id") %>%
                                              .[which(traintestmatchdf$nnarraytest[1,]  == {
                                                  train_post %>%
                                                      distinct_(.dots="patient_id") %>%
                                                      .[i,] %>% unlist %>% as.vector}
                                              ),]  %>% unlist %>% as.vector)){
                                        #message(paste0("PREDICTING FOR TEST = ",x))
                                        testpred <- testpred %>% 
                                            bind_rows(
                                                data.frame(
                                                    zsc = centiles.pred(plmr, type="z-scores",
                                                                        xname="time",
                                                                        data=matchmodel,
                                                                        xvalues=test_post %>%
                                                                            filter(.data$patient_id %in% x) %>%
                                                                            dplyr::select(.data$time) %>%
                                                                            unlist %>% as.vector,
                                                                        yval=test_post %>%
                                                                            filter(.data$patient_id %in% x) %>%
                                                                            dplyr::select_(outcome) %>%
                                                                            unlist %>% as.vector
                                                    ),
                                                    test_id = rep(x, length(test_post %>%
                                                                                filter(.data$patient_id %in% x) %>%
                                                                                dplyr::select(.data$time) %>%
                                                                                unlist %>% as.vector
                                                    )),
                                                    time = test_post %>%
                                                        filter(.data$patient_id %in% x) %>%
                                                        dplyr::select(.data$time) %>%
                                                        unlist %>% as.vector
                                                ))
                                    }
                                    
                                    return(testpred)
                                    #return(NULL)
                                },
                                finally={
                                    #message("PROCESSING GAMLSS")
                                    #message("GAMLSS EOL SUCCESS")
                                }
                                
                            )
                            return(out)
                        }
                        
                        zsc <- zscf()
                        if(!is.data.frame(zsc)) {
                            print(zsc)
                            print(str(zsc))
                            message(zsc)
                            message(str(zsc))
                            
                            message('Something wrong with Z-SCORE prediction. Prediction not included')
                            zsc_list[[cnt]] <- data.frame(zsc = NA,
                                                          test_id = NA,
                                                          time = NA
                            )
                        } else {
                            zsc_list[[cnt]] <- zsc
                        }
                        
                        # - - - - --  - - - - - - -#
                        # Calculate and store Test-set bia, coverage, and precision
                        # - - - - --  - - - - - - -#
                        
                        # -- IQR for the test data set
                        
                        # iqr$iqr<-iqr$C75-iqr$C25
                        # iqrvec[cnt]<-mean(iqr$iqr)
                        
                        # select the test id's that corresopnd to current train patient
                        targetid<-test_post %>%
                            distinct_(.dots="patient_id") %>%
                            .[which(traintestmatchdf$nnarraytest[1,]  == {
                                train_post %>%
                                    distinct_(.dots="patient_id") %>%
                                    .[i,] %>% unlist %>% as.vector}
                            ),]  %>% unlist %>% as.vector
                        # targetrec<-test_post[which(test_post$patient_id %in% targetid), ] # ge tthe trainging set post op data
                        
                        # bias<-merge(iqr,targetrec, by=time_elapsed)
                        # bias$diff<-bias$C50-bias[,outcome]
                        
                        #store mean of bias and rmse in biasvec, rmsevec
                        # biasvec[cnt] <- mean(bias$diff)
                        # rmsevec[cnt] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))
                        
                        # coverage prob based on IQR 50%
                        # bias$cov1<-bias[,outcome]-bias$C25
                        # bias$cov2<-bias$C75-bias[,outcome]
                        # bias$cov3<-ifelse(bias$cov1>0,1,0)
                        # bias$cov4<-ifelse(bias$cov2>0,1,0)
                        # bias$cov5<-bias$cov3+bias$cov4
                        # bias$coverage<-ifelse(bias$cov5>1,1,0)
                        
                        # store mean of the n coverage in vector
                        # coveragevec[cnt]<-mean(bias$coverage)
                        
                        if(any(iqr$C50 > thresh_val)){
                            crazymatch[[cnt]] <- matchmodel
                        } else {
                            crazymatch[[cnt]] <- NA
                        }
                        
                        #-- precision
                        # precisionvec[[cnt]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 
                        #-- store Testing C50
                        dfList_test[[cnt]] <- test_post[which(test_post$patient_id %in% targetid), c("patient_id",time_elapsed,outcome)] %>%
                            #train <- post[train <- post$patient <- id == ord <- data$id[c(i)],c("patient <- id",time <- elapsed,outcome)] %>% 
                            left_join(
                                data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                                by = "time"
                            )
                        #left_join(
                        #data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, patient_id = iqr[,patient_id]) ,
                        #by = "patient_id"
                        #)
                        
                        # #-- store Training C50
                        # dfList[[cnt]] <- train_post[which(train_post$patient_id %in% train[i,patid][[1]]), c("patient_id",time_elapsed,outcome)] %>%
                        #     left_join(
                        #         data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                        #         by = "time"
                        #     )
                        #left_join(
                        #data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, patient_id = iqr[,]) ,
                        #by = "patient_id"
                        #)
                        # -- store all centile for later test set merging
                        # centilepred[[cnt]] <- cbind(train[i,patid][[1]], iqr$time, iqr$C50)
                        
                    } else { # loocv = TRUE
                        
                        # - - - -
                        # z-scores vector
                        # - - - - 
                        zscf <- function(){
                            
                            out <- tryCatch(
                                {
                                    message("ZSCORE PREDICTION LOOCV TRY PART")
                                    trainzsc <- data.frame(
                                        zsc = centiles.pred(plmr, type="z-scores",
                                                            xname=time_elapsed,
                                                            data=matchmodel,
                                                            xvalues=train_post[train_post$patient_id %in% train[i,patid][[1]],time_elapsed][[1]],
                                                            yval=train_post[train_post$patient_id %in% train[i,patid][[1]],outcome][[1]]),
                                        train_id = rep(train[i,patid][[1]], length(train_post[train_post$patient_id %in% train[i,patid][[1]],time_elapsed][[1]]))
                                        ,
                                        time = train_post[train_post$patient_id %in% train[i,patid][[1]],time_elapsed][[1]]
                                    )
                                    return(trainzsc)
                                    
                                },
                                error = function(cond)
                                {
                                    #message(paste("ERROR in ZSCORE PREDICTION"))
                                    #message("ORIGINAL ERROR:")
                                    message(cond)
                                    return(NA)
                                },
                                warning=function(cond)
                                {
                                    #message(paste("GAMLSS WARNING"))
                                    #message("ORIIGNAL WARNING:")
                                    message(cond)
                                    trainzsc <- data.frame(
                                        zsc = centiles.pred(plmr, type="z-scores",
                                                            xname=time_elapsed,
                                                            data=matchmodel,
                                                            xvalues=train_post[train_post$patient_id %in% train[i,patid][[1]],time_elapsed][[1]],
                                                            yval=train_post[train_post$patient_id %in% train[i,patid][[1]],outcome][[1]]),
                                        train_id = rep(train[i,patid][[1]], length(train_post[train_post$patient_id %in% train[i,patid][[1]],time_elapsed][[1]]))
                                        ,
                                        time = train_post[train_post$patient_id %in% train[i,patid][[1]],time_elapsed][[1]]
                                    )
                                    return(trainzsc)
                                    
                                    #return(NULL)
                                },
                                finally={
                                    #message("PROCESSING GAMLSS")
                                    #message("GAMLSS EOL SUCCESS")
                                }
                            )
                            return(out)
                            
                        }
                        zsc <- zscf()
                        if(!is.data.frame(zsc)) {
                            
                            message('Something wrong with Z-SCORE prediction. Prediction not included')
                            #zsc <- list[[cnt]] <- list(train = NA)
                            zsc_list[[cnt]] <- data.frame(zsc = NA,
                                                          test_id = NA,
                                                          time = NA
                            )
                            
                        } else {
                            zsc_list[[cnt]] <-  zsc
                            #zsc <- list[[cnt]] <- list(train = zsc)
                            # below is when we want to also include the time point for the individual which is probably iportant... 
                            #zsc <- list[[cnt]] <- list(zsc = zsc, time = train <- post[train <- post$patient <- id %in% ord <- data$id[c(i)],time <- elapsed][[1]])
                        }
                        
                        message(paste0("Getting MLE  predictions: "))
                        
                        # - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - #
                        # code for storing either LOOCV result or extracting prediction #
                        # - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - #
                        
                        # iqr$iqr<-iqr$C75-iqr$C25
                        
                        # iqrvec[cnt]<-mean(iqr$iqr)
                        targetid<-train[i,patid][[1]]
                        # targetrec<-train_post[which(train_post$patient_id %in% targetid), ]
                        
                        # bias<-merge(iqr,targetrec, by=time_elapsed)
                        # bias$diff<-bias$C50-bias[,outcome]
                        
                        #store mean of bias and rmse in biasvec, rmsevec
                        # biasvec[cnt] <- mean(bias$diff)
                        # rmsevec[cnt] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))
                        
                        # coverage prob based on IQR 50%
                        # bias$cov1<-bias[,outcome]-bias$C25
                        # bias$cov2<-bias$C75-bias[,outcome]
                        # bias$cov3<-ifelse(bias$cov1>0,1,0)
                        # bias$cov4<-ifelse(bias$cov2>0,1,0)
                        # bias$cov5<-bias$cov3+bias$cov4
                        # bias$coverage<-ifelse(bias$cov5>1,1,0)
                        
                        # store mean of the n coverage in vector
                        # coveragevec[cnt]<-mean(bias$coverage)
                        
                        # coverage prob based on IQR 95%
                        # bias$cov1<-bias[,outcome]-bias$C2.5
                        # bias$cov2<-bias$C97.5 -bias[,outcome]
                        # bias$cov3<-ifelse(bias$cov1>0,1,0)
                        # bias$cov4<-ifelse(bias$cov2>0,1,0)
                        # bias$cov5<-bias$cov3+bias$cov4
                        # bias$coverage<-ifelse(bias$cov5>1,1,0)
                        
                        # store mean of the n coverage in vector
                        # coveragevec95a[cnt]<-mean(bias$coverage)
                        #-- store Testing C50
                        # dfList_test[[cnt]] <- test_post[which(test_post$patient_id %in% targetid), c("patient_id",time_elapsed,outcome)] %>%
                        #     #train <- post[train <- post$patient <- id == ord <- data$id[c(i)],c("patient <- id",time <- elapsed,outcome)] %>% 
                        #     left_join(
                        #         data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                        #         by = "time"
                        #     )
                        #left_join(
                        #data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, patient_id = iqr[,patient_id]) ,
                        #by = "patient_id"
                        #)
                        
                        #-- store Training C50
                        dfList[[cnt]] <- train_post[which(train_post$patient_id %in% train[i,patid][[1]]), c("patient_id",time_elapsed,outcome)] %>%
                            left_join(
                                data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                                by = "time"
                            )
                        #left_join(
                        #data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, patient_id = iqr[,]) ,
                        #by = "patient_id"
                        #)
                        
                        #-- precision potentially remove later because this is 7.5mb per n so 7.5*14 ~ 100MB
                        # precisionvec[[cnt]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 
                        
                    }
                    
                    
                }
                message(paste0("Current count is: ",cnt))
                message(paste0("and Current n: ",n))
                
            }
            
            # - - - - - - - - - - - - - - - - - - - - - - - - -#  
            # Code to create dataframe of the time and C75-C25 #
            # - - - - - - - - - - - - - - - - - - - - - - - - -#  
            
            message(paste0("Merging measures: "))
            if(loocv!=TRUE){
                message(paste0("Constructing Prediction Data"))
                
                nn_arr  <- list(
                    pred_train = dfList, 
                    pred_test = dfList_test, 
                    # biasvec = Filter(Negate(is.na),Filter(Negate(is.na), biasvec)),
                    # coveragevec = Filter(Negate(is.na),Filter(Negate(is.na), coveragevec)),
                    # centilerange = centilepred, 
                    # precisionvec=precisionvec, 
                    zsc_list = zsc_list,
                    # rmse = rmsevec,
                    crazymatch = crazymatch)
                
            } else {
                
                # Get rid of NA or NULL values that are created
                All_list <- list(
                    #                          Filter(Negate(is.null),Filter(Negate(is.null), fin)),
                    dfList,
                    dfList_test,
                    # Filter(Negate(is.na),Filter(Negate(is.na), biasvec)),
                    # Filter(Negate(is.na),Filter(Negate(is.na), coveragevec)),
                    # Filter(Negate(is.na),Filter(Negate(is.na), coveragevec95a)),
                    # #Filter(Negate(is.na),Filter(Negate(is.na), coveragevec95)),
                    # Filter(Negate(is.na),Filter(Negate(is.na), iqrvec)),
                    # rmsevec,
                    zsc_list,
                    # precisionvec = precisionvec,
                    misses 
                )
                
                # name the list objects
                message(paste0("Assigning names"))
                #names(All <- list) <- c("bias","iqrcoverage","coverage95c","coverage95m","iqr","rmse", "dropped <- cases")
                names(All_list) <- c(
                    "pred_train","pred_test",
                    # "dfList","dfList_test",
                    # "bias","iqrcoverage","coverage95c","iqr","rmse",
                    "zscore",
                    # "precisionvec", 
                    "dropped_cases")
                
                message(paste0("Putting in the list in array slot: ",n))
                # store the result of the n nearest <- n result
                nn_arr[[which(nearest == n)]] = All_list
                
                # rename the list items to indicate nearestn
                names(nn_arr)[which(nearest == n)] <- paste0('nearest_',n)
                
                
                # select of the n result, which has the smallest mean 'rmse' or selection criteria chosen by user
                
            }
            
        }
        
        # after loocv finishes for all nearest <- n, print which n is the lowest
        
        # return either LOOCV array or final prediction array
        return(nn_arr)
        
    }
    
}
