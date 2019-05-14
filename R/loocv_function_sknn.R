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
#' @param \dots Passed down to \code{gamlss}
#' @return There are many possible return values 
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
                           biasm="raw",
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
    # REFERENCE GAMLSS SETTING
    # - - - - - - - - - - - - - - - - - - - - - # 
    if(cs){
        if(is.null(dfspec)) {
            spl=paste0("cs(")
            spls=paste0("cs(")
            spln=paste0("cs(")
            splt=paste0("cs(")
        } else {
            spl=paste0("cs(df=",d_f_m,",")
            spls=paste0("cs(df=",d_f_s,",")
            spln=paste0("cs(df=",d_f_n,",")
            splt=paste0("cs(df=",d_f_t,",")
        }
    } else {
        spl="pb("
    }

    # - - - - - - - - - - - - - - - - - - - - - # 
    # Setting distribution for GAMLSS
    # If NULL then use NO as default
    # - - - - - - - - - - - - - - - - - - - - - # 
    if(is.null(dist_fam)) {

        # - - - - - - - - - - - - - - - - - - - - - # 
        # Reference gamlss fitting 
        # - - - - - - - - - - - - - - - - - - - - - # 
        message("Distribution of GAMLSS not chosen")
        message("Will determine based on Normal, Gamma, and Box Cox Cole and Green")

        ref1<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m,")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    data=na.omit(train_post), family=NO,
                    #gamlss control
                    ...)
        ref2<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m,")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    data=na.omit(train_post), family=GA,
                    #gamlss control
                    ...)
        ref3<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = as.formula(paste0(" ~ ",spln, time_elapsed, ")")),
                    data=na.omit(train_post), family=BCCGo,
                    #gamlss control
                    ...)
        ref4<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = as.formula(paste0(" ~ ",spln, time_elapsed, ")")),
                    data=na.omit(train_post), family=BCTo,
                    #gamlss control
                    ...)
        ref5<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = as.formula(paste0(" ~ ",spln, time_elapsed, ")")),
                    data=na.omit(train_post), family=BCPEo,
                    #gamlss control
                    ...)
        # here the BCCGo distribution with tau is intentionally left out bc the fit is good but not as when we jsut do the nu
        #ref5<-gamlss(eval(as.formula(paste0(outcome," ~ ", "pb(", time_elapsed, ")"))),
                    #sigma.formula = as.formula(paste0(" ~ ", "pb(", time_elapsed, ")")),
                    #nu.formula = as.formula(paste0(" ~ ", "pb(", time_elapsed, ")")),
                    #tau.formula = as.formula(paste0(" ~ ", "pb(", time_elapsed, ")")),
                    #data=na.omit(train_post), family=BCCGo)

        ref <- list(ref1,ref2,ref3,ref4,ref5)[which.min(lapply(list(ref1,ref2,ref3,ref4,ref5), function(x){x$aic}))][[1]]
    } else if(length(dist_fam()$parameters) == 2) {
        ref<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    data=data.frame(na.omit(train_post)), family=dist_fam,
                    #gamlss control
                    ...)
    } else if(length(dist_fam()$parameters) == 3) {
        ref<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = ~1,
                    data=na.omit(train_post), family=dist_fam,
                    #gamlss control
                    ...)
    } else if(length(dist_fam()$parameters) == 4) {
        ref<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = ~1,
                    tau.formula = ~1,
                    data=na.omit(train_post), family=dist_fam,
                    ...
                    )
    }
    gamlss_dist <- ref$family[2]


    # - - - - - - - - - - - - - - - - - - - - - # 
    # Full training post operative data fitting (i.e reference data) with GAMLSS 
    # - - - - - - - - - - - - - - - - - - - - - # 

    # - - - - - - - - - - - - - - - - - - - - - # 
    # Specify time_window if NULL for GAMLSS centiles prediction
    # - - - - - - - - - - - - - - - - - - - - - # 
    if(is.null(time_window)){
        mint <- min(min(train_post[, time_elapsed]), min(test_post[, time_elapsed]))
        maxt <- max(max(train_post[, time_elapsed]), max(test_post[, time_elapsed]))
        iqrfull <- centiles.pred(ref, type='centiles', xname = time_elapsed, xvalues=c(mint:maxt),
                                data = na.omit(train_post),
                                cent=c(25,75), plot=FALSE)
    } else {
        iqrfull <- centiles.pred(ref, type="centiles", xname = time_elapsed, xvalues=c(time_window[1]:time_window[2]),
                                 data = na.omit(train_post),
                                 cent=c(25,75), plot=FALSE)
    }
    iqrfull$iqr<-iqrfull$C75-iqrfull$C25

    # - - - - - - - - - - - - - - - - - - - - - # 
    # NEAREST NEIGHBOR MATCHING
    # - - - - - - - - - - - - - - - - - - - - - # 

    # Outerloop will iterate through a list of nearest_n numbers
    # Innerloop will iterate through all patients and store the result of the whole iteration into an array


    # - - - - - - - - - - - - - - - - - - - - - # 
    # GENERATE ARRAY OF MATCHED TRAINING AND TESTING ID
    # - - - - - - - - - - - - - - - - - - - - - # 

    temp <- matchIdExtractsknn(
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

    pat_level_func <- function(nearest=nearest_n, loocv=loocv, parallel=parallel){

        if(!is.null(parallel)){
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Parallel solution
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            registerDoParallel(cores = parallel)
            #registerDoSNOW(makeCluster(parallel, type = "SOCK"))
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Loop based on the nearest_n specified as user input 
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            loocvres <- foreach(n=nearest,
                                .export=c("train_post", "test_post","spl","ptr_m", "d_f_m", "d_f_s",
                                           "d_f_n", "d_f_t", "dist_fam", "spls","spln","splt",
                                          "outcome","time_elapsed", "plot", "matchprobweight","time_window",
                                          "interval","thresh_val","printtrace","userchoose", "ref",
                                          "loocv","mint","maxt", "temp"),                                                                                                              .packages=c("dplyr","gamlss"),
                                .combine=list,
                                .multicombine = TRUE
                                ) %dopar% {
                #for(n in nearest){
                cnt = 0
                misses = 0

                dfList <- list()                #-- list to store training data's predicted C50 from gamlss model
                dfList_test <- list()           #-- list to store testing data's predicted C50 from gamlss model
                centilepred <- list()           #-- store all centile for later test set merging
                biasvec<-vector()               #-- store mean of bias
                rmsevec <- vector()             #-- store rmse vector
                coveragevec<-vector()           #-- store coverage vector
                coveragevec95a<-vector()        #-- store mean of the n coverage in vector
                #coveragevec95<-vector()
                iqrvec<-vector()
                ninefiveconf <- data.frame()
                precisionvec <- list()
                crazymatch <- list()
                zsc_list <- list()

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # iterate through all patients and match patients according to order. ord_data is the training data
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                if(is.null(interval)){
                    patlistint <- seq(1,nrow(train[,patid]))
                } else {
                    patlistint <- seq(1,nrow(train[,patid]), by=interval)
                }
                for (i in patlistint) {
                    cnt = cnt + 1

                    # - - - - - - - - - - - - - - - - - - - - - - #
                    # Matching cases using absolute value of difference between Fitted values
                    # - - - - - - - - - - - - - - - - - - - - - - #

                    matchmodel <- train_post %>%
                        filter(.data$patient_id %in% (temp$nnarraytrain[,i] %>%
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
                                                              distinct_(.dots=patid) %>%
                                                              .[which(temp$nnarraytest[1,]  == {
                                                                          train_post %>%
                                                                              distinct_(.dots=patid) %>%
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
                                                              distinct_(.dots=patid) %>%
                                                              .[which(temp$nnarraytest[1,]  == {
                                                                          train_post %>%
                                                                              distinct_(.dots=patid) %>%
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

                            iqr$iqr<-iqr$C75-iqr$C25
                            iqrvec[cnt]<-mean(iqr$iqr)

                            # select the test id's that corresopnd to current train patient
                            targetid<-test_post %>%
                                distinct_(.dots=patid) %>%
                                .[which(temp$nnarraytest[1,]  == {
                                            train_post %>%
                                                distinct_(.dots="patient_id") %>%
                                                .[i,] %>% unlist %>% as.vector}
                                ),]  %>% unlist %>% as.vector

                            # get the trainging set post op data
                            targetrec<-test_post[which(test_post$patient_id %in% targetid), ] 

                            bias<-merge(iqr,targetrec, by=time_elapsed)
                            bias$diff<-bias$C50-bias[,outcome]

                            #store mean of bias and rmse in biasvec, rmsevec
                            biasvec[cnt] <- mean(bias$diff)
                            rmsevec[cnt] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))

                            # coverage prob based on IQR 50%
                            bias$cov1<-bias[,outcome]-bias$C25
                            bias$cov2<-bias$C75-bias[,outcome]
                            bias$cov3<-ifelse(bias$cov1>0,1,0)
                            bias$cov4<-ifelse(bias$cov2>0,1,0)
                            bias$cov5<-bias$cov3+bias$cov4
                            bias$coverage<-ifelse(bias$cov5>1,1,0)

                            # store mean of the n coverage in vector
                            coveragevec[cnt]<-mean(bias$coverage)


                            if(any(bias$C50 > thresh_val)){
                                crazymatch[[cnt]] <- matchmodel
                            } else {
                                crazymatch[[cnt]] <- NA
                            }

                            #-- precision
                            precisionvec[[cnt]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 
                            #-- store Testing C50
                            dfList_test[[cnt]] <- test_post[which(test_post$patient_id %in% targetid), c("patient_id",time_elapsed,outcome)] %>%
                                #train_post[train_post$patient_id == ord_data$id[c(i)],c("patient_id",time_elapsed,outcome)] %>% 
                                left_join(
                                          data.frame(time=iqr[,time_elapsed],c50 = iqr$C50) ,
                                          by = "time"
                                )

                            #-- store Training C50
                            dfList[[cnt]] <- train_post[which(train_post$patient_id %in% train[i, patid][[1]]), c("patient_id",time_elapsed,outcome)] %>%
                                left_join(
                                          data.frame(time=iqr[,time_elapsed],c50 = iqr$C50) ,
                                          by = "time"
                                )
                            #-- store all centile for later test set merging
                            centilepred[[cnt]] <- cbind(train[i, patid][[1]], iqr$time, iqr$C50)

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
                                #zsc_list[[cnt]] <- list(zsc = zsc, time = train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]])
                            }

                            message(paste0("Getting MLE  predictions: "))

                            # - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - #
                            # code for storing either LOOCV result or extracting prediction #
                            # - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - #

                            iqr$iqr<-iqr$C75-iqr$C25

                            iqrvec[cnt]<-mean(iqr$iqr)
                            targetid<-train[i,patid][[1]]
                            targetrec<-train_post[which(train_post$patient_id %in% targetid), ]

                            bias<-merge(iqr,targetrec, by=time_elapsed)
                            bias$diff<-bias$C50-bias[,outcome]

                            #store mean of bias and rmse in biasvec, rmsevec
                            biasvec[cnt] <- mean(bias$diff)
                            rmsevec[cnt] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))

                            # coverage prob based on IQR 50%
                            bias$cov1<-bias[,outcome]-bias$C25
                            bias$cov2<-bias$C75-bias[,outcome]
                            bias$cov3<-ifelse(bias$cov1>0,1,0)
                            bias$cov4<-ifelse(bias$cov2>0,1,0)
                            bias$cov5<-bias$cov3+bias$cov4
                            bias$coverage<-ifelse(bias$cov5>1,1,0)

                            # store mean of the n coverage in vector
                            coveragevec[cnt]<-mean(bias$coverage)

                            # coverage prob based on IQR 95%
                            bias$cov1<-bias[,outcome]-bias$C2.5
                            bias$cov2<-bias$C97.5 -bias[,outcome]
                            bias$cov3<-ifelse(bias$cov1>0,1,0)
                            bias$cov4<-ifelse(bias$cov2>0,1,0)
                            bias$cov5<-bias$cov3+bias$cov4
                            bias$coverage<-ifelse(bias$cov5>1,1,0)

                            # store mean of the n coverage in vector
                            coveragevec95a[cnt]<-mean(bias$coverage)

                            #-- precision potentially remove later because this is 7.5mb per n so 7.5*14 ~ 100MB
                            precisionvec[[cnt]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 


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
                                    biasvec = Filter(Negate(is.na),Filter(Negate(is.na), biasvec)),
                                    coveragevec = Filter(Negate(is.na),Filter(Negate(is.na), coveragevec)),
                                    centilerange = centilepred, 
                                    precisionvec=precisionvec, 
                                    zsc_list = zsc_list,
                                    rmse = rmsevec,
                                    crazymatch = crazymatch)


                } else {

                    # Get rid of NA or NULL values that are created
                    nn_arr <- list(
                                   #                          Filter(Negate(is.null),Filter(Negate(is.null), fin)),
                                   Filter(Negate(is.na),Filter(Negate(is.na), biasvec)),
                                   Filter(Negate(is.na),Filter(Negate(is.na), coveragevec)),
                                   Filter(Negate(is.na),Filter(Negate(is.na), coveragevec95a)),
                                   #Filter(Negate(is.na),Filter(Negate(is.na), coveragevec95)),
                                   Filter(Negate(is.na),Filter(Negate(is.na), iqrvec)),
                                   rmsevec,
                                   zsc_list,
                                   precisionvec = precisionvec,
                                   misses 
                    )

                    # name the list objects
                    message(paste0("Assigning names"))
                    #names(All_list) <- c("bias","iqrcoverage","coverage95c","coverage95m","iqr","rmse", "dropped_cases")
                    names(nn_arr) <- c("bias","iqrcoverage","coverage95c","iqr","rmse","zscore","precisionvec", "dropped_cases")
                }
                nn_arr

            }

            # after loocv finishes for all nearest_n, print which n is the lowest


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
                        centilepred <- list()           #-- store all centile for later test set merging
                        biasvec<-vector()               #-- store mean of bias
                        rmsevec <- vector()             #-- store rmse vector
                        coveragevec<-vector()           #-- store coverage vector
                        coveragevec95a<-vector()        #-- store mean of the n coverage in vector
                        #coveragevec95<-vector()
                        iqrvec<-vector()
                        ninefiveconf <- data.frame()
                        precisionvec <- list()
                        crazymatch <- list()
                        zsc_list <- list()
            
                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        # iterate through all patients and match patients according to order. ord <- data is the training data
                        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        if(is.null(interval)){
                            patlistint <- seq(1,nrow(train[,patid]))
                        } else {
                            patlistint <- seq(1,nrow(train[,patid]), by=interval)
                        }

                        for (i in patlistint) {
                            cnt = cnt + 1

                            # - - - - - - - - - - - - - - - - - - - - - - #
                            # Matching cases using absolute value of difference between Fitted values
                            # Matching by probability weighting also possible
                            # - - - - - - - - - - - - - - - - - - - - - - #

                            matchmodel <- train_post %>%
                                filter(.data$patient_id %in% (temp$nnarraytrain[,i] %>%
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
                                                            #message(paste0("PREDICTING FOR TRAIN = ",train[i,patid][[1]]))
                                                            testpred=data.frame()
                                                            # - - - - - - - - - - - - - - - - - - - - - - #
                                                            # Iterate through each testing obs(x) closest to the current train obs(i) 
                                                            # get predicted values using centles.pred()
                                                            # - - - - - - - - - - - - - - - - - - - - - - #

                                                            for(x in (test_post %>%
                                                                      distinct_(.dots=patid) %>%
                                                                      .[which(temp$nnarraytest[1,]  == {
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
                                                                      distinct_(.dots=patid) %>%
                                                                      .[which(temp$nnarraytest[1,]  == {
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

                                    iqr$iqr<-iqr$C75-iqr$C25
                                    iqrvec[cnt]<-mean(iqr$iqr)

                                    # select the test id's that corresopnd to current train patient
                                    targetid<-test_post %>%
                                        distinct_(.dots=patid) %>%
                                        .[which(temp$nnarraytest[1,]  == {
                                                    train_post %>%
                                                        distinct_(.dots="patient_id") %>%
                                                        .[i,] %>% unlist %>% as.vector}
                                        ),]  %>% unlist %>% as.vector
                                    targetrec<-test_post[which(test_post$patient_id %in% targetid), ] # ge tthe trainging set post op data

                                    bias<-merge(iqr,targetrec, by=time_elapsed)
                                    bias$diff<-bias$C50-bias[,outcome]

                                    #store mean of bias and rmse in biasvec, rmsevec
                                    biasvec[cnt] <- mean(bias$diff)
                                    rmsevec[cnt] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))

                                    # coverage prob based on IQR 50%
                                    bias$cov1<-bias[,outcome]-bias$C25
                                    bias$cov2<-bias$C75-bias[,outcome]
                                    bias$cov3<-ifelse(bias$cov1>0,1,0)
                                    bias$cov4<-ifelse(bias$cov2>0,1,0)
                                    bias$cov5<-bias$cov3+bias$cov4
                                    bias$coverage<-ifelse(bias$cov5>1,1,0)

                                    # store mean of the n coverage in vector
                                    coveragevec[cnt]<-mean(bias$coverage)

                                    if(any(bias$C50 > thresh_val)){
                                        crazymatch[[cnt]] <- matchmodel
                                    } else {
                                        crazymatch[[cnt]] <- NA
                                    }

                                    #-- precision
                                    precisionvec[[cnt]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 
                                    #-- store Testing C50
                                    dfList_test[[cnt]] <- test_post[which(test_post$patient_id %in% targetid), c("patient_id",time_elapsed,outcome)] %>%
                                        #train <- post[train <- post$patient <- id == ord <- data$id[c(i)],c("patient <- id",time <- elapsed,outcome)] %>% 
                                        left_join(
                                                  data.frame(time=iqr[,time_elapsed],c50 = iqr$C50) ,
                                                  by = "time"
                                        )
                                        #left_join(
                                                     #data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, patient_id = iqr[,patient_id]) ,
                                                     #by = "patient_id"
                                        #)

                                    #-- store Training C50
                                    dfList[[cnt]] <- train_post[which(train_post$patient_id %in% train[i,patid][[1]]), c("patient_id",time_elapsed,outcome)] %>%
                                        left_join(
                                                  data.frame(time=iqr[,time_elapsed],c50 = iqr$C50) ,
                                                  by = "time"
                                        )
                                        #left_join(
                                                     #data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, patient_id = iqr[,]) ,
                                                     #by = "patient_id"
                                        #)
                                    #-- store all centile for later test set merging
                                    centilepred[[cnt]] <- cbind(train[i,patid][[1]], iqr$time, iqr$C50)
                                
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

                                    iqr$iqr<-iqr$C75-iqr$C25

                                    iqrvec[cnt]<-mean(iqr$iqr)
                                    targetid<-train[i,patid][[1]]
                                    targetrec<-train_post[which(train_post$patient_id %in% targetid), ]

                                    bias<-merge(iqr,targetrec, by=time_elapsed)
                                    bias$diff<-bias$C50-bias[,outcome]

                                    #store mean of bias and rmse in biasvec, rmsevec
                                    biasvec[cnt] <- mean(bias$diff)
                                    rmsevec[cnt] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))

                                    # coverage prob based on IQR 50%
                                    bias$cov1<-bias[,outcome]-bias$C25
                                    bias$cov2<-bias$C75-bias[,outcome]
                                    bias$cov3<-ifelse(bias$cov1>0,1,0)
                                    bias$cov4<-ifelse(bias$cov2>0,1,0)
                                    bias$cov5<-bias$cov3+bias$cov4
                                    bias$coverage<-ifelse(bias$cov5>1,1,0)

                                    # store mean of the n coverage in vector
                                    coveragevec[cnt]<-mean(bias$coverage)

                                    # coverage prob based on IQR 95%
                                    bias$cov1<-bias[,outcome]-bias$C2.5
                                    bias$cov2<-bias$C97.5 -bias[,outcome]
                                    bias$cov3<-ifelse(bias$cov1>0,1,0)
                                    bias$cov4<-ifelse(bias$cov2>0,1,0)
                                    bias$cov5<-bias$cov3+bias$cov4
                                    bias$coverage<-ifelse(bias$cov5>1,1,0)

                                    # store mean of the n coverage in vector
                                    coveragevec95a[cnt]<-mean(bias$coverage)

                                    #-- precision potentially remove later because this is 7.5mb per n so 7.5*14 ~ 100MB
                                    precisionvec[[cnt]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 
                                
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
                                            biasvec = Filter(Negate(is.na),Filter(Negate(is.na), biasvec)),
                                            coveragevec = Filter(Negate(is.na),Filter(Negate(is.na), coveragevec)),
                                            centilerange = centilepred, 
                                            precisionvec=precisionvec, 
                                            zsc_list = zsc_list,
                                            rmse = rmsevec,
                                            crazymatch = crazymatch)

                        } else {
                        
                            # Get rid of NA or NULL values that are created
                            All_list <- list(
                                             #                          Filter(Negate(is.null),Filter(Negate(is.null), fin)),
                                             Filter(Negate(is.na),Filter(Negate(is.na), biasvec)),
                                             Filter(Negate(is.na),Filter(Negate(is.na), coveragevec)),
                                             Filter(Negate(is.na),Filter(Negate(is.na), coveragevec95a)),
                                             #Filter(Negate(is.na),Filter(Negate(is.na), coveragevec95)),
                                             Filter(Negate(is.na),Filter(Negate(is.na), iqrvec)),
                                             rmsevec,
                                             zsc_list,
                                             precisionvec = precisionvec,
                                             misses 
                            )

                            # name the list objects
                            message(paste0("Assigning names"))
                            #names(All <- list) <- c("bias","iqrcoverage","coverage95c","coverage95m","iqr","rmse", "dropped <- cases")
                            names(All_list) <- c("bias","iqrcoverage","coverage95c","iqr","rmse","zscore","precisionvec", "dropped_cases")

                            message(paste0("Putting in the list in array slot: ",n))
                            # store the result of the n nearest <- n result
                            nn_arr[[which(nearest == n)]] = All_list

                            # rename the list items to indicate nearestn
                            names(nn_arr)[which(nearest == n)] <- paste0('nearest_',n)


                            # Plotting
                            if (plot==TRUE){

                                # create dataframe from compiled array for plotting

                                df <- data.table(mrmse = sapply(nn_arr, function(x) {
                                                                    mean(x[['rmse']], na.rm=TRUE)
                                                            }),
                                                 mcov = sapply(nn_arr, function(x) {
                                                                   mean(x[['iqrcoverage']], na.rm=TRUE)
                                                            }),
                                                 mcov95c = sapply(nn_arr, function(x) {
                                                                      mean(x[['coverage95c']], na.rm=TRUE)
                                                            }),
                                                 mzscore = sapply(nn_arr, function(x) {
                                                                      mean(x[['zscore']], na.rm=TRUE)
                                                            }),
                                                 mmis = sapply(nn_arr, function(x) {
                                                                   mean(x[['dropped_cases']], na.rm=TRUE)
                                                            })
                                )

                                # add timepoints 
                                # df[,('nearest <- n') := nearest <- n]

                                #plot(fin[,time <- elapsed], fin$final, type="l",col="red", ylim=range(0.5:1.5), ylab="Normalized IQR", xlab="Days following Surgery", lwd=3)

                                #Add the individual data
                                #for(i in 1:nearest <- n) {
                                #lines(fin[,time <- elapsed], fin[,paste0('d',i)], col="green")
                            }

                            # select of the n result, which has the smallest mean 'rmse' or selection criteria chosen by user

                        }
            
            }

            # after loocv finishes for all nearest <- n, print which n is the lowest

            # return either LOOCV array or final prediction array
            return(nn_arr)

        }

    }


    if(!is.null(parallel)){
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Parallel solution
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # run the LOOCV function 
        pat_level_funcm <-pat_level_func

        # - - - - - - - - - - - - - - - - - - - - - - #
        # Calculation of Weighted Z-score, Coverage, and Bias to select Optimal N #
        # - - - - - - - - - - - - - - - - - - - - - - #

        if(length(nearest_n) != 1){

            if(loocv){
                loocv_test_result <- pat_level_funcm(nearest=nearest_n, loocv=loocv, parallel= parallel)
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
                    opt_n_index <- loocvperf(loocv_test_result, train %>% dplyr::select(patid) %>% rename(id = 1),bias=biasm,nearest_n) %>%
                        dplyr::select(.data$totscore) %>% unlist(.) %>% which.min(.) %>% as.vector(.)
                }
            }
            names(loocv_test_result) <- paste0("nearest_",nearest_n)

            print(paste0(mean(loocv_test_result[[opt_n_index]]$rmse, na.rm=TRUE), " from ", names(loocv_test_result))[[opt_n_index]])
            print(paste0("Number of misses is: ",loocv_test_result[[opt_n_index]]$dropped_cases,' cases'))
            print(paste0("Distribution chosen for matched GAMLSS: ", gamlss_dist))
            print(paste0("Optimal Number of Matches is: ", nearest_n[[opt_n_index]]))
            # - - - - - - - - - - - - - - - - - - - - - - #
            # Run the result on test set
            # - - - - - - - - - - - - - - - - - - - - - - #
            predict_test_result <- pat_level_funcm(nearest=nearest_n[opt_n_index], loocv=FALSE, parallel=parallel)
            # extract prediction results from the optimum nearest n 
            return(list(pred_res = predict_test_result,
                        loocv_res =  loocv_test_result,
                        loocv_score = loocvperf(loocv_test_result, train %>% dplyr::select(patid) %>% rename(id = 1), bias=biasm, nearest_n),
                        nearest_n=nearest_n[opt_n_index]))
        } else {
            if(loocv){
                loocv_test_result <- pat_level_funcm(nearest=nearest_n, loocv=loocv, parallel= parallel)
                predict_test_result <- pat_level_funcm(nearest=nearest_n, loocv=FALSE, parallel=parallel)
                retlist <- list(pred_res = predict_test_result,
                                loocv_res =  list(loocv_test_result),
                                nearest_n=nearest_n)
                names(retlist$loocv_res) <- c(paste0("nearest_",nearest_n))
                retlist$loocv_score <- loocvperf(retlist$loocv_res, train %>% dplyr::select(patid) %>% rename(id = 1), bias=biasm, nearest_n)
                return(retlist)
            } else {
                predict_test_result <- pat_level_funcm(nearest=nearest_n, loocv=FALSE, parallel=parallel)
                return(list(pred_res = predict_test_result,
                            nearest_n=nearest_n))
            }
        }

    } else {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Non Parallel solution
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # run the LOOCV function 
        pat_level_funcm <-pat_level_func

        # if we're doing this on training data
        loocv_test_result <- pat_level_funcm(nearest=nearest_n, loocv=TRUE, parallel=parallel)

        # - - - - - - - - - - - - - - - - - - - - - - #
        # Calculation of Weighted Z-score, Coverage, and Bias to select Optimal N #
        # - - - - - - - - - - - - - - - - - - - - - - #

        #opt <- n <- index <- which.min(as.vector(sapply(loocv <- test <- result, function(x) { mean(x$rmse, na.rm=TRUE)})))
        if(length(nearest_n) != 1){

            if(!is.null(userchoose)){
                opt_n_index <- which(nearest_n == userchoose)
            } else {
                opt_n_index <- loocvperf(loocv_test_result, train %>% dplyr::select(patid) %>% rename(id = 1),bias=biasm, nearest_n) %>% 
                    dplyr::select(.data$totscore) %>% unlist(.) %>% which.min(.) %>% as.vector(.)
            }

            print(paste0(mean(loocv_test_result[[opt_n_index]]$rmse, na.rm=TRUE), " from ", names(loocv_test_result))[[opt_n_index]])
            print(paste0("Number of misses is: ",loocv_test_result[[opt_n_index]]$dropped_cases,' cases'))
            print(paste0("Distribution chosen for matched GAMLSS: ", gamlss_dist))
            print(paste0("Optimal Number of Matches is: ", nearest_n[[opt_n_index]]))
            predict_test_result <- pat_level_funcm(nearest=nearest_n[opt_n_index], loocv=FALSE, parallel=parallel)
            # extract prediction results from the optimum nearest n 
            return(list(pred_res = predict_test_result,
                        loocv_res =  loocv_test_result,
                        loocv_score = loocvperf(loocv_test_result, train %>% dplyr::select(patid) %>% rename(id = 1), bias=biasm, nearest_n),
                        nearest_n=nearest_n[opt_n_index]))

        } else {

            predict_test_result <- pat_level_funcm(nearest=nearest_n, loocv=FALSE, parallel=parallel)
            return(list(pred_res = predict_test_result,
                        loocv_res =  loocv_test_result,
                        nearest_n=nearest_n))
        
        }

    }

}
