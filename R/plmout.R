#' Obtaining output of updated patient level model predictions 
#' 
#' The function primarily takes the output of the updated reference gamlss()
#' and outputs results of the performance measures based on a fitted reference model (This happens within the \code{\link{pat_level_func}}).
#' 
#' @param plmr          Fitted gamlss object updated with patient level data
#' based on \code{\link{fitrefgamlss}} and \code{update()}. 
#' @param i             Patient id indicator (type=numeric).
#' @param time_window   vector of numbers for `centiles.pred()`, `xvalues` argument 
#' @param mint          Numeric value indicating minimum time for predicting values
#' @param maxt          Numeric value indicating maximum time for predicting values
#' @param perfout       List of performance measures to be filled. 
#' @param loocv - Logical (\code{TRUE/FALSE}) that specifies whether 
#' or not to perform leave-one-out cross validation or just output 
#' predictions without hyperparameter tuning. If \code{loocv=FALSE}, then
#' users need to specify the value of the nearest_n 
#' @param traintestmatchdf  Matched train test dataframe based on \code{\link{matchTestDataGen}}
#' @param matchmodel    Dataset generated from the \code{\link{matchTrainDataGen}} 
#' @param time_elapsed  Name of the time variable. (type=string)
#' @param ord_data Data frame. Specifically, training data with patient_id ordered based on fitted distal outcome value using predicted mean matching.
#' Generated using \code{\link{preproc}}. Example, \code{x <- preproc()}, 
#' then \code{x$train_o} would be used for this parameter.
#' @param train_post - Data frame that contains the post-baseline observations from the training dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param test_post Data frame that contains the post-baseline observations from the testing dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param outcome    - Name of the outcomes variable (type=string)
#' @param thresh_val - Numeric value indicating value of bias to ignore (not include in output) in terms of the leave-one-out cross validation process. The default is set to \code{thresh_val = 10000}
#' 
#' @return An array of performance measures for each individuals in the training and testing set 
#' 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# Function that outputs performance measaures for individual predictions ----
# - - - - - - - - - - - - - - - - - - - -#
plmout <- function(
                   plmr, # updated patient level model
                   perfout=perfout,
                   i=i, time_window=time_window, mint,maxt,
                   loocv=loocv,matchmodel=matchmodel,
                   traintestmatchdf=traintestmatchdf,
                   outcome=outcome, time_elapsed=time_elapsed, 
                   ord_data=ord_data, train_post=train_post, 
                   test_post=test_post, thresh_val = thresh_val
                   ) {

    # - - - - - - - - - - - - - - - - - - - - - - #
    # 1. GET PREDICTIONS from:
    # plmr(), time_elapsed, mint,maxt,matchmodel
    # - - - - - - - - - - - - - - - - - - - - - - #
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
    # TEST DATA PREDICTION: 
    # - - - - - - - - - - - - - - - - - - - - - - #
    if(loocv!=TRUE){

        # - - - Zscore 
        zsc <- predZsc(
                       plmr=plmr,
                       traintestmatchdf=traintestmatchdf,
                       matchmodel=matchmodel,
                       ord_data=ord_data,
                       train_post=train_post, 
                       test_post=test_post, 
                       loocv=loocv,
                       time_elapsed=time_elapsed,
                       outcome = outcome,
                       i=i
        )

        if(!is.data.frame(zsc)) {
            message('Something wrong with Z-SCORE prediction. Prediction not included')
            perfout$zscore[[i]] <- data.frame(zsc = NA,
                                          test_id = traintestmatchdf[traintestmatchdf$train_id %in% ord_data$id[c(i)], "test_id"],
                                          time = NA
            )
        } else {
            perfout$zscore[[i]] <- zsc
        }

        # -- IQR values
        # iqr$iqr<-iqr$C75-iqr$C25
        # perfout$iqrvec[i]<-mean(iqr$iqr)

        targetid<-traintestmatchdf[traintestmatchdf$train_id %in% ord_data$id[c(i)], "test_id"] # select the test id's that corresopnd to current train patient
        # targetrec<-test_post[which(test_post$patient_id %in% targetid), ] # ge tthe trainging set post op data

        # bias<-merge(iqr,targetrec, by=time_elapsed)
        # bias$diff<-bias$C50-bias[,outcome]

        # Bais RAW
        # perfout$biasvec[i] <- mean(bias$diff)
        # Bais RMSE
        # perfout$rmsevec[i] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))

        # 50% Coverage
        # bias$cov1<-bias[,outcome]-bias$C25
        # bias$cov2<-bias$C75-bias[,outcome]
        # bias$cov3<-ifelse(bias$cov1>0,1,0)
        # bias$cov4<-ifelse(bias$cov2>0,1,0)
        # bias$cov5<-bias$cov3+bias$cov4
        # bias$coverage<-ifelse(bias$cov5>1,1,0)

        # Mean of 50% Coverage
        # perfout$coveragevec[i]<-mean(bias$coverage)


        # Crazy Matches 
        if(any(iqr$C50 > thresh_val)){
            perfout$crazymatch[[i]] <- matchmodel
        } else {
            perfout$crazymatch[[i]] <- NA
        }

        #-- Precision
        # perfout$precisionvec[[i]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 
        #-- Store the Test Predicted Values (i.e.C50)
        perfout$dfList_test[[i]] <- test_post[which(test_post$patient_id %in% targetid), c("patient_id",time_elapsed,outcome)] %>%
            left_join(
                      data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                      by = "time"
            )

        #-- Store the Train Predicted Values
        perfout$dfList[[i]] <- train_post[which(train_post$patient_id %in% ord_data$id[c(i)]), c("patient_id",time_elapsed,outcome)] %>%
            left_join(
                      data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25 = iqr$C25, c75 = iqr$C75) ,
                      by = "time"
            )
        #-- All centile 
        # perfout$centilepred[[i]] <- cbind(ord_data$id[c(i)], iqr$time, iqr$C50)


        #-- Final Output
        return(
               perfout
        )

    } else { # loocv = TRUE

        # - - - -
        # z-scores vector
        # - - - - 

        zsc <- predZsc(
                       plmr=plmr,
                       traintestmatchdf=traintestmatchdf,
                       matchmodel=matchmodel,
                       ord_data=ord_data,
                       train_post=train_post, 
                       test_post=test_post, 
                       loocv=loocv,
                       time_elapsed=time_elapsed,
                       outcome = outcome,
                       i=i
        )
        if(!is.data.frame(zsc)) {

            message('Something wrong with Z-SCORE prediction. Prediction not included')
            perfout$zscore[[i]] <- data.frame(zsc = NA,
                                          train_id = ord_data$id[c(i)],
                                          time = NA
            )

        } else {
            perfout$zscore[[i]] <-  zsc
        }

        #message(paste0("Getting MLE  predictions: "))

        # - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - #
        # code for storing either LOOCV result or extracting prediction #
        # - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - #

        # iqr$iqr<-iqr$C75-iqr$C25

        # perfout$iqrvec[i]<-mean(iqr$iqr)
        targetid<-ord_data$id[c(i)]
        # targetrec<-train_post[which(train_post$patient_id %in% targetid), ]

        # bias<-merge(iqr,targetrec, by=time_elapsed)
        # bias$diff<-bias$C50-bias[,outcome]

        #store mean of bias and rmse in biasvec, rmsevec
        # perfout$biasvec[i] <- mean(bias$diff)
        # perfout$rmsevec[i] <- sqrt(sum(na.omit(bias$diff)^2)/length(na.omit(bias$diff)))

        # coverage prob based on IQR 50%
        # bias$cov1<-bias[,outcome]-bias$C25
        # bias$cov2<-bias$C75-bias[,outcome]
        # bias$cov3<-ifelse(bias$cov1>0,1,0)
        # bias$cov4<-ifelse(bias$cov2>0,1,0)
        # bias$cov5<-bias$cov3+bias$cov4
        # bias$coverage<-ifelse(bias$cov5>1,1,0)

        # store mean of the n coverage in vector
        # perfout$coveragevec[i]<-mean(bias$coverage)

        # coverage prob based on IQR 95%
        # bias$cov1<-bias[,outcome]-bias$C2.5
        # bias$cov2<-bias$C97.5 -bias[,outcome]
        # bias$cov3<-ifelse(bias$cov1>0,1,0)
        # bias$cov4<-ifelse(bias$cov2>0,1,0)
        # bias$cov5<-bias$cov3+bias$cov4
        # bias$coverage<-ifelse(bias$cov5>1,1,0)

        # store mean of the n coverage in vector
        # perfout$coveragevec95[i]<-mean(bias$coverage)
        #-- Test Predicted Values (i.e.C50)
        # perfout$dfList_test[[i]] <- test_post[which(test_post$patient_id %in% targetid), c("patient_id",time_elapsed,outcome)] %>%
        #     left_join(
        #               data.frame(time=iqr[,time_elapsed],c50 = iqr$C50) ,
        #               by = "time"
        #     )

        #-- Train Predicted Values
        perfout$dfList[[i]] <- train_post[which(train_post$patient_id %in% ord_data$id[c(i)]), c("patient_id",time_elapsed,outcome)] %>%
            left_join(
                      data.frame(time=iqr[,time_elapsed],c50 = iqr$C50, c25=iqr$C25, c75=iqr$C75) ,
                      by = "time"
            )

        #-- precision potentially remove later because this is 7.5mb per n so 7.5*14 ~ 100MB
        # perfout$precisionvec[[i]] <-list(time= iqr[,time_elapsed], prec=iqr$iqr) 


        return(
               perfout
        )
    }

}
