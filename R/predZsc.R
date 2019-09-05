#' Predict Z-Score for Training and Testing Data
#' 
#' The function that predicts Z-score values initial GAMLSS distribution to the full data
#' in order to be used as a reference distribution for the following
#' patient like me GAMLSS fits
#' 
#' @param plmr  Fitted object from \code{\link{plm}}
#' @param traintestmatchdf  Matched train test dataframe based on \code{\link{matchTestDataGen}}
#' @param matchmodel - dataset generated from the \code{\link{matchTrainDataGen}} 
#' @param ord_data Idem, component \code{train_o}
#' @param train_post - datasets, typically the \code{train_post} list component 
#'  of the object produced by \code{\link{preproc}}.
#' @param test_post Idem, component \code{test_post}
#' @param loocv - Logical indicating whether or not to predict for train or test set.
#' @param time_elapsed  Name of the time variable. (type=string)
#' @param outcome       Name of the outcomes variable
#' @param i             Patient id indicator.
#' \code{TRUE} will result int zscores from the training set, otherwise zscores from test set
#' will be generated.
#' 
#' @return There are many possible return values 
#' 
#' @export
predZsc <- function(
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
                 ){
    if(!loocv){

        testpred=data.frame()
        # - - - - - - - - - - - - - - - - - - - - - - #
        # Iterate through each testing obs(x) closest to the current train obs(i) 
        # get predicted values using centles.pred()
        # - - - - - - - - - - - - - - - - - - - - - - #
        for(x in traintestmatchdf[traintestmatchdf$train_id %in% ord_data$id[c(i)], "test_id"]){
            #message(paste0("PREDICTING FOR TEST = ",x))
            testpred <- testpred %>% 
                bind_rows(
                          data.frame(
                                     zsc = centiles.pred(plmr, type="z-scores",
                                                         xname=time_elapsed,
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

        #out <- tryCatch(
                        #{
                            #message("ZSCORE PREDICTION TEST PART")
                            #message(paste0("PREDICTING FOR TRAIN = ",i))
                            #testpred=data.frame()
                            ## - - - - - - - - - - - - - - - - - - - - - - #
                            ## Iterate through each testing obs(x) closest to the current train obs(i) 
                            ## get predicted values using centles.pred()
                            ## - - - - - - - - - - - - - - - - - - - - - - #
                            #for(x in traintestmatchdf[traintestmatchdf$train_id %in% ord_data$id[c(i)], "test_id"]){
                                #message(paste0("PREDICTING FOR TEST = ",x))
                                #testpred <- testpred %>% 
                                    #bind_rows(
                                              #data.frame(
                                                         #zsc = centiles.pred(plmr, type="z-scores",
                                                                             #xname=time_elapsed,
                                                                             #data=matchmodel,
                                                                             #xvalues=test_post %>%
                                                                                 #filter(.data$patient_id %in% x) %>%
                                                                                 #dplyr::select(.data$time) %>%
                                                                                 #unlist %>% as.vector,
                                                                             #yval=test_post %>%
                                                                                 #filter(.data$patient_id %in% x) %>%
                                                                                 #dplyr::select_(outcome) %>%
                                                                                 #unlist %>% as.vector
                                                                             #),
                                                         #test_id = rep(x, length(test_post %>%
                                                                                 #filter(.data$patient_id %in% x) %>%
                                                                                 #dplyr::select(.data$time) %>%
                                                                                 #unlist %>% as.vector
                                                                             #)),
                                                         #time = test_post %>%
                                                             #filter(.data$patient_id %in% x) %>%
                                                             #dplyr::select(.data$time) %>%
                                                             #unlist %>% as.vector
                                                         #))
                            #}

                            #return(testpred)
                        #},
                        #error = function(cond)
                        #{
                            #message(paste("ERROR in ZSCORE PREDICTION"))
                            ##message("ORIGINAL ERROR:")
                            #message(cond)
                            #return(NA)
                        #},
                        #warning=function(cond)
                        #{
                            #message(paste("GAMLSS WARNING"))
                            ##message("ORIIGNAL WARNING:")
                            #message(cond)

                            #testpred=data.frame()
                            ## - - - - - - - - - - - - - - - - - - - - - - #
                            ## Iterate through each testing obs(x) closest to the current train obs(i) 
                            ## get predicted values using centles.pred()
                            ## - - - - - - - - - - - - - - - - - - - - - - #
                            #for(x in traintestmatchdf[traintestmatchdf$train_id %in% ord_data$id[c(i)], "test_id"]){
                                #message(paste0("PREDICTING FOR TEST = ",x))
                                #testpred <- testpred %>% 
                                    #bind_rows(
                                              #data.frame(
                                                         #zsc = centiles.pred(plmr, type="z-scores",
                                                                             #xname="time",
                                                                             #data=matchmodel,
                                                                             #xvalues=test_post %>%
                                                                                 #filter(.data$patient_id %in% x) %>%
                                                                                 #dplyr::select(.data$time) %>%
                                                                                 #unlist %>% as.vector,
                                                                             #yval=test_post %>%
                                                                                 #filter(.data$patient_id %in% x) %>%
                                                                                 #dplyr::select_(outcome) %>%
                                                                                 #unlist %>% as.vector
                                                                             #),
                                                         #test_id = rep(x, length(test_post %>%
                                                                                 #filter(.data$patient_id %in% x) %>%
                                                                                 #dplyr::select(.data$time) %>%
                                                                                 #unlist %>% as.vector
                                                                             #)),
                                                         #time = test_post %>%
                                                             #filter(.data$patient_id %in% x) %>%
                                                             #dplyr::select(.data$time) %>%
                                                             #unlist %>% as.vector
                                                         #))
                            #}
                            #return(testpred)
                            ##return(NULL)
                        #},
                        #finally={
                            ##message("PROCESSING GAMLSS")
                            ##message("GAMLSS EOL SUCCESS")
                        #}
        #)

        #return(out)

    } else {

        trainzsc <- data.frame(
                               zsc = centiles.pred(plmr, type="z-scores",
                                                   xname=time_elapsed,
                                                   data=matchmodel,
                                                   xvalues=train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]],
                                                   yval=train_post[train_post$patient_id %in% ord_data$id[c(i)],outcome][[1]]),
                               train_id = rep(ord_data$id[c(i)], length(train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]]))
                               ,
                               time = train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]]
        )
        return(trainzsc)

        #zscf <- function(){
            #out <- tryCatch(
                            #{
                                ##message("ZSCORE PREDICTION LOOCV TRY PART")
                                #trainzsc <- data.frame(
                                                       #zsc = centiles.pred(plmr, type="z-scores",
                                                                           #xname=time_elapsed,
                                                                           #data=matchmodel,
                                                                           #xvalues=train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]],
                                                                           #yval=train_post[train_post$patient_id %in% ord_data$id[c(i)],outcome][[1]]),
                                                       #train_id = rep(ord_data$id[c(i)], length(train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]]))
                                                       #,
                                                       #time = train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]]
                                #)
                                #return(trainzsc)

                            #},
                            #error = function(cond)
                            #{
                                #message(paste("ERROR in ZSCORE PREDICTION"))
                                ##message("ORIGINAL ERROR:")
                                #message(cond)
                                #return(NA)
                            #},
                            #warning=function(cond)
                            #{
                                #message(paste("GAMLSS WARNING"))
                                ##message("ORIIGNAL WARNING:")
                                #message(cond)
                                #trainzsc <- data.frame(
                                                       #zsc = centiles.pred(plmr, type="z-scores",
                                                                           #xname=time_elapsed,
                                                                           #data=matchmodel,
                                                                           #xvalues=train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]],
                                                                           #yval=train_post[train_post$patient_id %in% ord_data$id[c(i)],outcome][[1]]),
                                                       #train_id = rep(ord_data$id[c(i)], length(train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]]))
                                                       #,
                                                       #time = train_post[train_post$patient_id %in% ord_data$id[c(i)],time_elapsed][[1]]
                                #)
                                #return(trainzsc)
                                ##return(NULL)
                            #},
                            #finally={
                                ##message("PROCESSING GAMLSS")
                                ##message("GAMLSS EOL SUCCESS")
                            #}
            #)
            #return(out)
        #}
    
    }
}
