#' Generate Calibration Plot for Training or Testing Data 
#' 
#' Creates two types of plots.
#' \enumerate{
#'   \item \emph{Model performance plots} showing average bias, 
#'   coverage, 50 percent PI width (mean IQR diefference), 
#'   and a combined score of these statistics, at various 
#'   choices for the number of matches.
#'   \item \emph{Calibration plots} showing the distribution of the 
#'   observed outcomes at several predicted values. Separate plots
#'   are made for the training and test data.}
#'   
#' @param plotobj   An object produced by \code{\link{loocv_function}}
#' @param test_proc Preprocessed object from \code{\link{preproc}}
#' @param train     Logical indicating whether or not to generate
#' calibration plot for training set. Default is \code{TRUE}.
#' @param outcome   Name of the outcomes variable (type=string)
#' @param filt Logical (\code{TRUE/FALSE}) indicating whether or not to
#' filter the data in terms of performance values. This would be useful
#' if the user would want to exclude certain values in presenting the data
#' @param filter_exp  String. For filtering possible values of bias, precision, and coverage values that are out of range. (e.g. \code{"bias < 1.5"})
#' @param pred_sum  String value representing the summary used to depict
#' the predictions within the calibration. Usually \code{pred_sum = 'mean'} 
#' or \code{pred_sum = 'median'} would be a good choice to depict the 
#' summary statistic of predicted values across the deciles of observed values
#' @param obs_dist String value representing the summary used to depict
#' the observed value within the calibration plot. 
#' Usually \code{pred_sum = 'median'} woud be a good choice to depict the 
#' deciles of observed values in the calibration plot.
#' @param iqrfull Dataframe containing gamlss predictions which triggers the plotting of 
#' reference model prediction on the same plot as that of the patient like me 
#' predictions.
#' @param \dots   For specifying plotting options.
#' 
#' @return An object of class \code{ggplot} that outputs a calibration plot of observed vs. deciles of predicted values for training or testing data depending on \code{train=}.
#' 
#' @export
plot_func <- function(plotobj = plotobj, 
                      train = TRUE, 
                      test_proc = test_proc,
                      filt=filt, 
                      iqrfull=iqrfull, 
                      pred_sum = pred_sum, 
                      obs_dist=obs_dist,
                      outcome = outcome,
                      filter_exp=NULL){
  val <- "c50"
  if(train == TRUE){
    
     # For BS objects which will have null names for loocv_res
     if(is.null(names(plotobj$loocv_res))){

        nearest_n_idx <- which({unlist(lapply(plotobj$loocv_res, function(x) {
           unique(x$resdf$nearest_n)
}))} == plotobj$nearest_n)

        if(filt==TRUE){
           temp <- plotobj$loocv_res[[nearest_n_idx]]$resdf %>% filter_(filter_exp) # first one is c50
           colnames(temp)[3] <- outcome
           colnames(temp)[4] <- val
        } else {
           temp <- plotobj$loocv_res[[nearest_n_idx]]$resdf # first one is c50
           colnames(temp)[3] <- outcome
           colnames(temp)[4] <- val
        }

     } else {
        nearest_n_idx <- which(names(plotobj$loocv_res) == paste0("nearest_",plotobj$nearest_n))

        if(filt==TRUE){
           temp <- rbindlist(plotobj$loocv_res[[nearest_n_idx]]$pred_train) %>% filter_(filter_exp) # first one is c50
        } else {
           temp <- rbindlist(plotobj$loocv_res[[nearest_n_idx]]$pred_train) # first one is c50
        }

     }

    
    filtdf <- temp %>%
      bind_cols(
        dec = temp %>% 
          dplyr::select(c50) %>%
          as.data.frame() %>%
          ntile(., 10)# create deciles
      ) %>%
      left_join(
        temp %>%
          bind_cols(
            dec = temp %>% 
              dplyr::select(c50) %>%
              as.data.frame() %>%
              ntile(., 10)# create deciles
          ) %>%
          group_by(.data$dec) %>%
          summarise_(avg_val =  paste0("avg_val = ", pred_sum, "(",val,")")),
        by = "dec"
      ) 
    if(any(is.na(filtdf$dec))){
      print("NA's present in group")
      filtdf <- filtdf %>%
        filter(!is.na(.data$dec))
    }
    
    pred_tug <- temp  %>%
      bind_cols(
        dec = temp %>% 
          dplyr::select(c50) %>%
          as.data.frame() %>%
          ntile(., 10)# create deciles
      ) %>%
      group_by(.data$dec) %>%
      summarise_(avg_val = paste0("avg_val = ", pred_sum, "(",val,")"))
    
    
    #summarise(avg_val = mean(get(val))) 
    
    if(obs_dist == "gaussian" | obs_dist == "poisson") {
      print(paste0("Distribution within Each Decile of Predicted", toupper(outcome)))
      if(obs_dist == "gaussian") {
        dffitpois <- filtdf %>%
          group_by(.data$dec) %>%
          do(fitpois = glm(get(outcome) ~ 1, data= .))
        
        dfpoiscoef <- tidy(dffitpois, fitpois)
        
        # create poisson estimates and standard errors for each decile for plotting
        print(paste0("Plot DF creation"))
        plotdf <- pred_tug %>%
          left_join(dfpoiscoef %>% dplyr::select(.data$dec, .data$estimate, .data$std.error), by = "dec") %>%
          mutate(ul = .data$estimate+1.96*.data$std.error,
                 ll = .data$estimate-1.96*.data$std.error) 
      } else {
        dffitpois <- filtdf %>%
          group_by(.data$dec) %>%
          do(fitpois = glm(get(outcome) ~ 1, family=poisson(link="log"), data= .))
        
        dfpoiscoef <- tidy(dffitpois, fitpois)
        
        # create poisson estimates and standard errors for each decile for plotting
        print(paste0("Plot DF creation"))
        plotdf <- pred_tug %>%
          left_join(dfpoiscoef %>% dplyr::select(.data$dec, .data$estimate, .data$std.error), by = "dec") %>%
          mutate(
            ul = exp(.data$estimate+1.96*.data$std.error),
            ll = exp(.data$estimate-1.96*.data$std.error)) %>%
          mutate(estimate = exp(.data$estimate))
      }
      # set minimum and maximum based on the range of predicted and observed TUG values
      minc <- floor(min(plotdf$ll, plotdf$avg_val, na.rm=TRUE)) 
      maxc <- ceiling(max(plotdf$ul, plotdf$avg_val, na.rm=TRUE)) 
      
    } else { # for observed value as median and 95%CI for median
      
      dffitmed <- filtdf %>%
        group_by(.data$dec) %>%
        # group_by(.data$dec) %>%
        summarise(fitmed = MedianCI(!!dplyr::sym(outcome), conf.level= 0.95,
                                    method = "exact", R = 10000)) 
      
      dfmedcoef <- dffitmed %>%
        ungroup() %>%
        mutate(id = rep(c("median","lwr.ci","upr.ci"), length.= nrow(dffitmed))) %>%
        tidyr::pivot_wider(names_from = id, values_from = fitmed)
      
      plotdf <- pred_tug %>%
        # left_join(spread(dfmedcoef, .data$names, -.data$dec), by = "dec") %>%
        left_join(dfmedcoef, by = "dec") %>%
        rename(ul = .data$upr.ci,
               ll = .data$lwr.ci)
      # set minimum and maximum based on the range of predicted and observed TUG values
      minc <- floor(min(plotdf$ll, plotdf$avg_val, na.rm=TRUE)) 
      maxc <- ceiling(max(plotdf$ul, plotdf$avg_val, na.rm=TRUE)) 
      #minc <- floor(min(plotdf[,val], plotdf[,outcome], na.rm=TRUE)) 
      #maxc <- ceiling(max(boxplot.stats(plotdf[,val])$stats[5], boxplot.stats(plotdf[,outcome])$stats[5], na.rm=TRUE))*1.05
    }
    
    
    
    # plot observed vs. predicted TUG on decile of predicted TUG
    if(obs_dist == "gaussian" | obs_dist == "poisson"){
      if(obs_dist == "gaussian") {
        cptrain <-  ggplot(plotdf, aes(x = .data$avg_val, y = .data$estimate, 
                                       ymin = .data$estimate-1.96*.data$std.error, 
                                       ymax=.data$estimate+1.96*.data$std.error,
                                       colour="PLM"
        )) + 
          geom_pointrange() + 
          #xlim(minc, maxc) + ylim(minc,maxc) + 
          geom_abline(slope=1, intercept=0) + 
          xlab(paste0("Predicted ",toupper(outcome))) + 
          ylab(paste0("Observed ",toupper(outcome))) 
      } else {
        cptrain <-  ggplot(plotdf, aes(x = .data$avg_val, 
                                       y = .data$estimate, 
                                       ymin = .data$ll, 
                                       colour="PLM",
                                       ymax=.data$ul)) + 
          geom_pointrange() + 
          #xlim(minc, maxc) + ylim(minc,maxc) + 
          geom_abline(slope=1, intercept=0) + 
          xlab(paste0("Predicted ",toupper(outcome))) + 
          ylab(paste0("Observed ",toupper(outcome))) 
      }
      
      
    } else {
      # median and IQR
      cptrain <-  ggplot(plotdf, aes(x = .data$avg_val,
                                     ymin = .data$ll, ymax=.data$ul,
                                     colour="PLM",
                                     y = .data$median
      )) + geom_pointrange() + 
        #xlim(minc, maxc) + ylim(minc,maxc) + 
        geom_abline(slope=1, intercept=0) + 
        xlab(paste0("Predicted ",toupper(outcome))) + 
        ylab(paste0("Observed ",toupper(outcome))) 
      
      #cptrain <-  ggplot(plotdf, aes_string(x = "avg_val",
      #y = outcome,
      #group = "avg_val"
      #)) + geom_boxplot(outlier.colour=NA) + 
      ##xlim(minc, maxc) + ylim(minc,maxc) + 
      #geom_abline(slope=1, intercept=0) + 
      #xlab("Predicted TUG") + 
      #ylab("Observed TUG")
    }
    cptrainzsc <-  ggplot(plotdf, aes_string(x = "avg_val",
                                             y = "zsc",
                                             colour="PLM",
                                             group = "avg_val"
    )) + geom_boxplot(outlier.colour=NA) + 
      #xlim(minc, maxc) + ylim(minc,maxc) + 
      geom_hline(yintercept=0) + 
      xlab(paste0("Predicted ",toupper(outcome))) + 
      ylab("Z-score")
    
    mincz <- floor(min(plotdf$zsc, na.rm=TRUE)) 
    maxcz <- ceiling(max(plotdf$zsc, na.rm=TRUE)) 
    
    # correlation
    traincor <- cor(temp %>% 
                      dplyr::select(!!dplyr::sym(outcome), c50) %>%
                      as.matrix, method="spearman") %>%
      .[1,2] %>%
      round(3) %>%
      paste0("r^2 ==", .)
    
    # other performance measures
    
    #zscres <- loocvperf(plotobj$loocv_res, train_o = test_proc$train_o, bias="zsc", nearest_n=plotobj$loocv_score$nearest_n)
    #test_perf <- zscres[which.min(zscres$totscore),c("zscore","coverage","precision")]
   test_perf <- plotobj$loocv_score[plotobj$loocv_score$nearest_n == plotobj$nearest_n, c("rmse","cov","prec")]
    
    if(!is.null(iqrfull)){
      
      if(any(colnames(iqrfull) %in% "C50")){
        
        # Create Reference Plot 
        cptrainrefdf <- test_proc$train_post %>%
          dplyr::select(!!dplyr::sym(test_proc$varname[3]), !!dplyr::sym(test_proc$varname[1]), !!dplyr::sym(test_proc$varname[2])) %>% 
          left_join(
            iqrfull ,
            by="time"
          ) %>%
          mutate(dec= ntile(.data$C50, 10)) %>%
          group_by(.data$dec) %>%
          summarise_(avg_val =  paste0("avg_val = ", pred_sum, "(C50)")) %>%
          left_join(
            test_proc$train_post %>%
              dplyr::select(!!dplyr::sym(test_proc$varname[3]), !!dplyr::sym(test_proc$varname[1]), !!dplyr::sym(test_proc$varname[2])) %>% 
              left_join(
                iqrfull ,
                by="time"
              ) %>%
              mutate(dec= ntile(.data$C50, 10)
              ) %>% left_join(
                test_proc$train_post %>%
                  dplyr::select(!!dplyr::sym(test_proc$varname[3]), !!dplyr::sym(test_proc$varname[1]), !!dplyr::sym(test_proc$varname[2])) %>% 
                  left_join(
                    iqrfull ,
                    by="time"
                  ) %>%
                  mutate(dec= ntile(.data$C50, 10)) %>%
                  group_by(.data$dec) %>%
                  summarise_(avg_val =  paste0("avg_val = ", pred_sum, "(C50)"))
              ) %>%
              group_by(.data$dec) %>%
              do(fitmed = MedianCI(.[,outcome][[1]], conf.level = 0.95, 
                                   method = "exact", R = 10000, na.rm = TRUE)) %>%
              tidy(., fitmed) %>%
              spread(., .data$names, -.data$dec)
          ) %>%
          rename(ul = .data$upr.ci,
                 ll = .data$lwr.ci) 
        
        #cptrainrefplot <- cptrainrefdf %>%
        #ggplot(., aes(x = .data$avg_val,
        #ymin = .data$ll, ymax=.data$ul,
        #y = .data$median,
        #)) + geom_pointrange() + 
        #geom_abline(slope=1, intercept=0) + 
        #xlab(paste0("Predicted ",toupper("HDC"))) + 
        #ylab(paste0("Observed ",toupper("HDC"))) +
        #xlim(c(36, 50)) + ylim(c(36,50)) +
        #theme_bw()  + theme(legend.position="none", aspect.ratio = 1) 
        
        mincr <- floor(min(cptrainrefdf$ll, cptrainrefdf$avg_val, na.rm=TRUE)) 
        maxcr <- ceiling(max(cptrainrefdf$ul, cptrainrefdf$avg_val, na.rm=TRUE)) 
        
        return(list(cptrain + theme_bw()+ geom_pointrange(aes(x=cptrainrefdf$avg_val,
                                                              ymin = cptrainrefdf$ll, ymax=cptrainrefdf$ul,
                                                              colour = "REF",
                                                              y=cptrainrefdf$median)), 
                    min(minc,mincr), max(maxc, maxcr), 
                    traincor, 
                    cptrainzsc + theme_bw(), 
                    mincz, maxcz,tp=test_perf
                    #mincz, maxcz,tp=test_perf, 
                    #cptrainrefplot + theme_bw() 
        ))
        
      } else {
        stop("Reference Object Not Right Class for Creating Reference Prediction Chart. Require gamlss or lme4 obj")
      }
      
    } else {
      return(list(cptrain + theme_bw(), minc, maxc, traincor, cptrainzsc + theme_bw(), mincz, maxcz,tp=test_perf ))
    }
    #return(list(cptrain, minc, maxc))
    
  } else { # for test set we need to merge from the full data
    
    #-- add the z-score dataset 
    #tempzsc <- rbindlist(lapply(plotobj$pred_res$zsc_list, list))[[1]][!is.na(rbindlist(lapply(plotobj$pred_res$zsc_list, list))) & !is.infinite(unlist(rbindlist(lapply(plotobj$pred_res$zsc_list, list))))]
    #tempzsc <- rbindlist(lapply(plotobj$pred_res$zsc_list, function(x) {list( zsc = x)}), idcol = "patient_id")
    
    #-- get id of train_df that has closests y90 value to test df
    print("creating temp and test matching data")
    
    #temp <-  test_proc$test_post %>% 
      #dplyr::select_("patient_id", "time", outcome) %>%
      #left_join(
        #test_proc$test_o %>%
          #bind_cols(
            #train_id = test_proc$train_o$id[sapply(1:length(test_proc$test_o$pred), function(x) {
              #which.min(abs(test_proc$test_o$pred[x] - test_proc$train_o$Fitted))                                         })]
          #) %>% 
          #dplyr::select(.data$id, .data$train_id) %>%
          #rename(patient_id = .data$id),
        #by = "patient_id"
      #) %>%
      #rename(test_id = .data$patient_id)  %>%
      ##-- join centiles.pred from the unique train id's that were selected based on minimum difference between y90 of train and test
      #left_join(
        #data.frame(train_id = unique(test_proc$train_o$id[sapply(1:length(test_proc$test_o$pred), function(x) {
          #which.min(abs(test_proc$test_o$pred[x] - test_proc$train_o$Fitted))
        #})]) ) %>% 
          #full_join(
            #rbindlist(lapply(plotobj$pred_res$centilerange, data.frame)) %>%
              #rename(train_id = 1,
                     #time = 2,
                     #C50 = 3),
            #by = "train_id"
          #)
        #,
        #by = c("time","train_id") 
      #) 

  if(is.null(plotobj$pred_res$pred_test)) {

     # For BS
     if(class(plotobj$pred_res[[1]]) %in% c("numeric", "integer")){
        temp <- plotobj$pred_res %>% 
           dplyr::select(1:4) %>%
           rename(C50 = 4)

        colnames(temp)[3] <- outcome

     } else {
     # for pmmsknn 
        temp <- data.table::rbindlist(plotobj$pred_res[[1]]$pred_test) %>%
           dplyr::select(1:4) %>%
           rename(C50 = 4)
     }
  } else {
     temp <- data.table::rbindlist(plotobj$pred_res$pred_test) %>%
        dplyr::select(1:4) %>%
        rename(C50 = 4)
  }

    #-- bind with decile
    print("binding with decile")
    print("filt df made in test")
    filtdf <- temp %>%
      bind_cols(
        dec = temp %>%
          dplyr::select(.data$C50) %>%
          as.data.frame() %>%
          ntile(., 10)
      ) %>%
      left_join(
        temp %>%
          bind_cols(
            dec = temp %>%
              dplyr::select(.data$C50) %>%
              as.data.frame() %>%
              ntile(., 10)
          ) %>%
          group_by(.data$dec) %>%
          summarise_(avg_val = paste0("avg_val = ", pred_sum, "(C50)")),
        by = "dec"
      ) 
    if(any(is.na(filtdf$dec))){
      warning("NA's present in group")
      filtdf <- filtdf %>%
        dplyr::filter(!is.na(.data$dec))
    }
    
    print(paste0("Predicting TUG values"))
    
    pred_tug <-  temp %>%
      bind_cols(
        dec = temp %>%
          dplyr::select(.data$C50) %>%
          as.data.frame() %>%
          ntile(., 10)
      ) %>%
      group_by(.data$dec) %>%
      summarise_(avg_val = paste0("avg_val = ", pred_sum, "(C50)"))
    
    
    # use normal, poisson, or median/95% IQR
    if(obs_dist == "gaussian" | obs_dist == "poisson") {
      print(paste0("Distribution within Each Decile of Predicted TUG"))
      if(obs_dist == "gaussian"){
        dffitpois <- filtdf %>%
          group_by(.data$dec) %>%
          do(fitpois = glm(get(outcome)~ 1, data= .))
        dfpoiscoef <- tidy(dffitpois, fitpois)
        
        # create poisson estimates and standard errors for each decile for plotting
        print(paste0("Plot DF creation"))
        plotdf <- pred_tug %>%
          left_join(dfpoiscoef %>% dplyr::select(.data$dec, .data$estimate, .data$std.error), by = "dec") %>%
          mutate(ul = .data$estimate+1.96*.data$std.error,
                 ll = .data$estimate-1.96*.data$std.error) 
      } else {
        
        dffitpois <- filtdf %>%
          group_by(.data$dec) %>%
          do(fitpois = glm(get(outcome)~ 1,family=poisson(link="log"), data= .))
        
        dfpoiscoef <- tidy(dffitpois, fitpois)
        
        # create poisson estimates and standard errors for each decile for plotting
        print(paste0("Plot DF creation"))
        plotdf <- pred_tug %>%
          left_join(dfpoiscoef %>% dplyr::select(.data$dec, .data$estimate, .data$std.error), by = "dec") %>%
          mutate(ul = exp(.data$estimate+1.96*.data$std.error),
                 ll = exp(.data$estimate-1.96*.data$std.error)) %>%
          mutate(estimate = exp(.data$estimate))
      }
      
      
      # set minimum and maximum based on the range of predicted and observed TUG values
      minc <- floor(min(plotdf$ll, plotdf$avg_val, na.rm=TRUE)) 
      maxc <- ceiling(max(plotdf$ul, plotdf$avg_val, na.rm=TRUE)) 
      
    } else {
      
      dffitmed <- filtdf %>%
        group_by(.data$dec) %>%
        # group_by(.data$dec) %>%
        summarise(fitmed = MedianCI(!!dplyr::sym(outcome), conf.level= 0.95,
                                    method = "exact", R = 10000)) 
      
      dfmedcoef <- dffitmed %>%
        ungroup() %>%
        mutate(id = rep(c("median","lwr.ci","upr.ci"), length.= nrow(dffitmed))) %>%
        tidyr::pivot_wider(names_from = id, values_from = fitmed)
      
      plotdf <- pred_tug %>%
        # left_join(spread(dfmedcoef, .data$names, -.data$dec), by = "dec") %>%
        left_join(dfmedcoef, by = "dec") %>%
        rename(ul = .data$upr.ci,
               ll = .data$lwr.ci)
      # set minimum and maximum based on the range of predicted and observed TUG values
      minc <- floor(min(plotdf$ll, plotdf$avg_val, na.rm=TRUE)) 
      maxc <- ceiling(max(plotdf$ul, plotdf$avg_val, na.rm=TRUE)) 
      #plotdf <- filtdf
      
      ## set minimum and maximum based on the range of predicted and observed TUG values
      #minc <- floor(min(plotdf$C50, unlist(plotdf[,outcome]), na.rm=TRUE)) 
      #maxc <- ceiling(max(boxplot.stats(plotdf$C50)$stats[5], boxplot.stats(unlist(plotdf[,outcome]))$stats[5], na.rm=TRUE))*1.05
    }
    
    
    # plot observed vs. predicted TUG on decile of predicted TUG
    if(obs_dist == "gaussian" | obs_dist == "poisson"){
      if(obs_dist == "gaussian"){
        cptest <-  ggplot(plotdf, aes(x = .data$avg_val, 
                                      y = .data$estimate, 
                                      ymin = .data$estimate-1.96*.data$std.error, 
                                      colour = "PLM",
                                      ymax=.data$estimate+1.96*.data$std.error)) + 
          geom_pointrange() + 
          #xlim(minc, maxc) + ylim(minc,maxc) + 
          geom_abline(slope=1, intercept=0) + 
          xlab(paste0("Predicted ",toupper(outcome))) + 
          ylab(paste0("Observed ",toupper(outcome))) 
        
      } else {
        cptest <-  ggplot(plotdf, aes(x = .data$avg_val, 
                                      y = .data$estimate, 
                                      ymin = .data$ll, 
                                      colour = "PLM",
                                      ymax=.data$ul)) + 
          geom_pointrange() + 
          #xlim(minc, maxc) + ylim(minc,maxc) + 
          geom_abline(slope=1, intercept=0) + 
          xlab(paste0("Predicted ",toupper(outcome))) + 
          ylab(paste0("Observed ",toupper(outcome))) 
        
      }
      
    } else {
      
      cptest <-  ggplot(plotdf, aes(x = .data$avg_val,
                                    ymin = .data$ll, ymax=.data$ul,
                                    colour = "PLM",
                                    y = .data$median
                                    #y = estimate
      )) + 
        geom_pointrange() + 
        #xlim(minc, maxc) + ylim(minc,maxc) + 
        geom_abline(slope=1, intercept=0) + 
        xlab(paste0("Predicted ",toupper(outcome))) + 
        ylab(paste0("Observed ",toupper(outcome))) 
      
      #cptest <-  ggplot(plotdf, aes_string(x = "avg_val",
      #y = outcome,
      #group = "avg_val"
      #)) + 
      #geom_boxplot(outlier.colour=NA) + 
      ##xlim(minc, maxc) + ylim(minc,maxc) + 
      #geom_abline(slope=1, intercept=0) + 
      #xlab("Predicted TUG") + 
      #ylab("Observed TUG")
      
    }
    cptestzsc <-  ggplot(plotdf, aes_string(x = "avg_val",
                                            y = "zsc",
                                            colour = "PLM",
                                            group = "avg_val"
    )) + geom_boxplot(outlier.colour=NA) + 
      #xlim(minc, maxc) + ylim(minc,maxc) + 
      geom_hline(yintercept=0) + 
      xlab(paste0("Predicted ",toupper(outcome))) + 
      ylab("Z-score")
    
    mincz <- floor(min(plotdf$zsc, na.rm=TRUE)) 
    maxcz <- ceiling(max(plotdf$zsc, na.rm=TRUE)) 
    
    # correlation
    testcor <- cor(temp %>% 
                     dplyr::select_(outcome, "C50") %>%
                     as.matrix, method="spearman", use="complete.obs") %>%
      .[1,2] %>%
      round(3) %>%
      paste0("r^2 ==", .)
    
    # other performance measures
    #test_perf <- extvalid(plotobj, test_proc)
    test_perf <- plotobj$test_score
    
    # Reference plot plotting
    if(!is.null(iqrfull)){
      
      if(any(colnames(iqrfull) %in% "C50")){
        
        print("creating reference plot")
        # Create Reference Plot 
        cptestrefdf <- test_proc$test_post %>%
          dplyr::select(!!dplyr::sym(test_proc$varname[3]), !!dplyr::sym(test_proc$varname[1]), !!dplyr::sym(test_proc$varname[2])) %>% 
          left_join(
            iqrfull ,
            by="time"
          ) %>%
          mutate(dec= ntile(.data$C50, 10)) %>%
          group_by(.data$dec) %>%
          summarise_(avg_val =  paste0("avg_val = ", pred_sum, "(C50)")) %>%
          left_join(
            test_proc$test_post %>%
              dplyr::select(!!dplyr::sym(test_proc$varname[3]), !!dplyr::sym(test_proc$varname[1]), !!dplyr::sym(test_proc$varname[2])) %>% 
              left_join(
                iqrfull ,
                by="time"
              ) %>%
              mutate(dec= ntile(.data$C50, 10)
              ) %>% left_join(
                test_proc$test_post %>%
                  dplyr::select(!!dplyr::sym(test_proc$varname[3]), !!dplyr::sym(test_proc$varname[1]), !!dplyr::sym(test_proc$varname[2])) %>% 
                  left_join(
                    iqrfull ,
                    by="time"
                  ) %>%
                  mutate(dec= ntile(.data$C50, 10)) %>%
                  group_by(.data$dec) %>%
                  summarise_(avg_val =  paste0("avg_val = ", pred_sum, "(C50)"))
              ) %>%
              group_by(.data$dec) %>%
              do(fitmed = MedianCI(.[,outcome][[1]], conf.level = 0.95, 
                                   method = "exact", R = 10000, na.rm = TRUE)) %>%
              tidy(., fitmed) %>%
              spread(., .data$names, -.data$dec),
            by = "dec"
          ) %>%
          rename(ul = .data$upr.ci,
                 ll = .data$lwr.ci) 
        
        #cptestrefplot <- cptestrefdf %>%
        #ggplot(., aes(x = .data$avg_val,
        #ymin = .data$ll, ymax=.data$ul,
        #colour = "REF",
        #y = .data$median,
        #)) + geom_pointrange() + 
        ##xlim(minc, maxc) + ylim(minc,maxc) + 
        #geom_abline(slope=1, intercept=0) + 
        #xlab(paste0("Predicted ",toupper("HDC"))) + 
        #ylab(paste0("Observed ",toupper("HDC"))) +
        #xlim(c(36, 50)) + ylim(c(36,50)) +
        #theme_bw()  + theme(legend.position="none", aspect.ratio = 1) 
        
        
        mincr <- floor(min(cptestrefdf$ll, cptestrefdf$avg_val, na.rm=TRUE)) 
        maxcr <- ceiling(max(cptestrefdf$ul, cptestrefdf$avg_val, na.rm=TRUE)) 
        
        return(
          list(cptest + theme_bw() + 
                 geom_pointrange(
                   aes(x=cptestrefdf$avg_val,
                       ymin = cptestrefdf$ll, 
                       ymax=cptestrefdf$ul,
                       y=cptestrefdf$median, 
                       colour="REF")
                 )
               , 
               min(minc,mincr), max(maxc, maxcr), 
               testcor, 
               cptestzsc + theme_bw(), 
               mincz, maxcz,tp=test_perf
               #mincz, maxcz,tp=test_perf, 
               #cptestrefplot + theme_bw() 
          ))
      } else {
        stop("Reference Object Not Right Class for Creating Reference Prediction Chart. Require gamlss or lme4 obj")
      }
      
    } else {
      return(list(cptest + theme_bw(), 
                  minc, maxc, testcor, 
                  cptestzsc + theme_bw(), 
                  mincz, maxcz, 
                  tp=test_perf
      )
      )
    }
    
  }
}
