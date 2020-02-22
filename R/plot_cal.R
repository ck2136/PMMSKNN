#' Plot for statistical validity and calibration
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
#' @param plotobj   - An object produced by \code{\link{loocv_function}}
#' @param test_proc - Preprocessed object from \code{\link{preproc}}
#' @param outcome   - Name of the outcomes variable (type=string)
#' @param filt Logical (\code{TRUE/FALSE}) indicating whether or not to
#' filter the data in terms of performance values. This would be useful
#' if the user would want to exclude certain values in presenting the data
#' @param pred_sum  - String value representing the summary used to depict
#' the predictions within the calibration. Usually \code{pred_sum = 'mean'} 
#' or \code{pred_sum = 'median'} would be a good choice to depict the 
#' summary statistic of predicted values across the deciles of observed values
#' @param obs_dist - String value representing the summary used to depict
#' the observed value within the calibration plot. 
#' Usually \code{pred_sum = 'median'} woud be a good choice to depict the 
#' deciles of observed values in the calibration plot.
#' @param loocv Logical indicating the type of plot: 
#' Model performance plot (if \code{loocv = TRUE}, default), 
#' or or calibration plot (if \code{loocv = FALSE}). 
#' @param filter_exp - String. For filtering possible values of bias, precision, and coverage values that are out of range. (e.g. \code{"bias < 1.5"})
#' @param plot_cal_zscore - Logical (\code{TRUE/FALSE}) indicating whether to plot zscore calibration 
#' @param wtotplot - Logical (\code{TRUE/FALSE}) indicating wehter to include a weighted total score plot
#' that indicates the optimal n match based on equally weighting bias, coverage and precision
#' @param plotvals - Logical (\code{TRUE/FALSE}) indicating whether to plot bias, coverage, and precision values onto the calibration plot
#' @param iqrfull - Dataframe containing gamlss predictions which triggers the plotting of reference model prediction on the same plot as that of the patient like me predictions.
#' @param \dots   - For specifying plotting options.
#' 
#' @return An object of class \code{ggplot} that outputs a calibration plot of observed vs. deciles of predicted values.
#' 
#' @export
plot_cal <- function(plotobj,
                     test_proc=test_proc,
                     outcome = "tug",
                     filt=FALSE,
                     pred_sum="mean",
                     obs_dist="median",
                     #plot_by=seq(10,150,5),
                     loocv=TRUE,
                     filter_exp = NULL,
                     plot_cal_zscore=FALSE,
                     wtotplot=FALSE,
                     plotvals=FALSE,
                     iqrfull=NULL,
                     bs=FALSE,
                     ...) {
  
  # - - - - - - - - - - - - - - - - - - - - - - #
  # Instantiate all plot objects for viewing
  # - - - - - - - - - - - - - - - - - - - - - - #
  
  #-- NON CALIBRATION plots only
  if(loocv){
    
    # - - - - - - - - - - - - - - - - - - - - - - #
    # RMSE/Coverage Plot function for test and train
    # - - - - - - - - - - - - - - - - - - - - - - #
    
    # For brokenstick object it's simple data manipulation of loocv_score
    if(bs){
      tmp1 <- plotobj$loocv_score %>%
        tidyr::pivot_longer(rmse:prec, names_to = "measure") %>% rename(nearest_n = 1)
      perfdf <- plotobj$loocv_score
    
    } else {
      nearest_n =as.numeric(regmatches(names(plotobj$loocv_res), regexpr("\\d+",names(plotobj$loocv_res)))) 
      
      perfdf <- loocv_perf(
        plotobj$loocv_res,
        outcome=outcome,
        nearest_n=nearest_n,
        perf_round_by=4
      ) 
      
      tmp1 <- perfdf %>%
        tidyr::pivot_longer(rmse:prec, names_to = "measure") %>% rename(nearest_n = 1)
    }
    
    # tmp1 <-listtodf(plotobj$loocv_res) 
    
    train_bias <- ggplot(tmp1 %>%
                           filter(
                             #abs(value) < 50,
                             #measure == 'bias' | measure == 'rmse' | measure == 'zscore') 
                             #measure == 'bias' | measure == 'zscore') 
                             .data$measure == 'rmse') 
    ) + 
      xlab("Matches (N)") + ylab("RMSE") +
      geom_point(aes(x=.data$nearest_n, y=.data$value, colour=.data$measure)) + 
      #geom_smooth(aes(x=.data$nearest_n, y=.data$value, colour = .data$measure), 
                  #method="gam",formula = y ~ s(x, bs="cs", k=splinek ), se=FALSE) + 
      theme_bw()  + theme(legend.position="none", aspect.ratio = 1) +
      geom_hline(yintercept = 0) 
      #annotate("text", x = median(tmp1$nearest_n), y = 0, vjust = -1, label = "0 Bias")
    
    # Coverage (Excluding extreme measures)
    train_cov <- ggplot(tmp1 %>%
                          filter(
                            #nearest_n > 10,
                            .data$measure == 'cov') 
                        #measure == 'iqrcoverage' | measure == 'coverage95c' ) 
    ) + 
      geom_point(aes(x=.data$nearest_n, y=.data$value), colour="blue") + 
      #geom_smooth(aes(x=.data$nearest_n, y=.data$value), colour = "blue", 
                  #method="gam",formula = y ~ s(x, bs="cs", k=splinek ), se=FALSE) + 
                       xlab("Matches (N)") + ylab("Coverage (50%)") +
      ylim(min(tmp1 %>% 
                 filter(.data$measure == 'cov') %>%
                 dplyr::select(.data$value) %>%
                 unlist %>%
                 as.vector) * 0.95 ,
           max(tmp1 %>%
                 filter(.data$measure == 'cov') %>%
                 dplyr::select(.data$value) %>%
                 unlist %>%
                 as.vector) * 1.05)+
      #ylim(0.3,1)+
      #scale_colour_manual(labels=c("95% IQR Coverage","50% IQR Coverage"), values=c("blue","red")) +
      scale_colour_manual(labels=c("50% IQR difference"), values=c("blue")) +
      theme_bw() + theme(legend.position="none", aspect.ratio = 1) +
      geom_hline(yintercept = 0.50) 
      #annotate("text", x = median(tmp1$nearest_n),y = .50, vjust = -0.1, label = "Coverage")
    
    # - - - - - - - - - - - - - - - - - - - - - - #
    # Precision Plot: Mean IQR dif by Nearest N
    # - - - - - - - - - - - - - - - - - - - - - - #
    pppm <- ggplot(tmp1 %>%
                     filter(
                       #nearest_n > 10,
                       .data$measure == 'prec') 
                   #measure == 'iqrcoverage' | measure == 'coverage95c' ) 
    ) +
      geom_point(aes(x=.data$nearest_n, y=.data$value), colour="green") + 
      #geom_smooth(aes(x=.data$nearest_n, y=.data$meaniqrdif), colour="green",
                 #method="gam",formula = y ~ s(x, bs="cs", k=splinek ), se=FALSE) + 
      xlab("Matches (N)") + ylab("Mean IQR difference") +
      ylim(min(tmp1 %>% 
                 filter(.data$measure == 'prec') %>%
                 dplyr::select(.data$value) %>%
                 unlist %>%
                 as.vector) * 0.95 ,
           max(tmp1 %>% 
                 filter(.data$measure == 'prec') %>%
                 dplyr::select(.data$value) %>%
                 unlist %>%
                 as.vector) * 1.05)+
      #ylim(0.3,1)+
      #scale_colour_manual(labels=c("95% IQR Coverage","50% IQR Coverage"), values=c("blue","red")) +
      #scale_colour_manual(labels=c("50% IQR Coverage"), values=c("blue")) +
      theme_bw() + theme(legend.position="none", aspect.ratio = 1) +
      geom_hline(yintercept = max(perfdf$prec)) 
      #annotate("text", x=median(ppdf_means$nearest_n), y = max(ppdf_means %>% dplyr::select(.data$meaniqrdif) %>% unlist %>% as.vector), vjust = -1, label = "Max IQR Difference")
    
    
  # - - - - - - - - - - - - - - - - - - - - - - #
  # Return plot objects
  # - - - - - - - - - - - - - - - - - - - - - - #
  if(wtotplot){
      # - - - - - - - - - - - - - - - - - - - - - - #
      # Weighted Total Score Plot included
      # - - - - - - - - - - - - - - - - - - - - - - #
      wtspdf <- loocvperf(plotobj$loocv_res, test_proc$train_o)

      wtsp <- ggplot(wtspdf) + 
          geom_point(aes(x=.data$nearest_n, y=.data$totscore)) + 
          #geom_smooth(aes(x=.data$nearest_n, y=.data$totscore), method="gam",
                 #formula = y ~ s(x, bs="cs", k=splinek ), se=FALSE) + 
          xlab("Matches (N)") + ylab("Weighted Total Score") +
          theme_bw() + theme(legend.position="none", aspect.ratio = 1)

      return(plot_grid(train_bias, train_cov, pppm, wtsp, labels="AUTO", ncol=2))
  } else {
      # - - - - - - - - - - - - - - - - - - - - - - #
      # Weighted Total Score Plot Not included
      # - - - - - - - - - - - - - - - - - - - - - - #
      return(plot_grid(train_bias, train_cov, pppm,  labels="AUTO", ncol=3))
  }
    
  } else {
    
    print("creating training calibration plot")
    cptrainlist = plot_func(plotobj = plotobj, 
                            test_proc = test_proc,
                            train=TRUE, 
                            filt=filt, 
                            iqrfull=iqrfull, 
                            pred_sum=pred_sum, 
                            obs_dist=obs_dist,
                           outcome=outcome
    )
    print("creating testing calibration plot")
    cptestlist = plot_func(plotobj = plotobj, 
                           train=FALSE, 
                           test_proc = test_proc,
                           filt=filt, 
                           iqrfull=iqrfull, 
                           pred_sum=pred_sum, 
                           obs_dist=obs_dist, 
                           outcome=outcome
    )
    
    if(plot_cal_zscore==FALSE){
      
      minc <- floor(min(cptrainlist[[2]], cptestlist[[2]], na.rm=TRUE))
      maxc <- ceiling(max(cptrainlist[[3]], cptestlist[[3]], na.rm=TRUE))
      
      # PLOT BIAS, PRECISION, COVERAGE
      if(plotvals){
          #labels
          train_zs_lab <- paste0("zscore == ", round(cptrainlist$tp$zscore,3))
          train_cov_lab <- paste0("coverage == ", round(cptrainlist$tp$coverage,3))
          train_prec_lab <- paste0("precision == ", round(cptrainlist$tp$precision,3))

          test_zs_lab <- paste0("zscore == ", round(cptestlist$tp$zscore, 3))
          test_cov_lab <- paste0("coverage == ", round(cptestlist$tp$coverage, 3))
          test_prec_lab <- paste0("precision == ", round(cptestlist$tp$precision,3))

          # Calibration plots
          cptrain <- cptrainlist[[1]] + xlim(minc, maxc) + ylim(minc,maxc) + 
              #geom_text(aes(label=paste0(cptrainlist[[4]]), y=minc+(maxc-minc)/10+1, x=(minc+maxc)*0.6), parse= TRUE, color="red") +
              geom_text(aes(label=paste0(train_zs_lab), y=minc+(maxc-minc)/10+(maxc-minc)/10, x=(minc+maxc)*0.6), parse= TRUE, color="red") +
              geom_text(aes(label=paste0(train_cov_lab), y=minc+(maxc-minc)/10, x=(minc+maxc)*0.6), parse= TRUE, color="blue") +
              geom_text(aes(label=paste0(train_prec_lab), y=minc+(maxc-minc)/10-(maxc-minc)/10, x=(minc+maxc)*0.6), parse= TRUE, color="green")

          cptest <- cptestlist[[1]] + xlim(minc, maxc) + ylim(minc,maxc) + 
              #geom_text(aes(label=paste0(cptestlist[[4]]), y=minc+(maxc-minc)/10+1, x=(minc+maxc)*0.6), parse= TRUE, color="red")+
              geom_text(aes(label=paste0(test_zs_lab), y=minc+(maxc-minc)/10+(maxc-minc)/10, x=(minc+maxc)*0.6), parse= TRUE, color="red") +
              geom_text(aes(label=paste0(test_cov_lab), y=minc+(maxc-minc)/10, x=(minc+maxc)*0.6), parse= TRUE, color="blue") +
              geom_text(aes(label=paste0(test_prec_lab), y=minc+(maxc-minc)/10-(maxc-minc)/10, x=(minc+maxc)*0.6), parse= TRUE, color="green")

          return(plot_grid(cptrain + theme(aspect.ratio = 1), cptest + theme(aspect.ratio = 1), labels = "AUTO", vjust = 3))
      } else {

          minc <- floor(min(cptrainlist[[2]], cptestlist[[2]], na.rm=TRUE))
          maxc <- ceiling(max(cptrainlist[[3]], cptestlist[[3]], na.rm=TRUE))
          cptrain <- cptrainlist[[1]] + 
              xlim(minc, maxc) + ylim(minc,maxc) + 
              theme(aspect.ratio = 1) 
              #scale_colour_manual(name="", values=c("REF"="red", "PLM"="blue"),
                                  #guide = guide_legend(fill = NULL, colour=NULL))

          cptest <- cptestlist[[1]] + 
              xlim(minc, maxc) + ylim(minc,maxc) + 
              theme(aspect.ratio = 1) 
              #scale_colour_manual(name="", values=c("REF"="red", "PLM"="blue"),
                                  #guide = guide_legend(fill = NULL, colour=NULL))

          cpfin <- plot_grid(cptrain  + theme(legend.position = "none"), 
                           cptest + theme(legend.position = "none"), 
                           align = "vh"
                           )

          #legend <- get_legend(cptrain)
          return(plot_grid(cpfin))
          #return(plot_grid(cpfin, legend, rel_widths=c(3,0.3)))
          #if(!is.null(iqrfull)){
              #cptrainref <- cptrainlist[[9]] + xlim(minc, maxc) + ylim(minc,maxc) 
              #cptestref <- cptestlist[[9]] + xlim(minc, maxc) + ylim(minc,maxc) 
              #return(plot_grid(cptrain + theme(aspect.ratio = 1), 
                               #cptest + theme(aspect.ratio = 1),
                               #cptrainref + theme(aspect.ratio = 1),
                               #cptestref + theme(aspect.ratio = 1),
                               #labels = "AUTO", label_y = 3, ncol=2))
          #} else {
              #return(plot_grid(cptrain + theme(aspect.ratio = 1), cptest + theme(aspect.ratio = 1), labels = "AUTO", label_y = 3, ncol=2))
          #}
      }
      
      
    } else {
      
      minc <- floor(min(cptrainlist[[6]], cptestlist[[6]], na.rm=TRUE))
      maxc <- ceiling(max(cptrainlist[[7]], cptestlist[[7]], na.rm=TRUE))
      
      # Calibration plots
      cptrain <- cptrainlist[[5]] + 
          xlim(minc, maxc) + ylim(minc,maxc) + 
          theme(aspect.ratio = 1) 
          #scale_colour_manual(name="", values=c("REF"="red", "PLM"="blue"),
                              #guide = guide_legend(fill = NULL, colour=NULL))

      cptest <- cptestlist[[5]] + 
          xlim(minc, maxc) + ylim(minc,maxc) + 
          theme(aspect.ratio =1) 
          #scale_colour_manual(name="", values=c("REF"="red", "PLM"="blue"),
                              #guide = guide_legend(fill = NULL, colour=NULL))


      #legend <- get_legend(cptrain)
      cpfin <- plot_grid(cptrain + theme(legend.position = "none"), 
                       cptest + theme(legend.position = "none"), 
                       align="vh"
      )

      return(plot_grid(cpfin))
      #return(plot_grid(cpfin, legend, rel_widths=c(3,0.3)))
      
      
    }
    
  }
  #return(plot_grid(cptrain, cptest, train_bias, train_cov, ppp, pppm, labels = "AUTO"))
  
}
