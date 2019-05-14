#' Plot for calibration of LMM models
#' 
#' Creates Calibration Plots:
#' \enumerate{
#'   \item \emph{Calibration plots} showing the distribution of the 
#'   observed outcomes at several predicted values. Separate plots
#'   are made for the training and test data.}
#' @param fit       An object produced by \code{lmer} function.
#' @param train     A dataframe of training observations with \code{outcome} column
#' included.
#' @param test      A dataframe of testing observations with \code{outcome} column
#' included.
#' @param outcome   A column name (i.e. string type) indicating \code{outcome} column
#' @param perfObj   An object from \code{\link{perfCalcLMM}}
#' @param fitopt    Number indicating which LMM prediction method to use (mostly for CI).
#' 1 = predictInterval() a; 2 = predictInterval() b; 3 = predictInterval() c; 
#' 4 = bootMer() parametric; 5 = bootMer() semiparametric
#' @param \dots     Used for specifying options in the predict() method
#' @return A list with components ..., or an object of class 
#' \code{ggplot} [??]
#' @export
plot_cal_lmer <- function(
        fit, # lmer model from training
        train, # training data
        test, # test data
        outcome,
        perfObj = NULL,
        fitopt = 1,
        ...
){

    if(is.null(perfObj)){
        # - - - - - - - - - - - - - - - - - - - - #
        # training portion
        # - - - - - - - - - - - - - - - - - - - - #
        pred <- predict(fit)

        plottrain <- calplotDFgen(
                                  train[,outcome][[1]],
                                  pred,
                                  n = 10
        )

        # - - - - - - - - - - - - - - - - - - - - #
        # TESTING PORTION
        # - - - - - - - - - - - - - - - - - - - - #
        pred <- predict(fit, test, allow.new.levels = TRUE)

        plottest <- calplotDFgen(
                                 test[,outcome][[1]],
                                 pred,
                                 n = 10
        )
    } else {
        if(fitopt < 1 | fitopt > 5) {
            stop('fit opt needs to be greater than 1 or less than 6')
        }
        pred <- perfObj[[fitopt+1]]$train$fit
        plottrain <- calplotDFgen(
                                  train[,outcome][[1]],
                                  pred,
                                  n = 10
        )
        pred <- perfObj[[fitopt+1]]$test$fit
        plottest <- calplotDFgen(
                                 test[,outcome][[1]],
                                 pred,
                                 n = 10
        )
    }

    # - - - - - - - - - - - - - - - - - - - - #
    # PLOT GENERATION
    # - - - - - - - - - - - - - - - - - - - - #

    # SET THE MIN AND MAX for X AND Y BASED ON THE TEST AND TRAIN PLOTS
    minc <- floor(min(plottrain[[2]], plottest[[2]], na.rm=TRUE))
    maxc <- ceiling(max(plottrain[[3]], plottest[[3]], na.rm=TRUE))

    # GENERATE PLOTS
    cptrain <-  ggplot(plottrain[[1]], aes(x = .data$avg_val,
                                   ymin = .data$ll, ymax=.data$ul,
                                   colour="PLM",
                                   y = .data$median
                                   )) + geom_pointrange() + 
                  #xlim(minc, maxc) + ylim(minc,maxc) + 
                  geom_abline(slope=1, intercept=0) + 
                  xlim(minc, maxc) + ylim(minc,maxc) + 
                  xlab(paste0("Predicted ",toupper(outcome))) + 
                  ylab(paste0("Observed ",toupper(outcome)))  + theme_bw()

    # CALIBRATION PLOT FOR TESTING
    cptest <-  ggplot(plottest[[1]], aes(x = .data$avg_val,
                                   ymin = .data$ll, ymax=.data$ul,
                                   colour="PLM",
                                   y = .data$median
                                   )) + geom_pointrange() + 
                  #xlim(minc, maxc) + ylim(minc,maxc) + 
                          xlim(minc, maxc) + ylim(minc,maxc) + 
                  geom_abline(slope=1, intercept=0) + 
                  xlab(paste0("Predicted ",toupper(outcome))) + 
                  ylab(paste0("Observed ",toupper(outcome)))   + theme_bw()

    # Final Return with bot plots
    return(
           plot_grid(cptrain + theme(aspect.ratio = 1, legend.position="none"), 
                     cptest + theme(aspect.ratio = 1, legend.position="none"), 
                     vjust = 3)
    )
}
