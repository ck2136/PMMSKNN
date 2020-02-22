#'@importFrom brokenstick brokenstick get_knots
#'@importFrom reshape2    melt
#'@importFrom broom       tidy
#'@importFrom gamlss      gamlss centiles.pred cs pb gamlss.control
#'@importFrom gamlss.dist BCCGo BCPEo BCTo GA NO
#'@importFrom dplyr       bind_rows filter filter_ mutate %>% group_by
#'                        n select select_ left_join distinct distinct_
#'                        summarise summarise_ arrange arrange_ 
#'                        rename rename_ bind_cols full_join ntile do 
#'                        if_else row_number group_by_ slice ungroup
#'                        tibble mutate_ top_n
#'@importFrom data.table  data.table rbindlist 
#'@importFrom cowplot     plot_grid get_legend
#'@importFrom DescTools   MedianCI
#'@importFrom lme4        lmer
#'@importFrom merTools    predictInterval
#'@importFrom tidyr       spread nest
#'@importFrom tidyselect  contains matches
#'@importFrom rlang       .data := parse_quo caller_env global_env
#'@importFrom ggplot2     aes aes_string facet_wrap geom_abline 
#'                        geom_boxplot geom_hline geom_line geom_point
#'                        geom_pointrange geom_text ggplot geom_smooth
#'                        scale_colour_manual theme theme_bw
#'                        xlab xlim ylab ylim annotate guide_legend
#'                        labs
#'@importFrom MASS        stepAIC mvrnorm
#'@importFrom doParallel  registerDoParallel stopImplicitCluster
#'@importFrom doSNOW      registerDoSNOW 
#'@importFrom parallel    makeCluster detectCores 
#'@importFrom future      plan multiprocess
#'@importFrom future.apply      future_lapply
#'@importFrom foreach     foreach %dopar%
#'@importFrom stats       as.formula complete.cases cor formula
#'                        glm lm median  na.omit pnorm poisson 
#'                        predict setNames update median rchisq
#'                        model.matrix coef quantile fitted qt sd
#'@importFrom utils       head str globalVariables
NULL

utils::globalVariables(c(".","fitmed", "fitpois","test_id","train_id","..","bs_obj"))

#' \pkg{PMMSKNN}: Sequential KNN Extended via Predicted Mean Matching.
#'
#' This package provides functions that allow personalized predictions of longitudinal outcome trajectory based on extending the sequential KNN algorithm via Predicted Mean Matching. Essentially individuals are matched according to a distal outcome of interest based on characteristics that explains the variation in the outcome. The trajectory of the matched population is used to generate predictions for each individuals.
#'
#' @section PMMSKNN functions:
#' The main functions are:
#' \tabular{ll}{
#'   \code{preproc()}        \tab Preprocess the data\cr
#'   \code{loocv_function()} \tab Evaluate statistical properties\cr
#'   \code{plot_cal()}       \tab Create calibration plots}
#' @docType package
#' @name PMMSKNN-pkg
#' @seealso \code{\link{brokenstick}}
#' @note
#' Development of this package was kindly supported by 
#' FUNDER under THIS GRANT NUMBER
#' @references
#' Kittelson et al. (2019). \emph{Development and testing of a 
#' neighbors-based prediction for physical function after total 
#' knee arthroplasty}. In preparation.
NULL
