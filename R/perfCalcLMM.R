#' Calculate Performance (bias, coverage, and precision) based on fitting a linear Mixed Model to longitudinal data
#' 
#' Creates Data Frame:
#' \enumerate{
#'   \item Bias (Difference in Obs vs Pred)
#'   \item Coverage (% time obs value in Pred Interval)
#'   \item Precision (Width of Pred Interval)}
#' @param fit       An object produced by \code{lmer} function.
#' @param n         Number of simulations for obtaining prediction 
#' and prediction Interval. Default = 1000.
#' @param train     A dataframe of training observations with \code{outcome} column
#' included.
#' @param test      A dataframe of testing observations with \code{outcome} column
#' included.
#' @param parallel  Boolean value (i.e. \code{TRUE/FALSE}) indicating whether or not to use parallel computing.
#' Default is set to FALSE. If \code{true}, then uses multicore capability.
#' @param outcome   A column name (i.e. string type) indicating \code{outcome} column
#' 
#' @return          A list of data frames. 1) \code{perfsum}: Data frame that summarizes the performance based on a variety of variance/confidence interval estimation of LMM; 2) \code{predictInterval1, predictInterval1a, predictInterval1b}: List of two data frames that contain the prediction (95% CI) for training and testing data; 3) \code{bootMerpar, bootMersp}: List of two data frames that contain the prediction (95% CI) for training and testing data using bootMer package
#' 
#' @export
perfCalcLMM <- function(
        fit, # lmer model from training
        n=1000,
        train, # training data
        test, # test data
        parallel=FALSE,
        outcome # outcome name
){

    if(parallel){
        par = "multicore"
    } else {
        par = "no"
    }
    # OPTION 1  predictInterval()
    PItrain1 <- predictInterval(merMod = fit, train,
                                level = 0.5, n.sims = n,
                                stat = "median", type = "linear.prediction",
                                include.resid.var = TRUE, .parallel = parallel)
    PItest1 <- predictInterval(merMod = fit, test,
                               level = 0.5, n.sims = n,
                               stat = "median", type = "linear.prediction", 
                               include.resid.var = TRUE, .parallel = parallel)

    # LESS CONSERVATIVE OPTION
    PItrain1a <- predictInterval(merMod = fit, train,
                                level = 0.5, n.sims = n,
                                stat = "median", type = "linear.prediction",
                                include.resid.var = FALSE, .parallel = parallel)
    PItest1a <- predictInterval(merMod = fit, test,
                               level = 0.5, n.sims = n,
                               stat = "median", type = "linear.prediction", 
                               include.resid.var = FALSE, .parallel = parallel)

    PItrain1b <- predictInterval(merMod = fit, train,
                                level = 0.5, n.sims = n,
                                stat = "median", type = "linear.prediction",
                                include.resid.var = FALSE, 
                                fix.intercept.variance=TRUE,
                                .parallel = parallel)
    PItest1b <- predictInterval(merMod = fit, test,
                               level = 0.5, n.sims = n,
                               stat = "median", type = "linear.prediction", 
                               fix.intercept.variance=TRUE,
                               include.resid.var = FALSE, .parallel = parallel)

    # Option 2 bootMer()
    mySumm1 <- function(.) {
        predict(., newdata=train, re.form=NULL)
    }
    mySumm2 <- function(.) {
        predict(., newdata=test, re.form=NULL, allow.new.levels=TRUE)
    }
    ####Collapse bootstrap into median, 50% PI
    sumBoot <- function(merBoot) {
        return(
               data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
                          lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.25, na.rm=TRUE))),
                          upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.75, na.rm=TRUE)))
               )
        )
    }

    ##lme4::bootMer() method 1
    boot1 <- lme4::bootMer(fit, mySumm1, nsim=n, use.u=FALSE, type="parametric", parallel = par)
    boot2 <- lme4::bootMer(fit, mySumm1, nsim=n, use.u=TRUE, type="semiparametric", parallel = par)

    PItrain2 <- sumBoot(boot1)
    PItrain3 <- sumBoot(boot2)

    boot1 <- lme4::bootMer(fit, mySumm2, nsim=n, use.u=FALSE, type="parametric", parallel = par)
    boot2 <- lme4::bootMer(fit, mySumm2, nsim=n, use.u=TRUE, type="semiparametric", parallel = par)

    PItest2 <- sumBoot(boot1)
    PItest3 <- sumBoot(boot2)

    # Final Return with bot plots
    return(
           list(
                perfsum = data.frame(
                                     Method = c("LMM: PI-a","LMM: PI-b","LMM: PI-c","LMM: bootMer-par","LMM: bootMer-semi"),
                                     train_bias = c(
                                                    mean(abs(train[,outcome][[1]] - PItrain1$fit)),
                                                    mean(abs(train[,outcome][[1]] - PItrain1a$fit)),
                                                    mean(abs(train[,outcome][[1]] - PItrain1b$fit)),
                                                    mean(abs(train[,outcome][[1]] - PItrain2$fit)),# boot-parametric
                                                    mean(abs(train[,outcome][[1]] - PItrain3$fit))# boot-semi
                                                    ),
                                     train_cov = c(
                                                   mean(ifelse(ifelse((PItrain1$upr - train[,outcome][[1]]) > 0, 1, 0) +  
                                                               ifelse((train[,outcome][[1]] - PItrain1$lwr) > 0, 1, 0) > 1, 1, 0)),
                                                   mean(ifelse(ifelse((PItrain1a$upr - train[,outcome][[1]]) > 0, 1, 0) +  
                                                               ifelse((train[,outcome][[1]] - PItrain1a$lwr) > 0, 1, 0) > 1, 1, 0)),
                                                   mean(ifelse(ifelse((PItrain1b$upr - train[,outcome][[1]]) > 0, 1, 0) +  
                                                               ifelse((train[,outcome][[1]] - PItrain1b$lwr) > 0, 1, 0) > 1, 1, 0)),
                                                   mean(ifelse(ifelse((PItrain2$upr - train[,outcome][[1]]) > 0, 1, 0) +  
                                                               ifelse((train[,outcome][[1]] - PItrain2$lwr) > 0, 1, 0) > 1, 1, 0)),
                                                   mean(ifelse(ifelse((PItrain3$upr - train[,outcome][[1]]) > 0, 1, 0) +  
                                                               ifelse((train[,outcome][[1]] - PItrain3$lwr) > 0, 1, 0) > 1, 1, 0))
                                                   ),
                                     train_prec = c(
                                                    mean(PItrain1$upr - PItrain1$lwr),
                                                    mean(PItrain1a$upr - PItrain1a$lwr),
                                                    mean(PItrain1b$upr - PItrain1b$lwr),
                                                    mean(PItrain2$upr - PItrain2$lwr),
                                                    mean(PItrain3$upr - PItrain3$lwr)
                                                    ),
                                     test_bias = c(
                                                   mean(abs(test[,outcome][[1]] - PItest1$fit)),
                                                   mean(abs(test[,outcome][[1]] - PItest1a$fit)),
                                                   mean(abs(test[,outcome][[1]] - PItest1b$fit)),
                                                   mean(abs(test[,outcome][[1]] - PItest2$fit)),# boot-parametric
                                                   mean(abs(test[,outcome][[1]] - PItest3$fit))# boot-semi
                                                   ),
                                     test_cov = c(
                                                  mean(ifelse(ifelse((PItest1$upr - test[,outcome][[1]]) > 0, 1, 0) +  
                                                              ifelse((test[,outcome][[1]] - PItest1$lwr) > 0, 1, 0) > 1, 1, 0)),
                                                  mean(ifelse(ifelse((PItest1a$upr - test[,outcome][[1]]) > 0, 1, 0) +  
                                                              ifelse((test[,outcome][[1]] - PItest1a$lwr) > 0, 1, 0) > 1, 1, 0)),
                                                  mean(ifelse(ifelse((PItest1b$upr - test[,outcome][[1]]) > 0, 1, 0) +  
                                                              ifelse((test[,outcome][[1]] - PItest1b$lwr) > 0, 1, 0) > 1, 1, 0)),
                                                  mean(ifelse(ifelse((PItest2$upr - test[,outcome][[1]]) > 0, 1, 0) +  
                                                              ifelse((test[,outcome][[1]] - PItest2$lwr) > 0, 1, 0) > 1, 1, 0)),
                                                  mean(ifelse(ifelse((PItest3$upr - test[,outcome][[1]]) > 0, 1, 0) +  
                                                              ifelse((test[,outcome][[1]] - PItest3$lwr) > 0, 1, 0) > 1, 1, 0))
                                                  ),
                                     test_prec = c(
                                                   mean(PItest1$upr - PItest1$lwr),
                                                   mean(PItest1a$upr - PItest1a$lwr),
                                                   mean(PItest1b$upr - PItest1b$lwr),
                                                   mean(PItest2$upr - PItest2$lwr),
                                                   mean(PItest3$upr - PItest3$lwr)
                                     )
                ),
                predictInterval1 = list(train = PItrain1, test = PItest1),
                predictInterval1a = list(train = PItrain1a, test = PItest1b),
                predictInterval1b = list(train = PItrain1a, test = PItest1b),
                bootMerpar = list(train = PItrain2, test = PItest2),
                bootMersp = list(train = PItrain3, test = PItest3) 
           )
           
    )
}
