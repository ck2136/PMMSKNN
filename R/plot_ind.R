#' Plot an individuals outcome trajectory using PMMSKNN-GAMLSS
#' 
#' Creates a single plot with the below specification:
#' \enumerate{
#'   \item \emph{Individuals data points} showing individual patient values 
#'   \item \emph{Matched indivduals' data points} showing the 
#'   data points from the nearest neighbors of the selected individual
#'   }
#'   
#' @param test_proc   Preprocessed object from \code{\link{preproc}}
#' @param outcome     Name of the outcomes variable (type=string)
#' @param idnum       Integer indicating the id of the individual
#' @param mtype       Integer value indicating matching type. Default is set to 1 which follows the
#' matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. \code{mtype} values are 
#' from \code{0} to \code{4}
#' @param n           Integer representing number of matches.
#' @param dist        Object of \code{\link{gamlss.dist}} specifying the distribution for gamlss
#' @param df_m        Numeric value that specifies the degrees of freedom for the cubic spline specified for the mean parameter of the distribution specified according to \code{dist_fam}
#' @param df_s        Numeric value that specifies the degrees of freedom for the cubic spline specified for the scale parameter of the distribution specified according to \code{dist_fam}
#' @param df_n        Numeric value that specifies the degrees of freedom for the cubic spline specified for the shape parameter, specifically the \eqn{\nu} parameter, of the distribution specified according to \code{dist_fam}
#' @param df_t        Numeric value that specifies the degrees of freedom for the cubic spline specified for the shape parameter, specifically the \eqn{\tau} parameter, of the distribution specified according to \code{dist_fam}
#' @param xvalues     Vector of values for the prediction centile to be plotted (default: 3:200)
#' @param xlab        String indicating x axis label for plot (default: "Time")
#' @param ylab        String indicating y axis label for plot (default: "Outcome")
#' @param \dots       Options to specify in the plotting of the percentile curves
#' 
#' @return            A list that contains 1) an object of class \code{ggplot} that outputs a gamlss predition curve for an individual in terms of their predicted values over time, 2) a data frame that contains the observations from the matched neighbors as well as the individual w idnum and 3) A data frame of centiles prediction results  (10th, 25th, 50th, 75th, and 90th)
#' 
#' @export
plot_ind <- function(
                     test_proc=test_proc,
                     outcome = "tug",
                     idnum = 1,
                     mtype = 1,
                     n=10,
                     dist=BCCGo,
                     df_m=2,
                     df_s=1,
                     df_n=1,
                     df_t=1,
                     xvalues=3:200,
                     xlab = "Time",
                     ylab = "Outcome",
                     ...) {
  # Main input: object from LOOCV function which spits out train data output
  # plotobj$pred_res$pred contains the training set predictions in dataframe
  # plotobj$pred_res

    #-- DATASET for patient with idnum
    matches <- matchIdExtractTest(test_proc = test_proc, 
                                  mtype = mtype, 
                                  i = idnum,
                                  n = n) 

    matchmodel <- test_proc$train_post %>% 
        filter(.data$patient_id %in% matches)

    plotdf <- test_proc$test_post %>%
        filter(.data$patient_id == idnum) %>%
        bind_rows(
                  matchmodel
        )

    #-- GAMLSS FITTING and CENTILE CURVE
    fit <- gamlss(as.formula(paste0(outcome, " ~ cs(", "time", ",df=", df_m,")")), 
                  sigma.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_s,")")), 
                  nu.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_n,")")), 
                  tau.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_t,")")), 
                  data = plotdf %>% filter(.data$train_test == 1), family=dist)

    iqrplm <- centiles.pred(fit, type="centiles", xname = "time", xvalues=xvalues,
                            data = plotdf %>% filter(.data$train_test == 1),
                            plot = FALSE, cent = c(10,25,50,75,90))

    #-- PLOT for patient idnum

    plot_l <- ggplot(plotdf %>% filter(.data$train_test == 1)) +
        geom_point(aes(x = eval(parse(text = "time")), y = eval(parse(text = outcome)))) +
        geom_point(data =plotdf %>% 
               filter(.data$train_test ==  2) , 
               aes(x = eval(parse(text = "time")), y = eval(parse(text = outcome)), colour="Individual")) + 
        geom_point(data = plotdf %>%
                               filter(.data$train_test == 1), 
               aes(x = eval(parse(text = "time")), y = eval(parse(text =outcome)), colour="Matches")) + 
        geom_line(data = iqrplm, aes(x = eval(parse(text = "time")), y = eval(parse(text="C50")))) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = "time")), y = eval(parse(text="C25")), colour="50% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = "time")), y = eval(parse(text="C75")), colour="50% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = "time")), y = eval(parse(text="C10")), colour="80% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = "time")), y = eval(parse(text="C90")), colour="80% IQR")) +
        scale_colour_manual(name="", values=c("Individual"="orange", "50% IQR"="green",
                                              "80% IQR"="red", "Matches"="purple"
                                              ),
                            guide = guide_legend(fill = NULL, colour=NULL))+
        theme_bw() + theme(aspect.ratio = 1) +
      xlab(xlab) + ylab(ylab)

    return(
      list(
        plot = plot_l,
        plotdf = plotdf,
        predres = iqrplm 
        )
    )
}
