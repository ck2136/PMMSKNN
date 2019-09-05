#' Plot Nth Percentile patient like me plots
#' 
#' Creates two plots.
#' \enumerate{
#'   \item \emph{Lower N_1-th plot} showing individual patient values 
#'   along with matched individuals and gamlss centile curves for the
#'   lower N-th vector value specified
#'   \item \emph{Upper N_2-th plot} showing individual patient values 
#'   along with matched individuals and gamlss centile curves for the
#'   upper N-th vector value specified
#'   }
#'   
#' @param plotobj   - An object produced by \code{\link{loocv_function}}
#' @param test_proc - Preprocessed object from \code{\link{preproc}}
#' @param outcome   - Name of the outcomes variable (type=string)
#' @param time_var  - Name of the time variable. (type=string)
#' @param nvec      - Vector of two proportions (e.g. \code{nvec = c(0.2,0.8)}).
#' These values indicate the lower and upper n-th percentile/100 person to plot.
#' @param mtype - Integer value indicating matching type. Default is set to 1 which follows the
#' matching of patients based on recommendation from \href{https://stefvanbuuren.name/fimd/sec-pmm.html}{van Buuren et al.}. \code{mtype} values are 
#' from \code{0} to \code{4}
#' @param n     - Integer representing number of matches.
#' @param dist  - Object of \code{\link{gamlss.dist}} specifying the distribution for gamlss
#' @param df_m - Numeric value that specifies the degrees of freedom for the cubic spline specified for the mean parameter of the distribution specified according to \code{dist_fam}
#' @param df_s - Numeric value that specifies the degrees of freedom for the cubic spline specified for the scale parameter of the distribution specified according to \code{dist_fam}
#' @param df_n - Numeric value that specifies the degrees of freedom for the cubic spline specified for the shape parameter, specifically the \eqn{\nu} parameter, of the distribution specified according to \code{dist_fam}
#' @param df_t - Numeric value that specifies the degrees of freedom for the cubic spline specified for the shape parameter, specifically the \eqn{\tau} parameter, of the distribution specified according to \code{dist_fam}
#' @param xvalues - Vector of values for the prediction centile to be plotted
#' @param \dots   - Options to specify in the plotting of the percentile curves.
#' 
#' @return An object of class \code{ggplot} that outputs a gamlss predition curve for 2 individuals in terms of their predicted values over time.
#' 
#' @export
plot_NthP_plm <- function(plotobj,
                     test_proc=test_proc,
                     outcome = "tug",
                     time_var = "time",
                     nvec=c(0.2,0.8),
                     mtype = 1,
                     n=10,
                     dist=BCCGo,
                     df_m=2,
                     df_s=1,
                     df_n=1,
                     df_t=1,
                     xvalues=3:200,
                     ...) {
  # Main input: object from LOOCV function which spits out train data output
  # plotobj$pred_res$pred contains the training set predictions in dataframe
  # plotobj$pred_res

    #-- DATASET for Lower N_1-th patient
    matches <- matchIdExtractTest(test_proc = test_proc, 
                                  mtype = mtype, 
                                  i = extractIdbyPerf(test_proc, nvec[1]),
                                  n = n) 

    matchmodel <- test_proc$train_post %>% 
        filter(.data$patient_id %in% matches)

    plotdf <- test_proc$test_post %>%
        filter(.data$patient_id == extractIdbyPerf(test_proc, nvec[1])) %>%
        bind_rows(
                  matchmodel
        )
    plotdf1 <- plotdf

    #-- Lower N_1-th  PORTION
    fit <- gamlss(as.formula(paste0("tug", " ~ cs(", "time", ",df=", df_m,")")), 
                  sigma.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_s,")")), 
                  nu.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_n,")")), 
                  tau.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_t,")")), 
                  data = plotdf %>% filter(.data$train_test == 1), family=dist)

    iqrplm <- centiles.pred(fit, type="centiles", xname = time_var, xvalues=xvalues,
                            data = plotdf %>% filter(.data$train_test == 1),
                            plot = FALSE, cent = c(10,25,50,75,90))

    #-- PLOT for lower N_1-th patient 

    plot_l <- ggplot(plotdf %>% filter(.data$train_test == 1)) +
        geom_point(aes(x = eval(parse(text = time_var)), y = eval(parse(text = outcome)))) +
        geom_point(data =plotdf %>% 
               filter(.data$train_test ==  2) , 
               aes(x = eval(parse(text = time_var)), y = eval(parse(text = outcome)), colour="Patient")) + 
        geom_point(data = plotdf %>%
                               filter(.data$train_test == 1), 
               aes(x = eval(parse(text = time_var)), y = eval(parse(text =outcome)), colour="Matches")) + 
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C50")))) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C25")), colour="50% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C75")), colour="50% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C10")), colour="80% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C90")), colour="80% IQR")) +
        scale_colour_manual(name="", values=c("Patient"="orange", "50% IQR"="green",
                                              "80% IQR"="red", "Matches"="purple"
                                              ),
                            guide = guide_legend(fill = NULL, colour=NULL))+
        theme_bw() + theme(aspect.ratio = 1)

        
    #-- DATASET for Upper N_2-th patient
    matches <- matchIdExtractTest(test_proc = test_proc, 
                                  mtype = mtype, 
                                  i = extractIdbyPerf(test_proc, nvec[2]),
                                  n = n) 

    matchmodel <- test_proc$train_post %>% 
        filter(.data$patient_id %in% matches)

    plotdf <- test_proc$test_post %>%
        filter(.data$patient_id == extractIdbyPerf(test_proc, nvec[2])) %>%
        bind_rows(
                  matchmodel
        )

    #-- Lower N_1-th  PORTION
    fit <- gamlss(as.formula(paste0("tug", " ~ cs(", "time", ",df=", df_m,")")), 
                  sigma.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_s,")")), 
                  nu.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_n,")")), 
                  tau.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_t,")")), 
                  data = plotdf %>% filter(.data$train_test == 1), family=dist)

    iqrplm <- centiles.pred(fit, type="centiles", xname = time_var, xvalues=xvalues,
                            data = plotdf %>% filter(.data$train_test == 1),
                            plot = FALSE, cent = c(10,25,50,75,90))

    #-- PLOT for lower N_1-th patient 

    plot_u <- ggplot(plotdf %>% filter(.data$train_test == 1)) +
        geom_point(aes(x = eval(parse(text = time_var)), y = eval(parse(text = outcome)))) +
        geom_point(data =plotdf %>% 
               filter(.data$train_test ==  2) , 
               aes(x = eval(parse(text = time_var)), y = eval(parse(text = outcome)), colour="Patient")) + 
        geom_point(data = plotdf %>%
                               filter(.data$train_test == 1), 
               aes(x = eval(parse(text = time_var)), y = eval(parse(text =outcome)), colour="Matches")) + 
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C50")))) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C25")), colour="50% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C75")), colour="50% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C10")), colour="80% IQR")) +
        geom_line(data = iqrplm, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C90")), colour="80% IQR")) +
        scale_colour_manual(name="", values=c("Patient"="orange", "50% IQR"="green",
                                              "80% IQR"="red", "Matches"="purple"
                                              ),
                            guide = guide_legend(fill = NULL, colour=NULL))+
        theme_bw() + theme(aspect.ratio = 1)




    # find minimum and maximum of x and y coordinates
    # usually the max of ref will be used for both x and y

    prow <- plot_grid(plot_l + theme(legend.position = "none") + 
                      xlim(0, max(max(plotdf[,time_var]),max(plotdf1[,time_var]))) +
                      ylim(0, max(max(plotdf[,outcome]),max(plotdf1[,outcome]))) +
                      labs(title = paste0("Patient Reference Chart \n ", nvec[1]*100,
                                          " Percentile Patient" , "(id = ", 
                                          extractIdbyPerf(test_proc, nvec[1]),")"), 
                           x = "Time (days)", y = toupper(outcome))
                      ,
                      plot_u + theme(legend.position = "none") +
                          xlim(0, max(max(plotdf[,time_var]),max(plotdf1[,time_var]))) +
                          ylim(0, max(max(plotdf[,outcome]),max(plotdf1[,outcome]))) +
                          labs(title = paste0("Patient Reference Chart \n ", nvec[2]*100,
                                              " Percentile Patient" , "(id = ", 
                                              extractIdbyPerf(test_proc, nvec[2]),")"), 
                               x = "Time (days)", y = toupper(outcome))
                          ,
                          align = "vh"
    )
    legend <- get_legend(plot_u)
    p <- plot_grid(prow, legend, rel_widths = c(3,.3))
    return(p)
}
