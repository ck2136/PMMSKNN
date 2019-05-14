#' Plot reference and patient like me plot
#' 
#' Creates two plots.
#' \enumerate{
#'   \item \emph{Reference plot} showing reference population values 
#'   along with fitted gamlss centile curves 
#'   \item \emph{Patient like me plot} showing the distribution of the 
#'   outcome of matched patients along with the selected patients' predicted
#'   centile curve. 
#'   }
#' @param plotobj An object produced by \code{\link{loocv_function}}
#' @param test_proc - An object produced by \code{\link{preproc}}
#' @param outcome Name of the outcome variable
#' @param time_var Name of the time variable
#' @param patnum Integer value of the patient number to plot alongside
#' the population reference chart
#' @param dist object of gamlss.dist specifying the distribution for gamlss
#' @param df_m  Number representing the degrees of freedome for the location parameter
#' @param df_s  Number representing the degrees of freedome for the sigma parameter
#' @param df_n  Number representing the degrees of freedome for the nu parameter
#' @param df_t  Number representing the degrees of freedome for the tau parameter
#' @param xvalues  Vector of values for the prediction centile to be plotted
#' @param \dots Not used
#' @return A list with components ..., or an object of class 
#' \code{ggplot} [??]
#' @export
plot_ref_plm <- function(plotobj,
                     test_proc=test_proc,
                     outcome = "tug",
                     time_var = "time",
                     patnum=219,
                     dist=NO,
                     df_m=2,
                     df_s=1,
                     df_n=1,
                     df_t=1,
                     xvalues=3:200,
                     ...) {
  # Main input: object from LOOCV function which spits out train data output
  # plotobj$pred_res$pred contains the training set predictions in dataframe
  # plotobj$pred_res

    #-- REFERENCE PORTION
    ref <- gamlss(as.formula(paste0("tug", " ~ cs(", "time", ",df=", df_m,")")), 
                  sigma.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_s,")")), 
                  nu.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_n,")")), 
                  tau.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_t,")")), 
                  data = test_proc$train_post, family=dist)

    iqrref <- centiles.pred(ref, type="centiles", xname = time_var, xvalues=xvalues,
                            data = test_proc$train_post, 
                            plot = FALSE, cent = c(10,25,50,75,90))

    #-- TRAINING MATCHES
    if(!any(test_proc$train_post$patient_id %in% patnum)){
        stop("patnum is not a id number in the training set! Select another patient number in the training set")
    }

    refplot <- ggplot(test_proc$train_post) +
        geom_point(aes(x = eval(parse(text = time_var)), y = eval(parse(text = outcome)))) +
        geom_point(data =test_proc$train_post %>% 
                   filter(test_proc$train_post$patient_id %in% patnum), 
               aes(x = eval(parse(text = time_var)), y = eval(parse(text = outcome)), colour= "Patient")) + 
        geom_line(data = iqrref, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C50")))) +
        geom_line(data = iqrref, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C25")), colour="50% IQR")) +
        geom_line(data = iqrref, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C75")), colour="50% IQR")) +
        geom_line(data = iqrref, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C10")), colour="80% IQR")) +
        geom_line(data = iqrref, aes(x = eval(parse(text = time_var)), y = eval(parse(text="C90")), colour="80% IQR")) +
        scale_colour_manual(name="", values=c("Patient"="orange", "50% IQR"="green",
                                              "80% IQR"="red", "Matches"="purple"
                                              ),
                            guide = guide_legend(fill = NULL, colour=NULL))+
        theme_bw() + theme(aspect.ratio = 1)
        


    #-- PATIENTLIKEME APPROACH
    # choose a patient at random

    # find closest match to patient in training to training
    matches <- test_proc$train_o %>% 
        bind_cols(diff = abs(test_proc$train_o[test_proc$train_o$id %in% patnum,"Fitted"] - test_proc$train_o$Fitted)) %>%
        arrange(diff) %>%
        dplyr::select(.data$id) %>%
        .[-1,] %>%
        head(n = plotobj$nearest_n) %>% unlist %>% as.vector

    iqrplm <- test_proc$train_post %>%
        filter(test_proc$train_post$patient_id %in% matches) %>%
        gamlss(as.formula(paste0(outcome, " ~ cs(", "time", ",df=", df_m,")")), 
                      sigma.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_s,")")), 
                      nu.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_n,")")), 
                      tau.formula = as.formula(paste0(" ~ cs(", "time", ",df=", df_t,")")), 
                      data = ., family=dist) %>%
        centiles.pred(., type="centiles", xname = time_var, xvalues=xvalues,
                      data = test_proc$train_post %>%
                          filter(test_proc$train_post$patient_id  %in% matches), plot = FALSE, cent = c(10,25,50,75,90)
        )

    plmplot <- ggplot(test_proc$train_post %>%
                      filter(test_proc$train_post$patient_id %in% matches)
                  ) +
        geom_point(aes(x = eval(parse(text = time_var)), y = eval(parse(text = outcome)))) +
        geom_point(data =test_proc$train_post %>% 
                   filter(test_proc$train_post$patient_id %in% patnum) , 
               aes(x = eval(parse(text = time_var)), y = eval(parse(text = outcome)), colour="Patient")) + 
        geom_point(data = test_proc$train_post %>%
                   filter(test_proc$train_post$patient_id %in% matches), 
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


    prow <- plot_grid(refplot + theme(legend.position = "none") + 
                      xlim(0, max(test_proc$train_post[,time_var])) +
                      ylim(0, max(test_proc$train_post[,outcome])) +
                      labs(title = paste0("Population Reference Chart \n ", "of N = ", nrow(test_proc$train_o)), 
                           x = "Time (days)", y = toupper(outcome))
                      ,
                      plmplot + theme(legend.position = "none") +
                          xlim(0, max(test_proc$train_post[,time_var])) +
                          ylim(0, max(test_proc$train_post[,outcome])) +
                          labs(title = paste0("Patients Like Me Chart \n ", "of m = ", length(matches)),
                               x = "Time (days)", y = toupper(outcome))
                          ,
                          align = "vh"
    )
    legend <- get_legend(plmplot)
    p <- plot_grid(prow, legend, rel_widths = c(3,.3))
    return(p)
}
