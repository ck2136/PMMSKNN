#' Generate a DataFrame for producing a calibration plot (for Continuous Value Predictions)
#' 
#' This is a utility function that takes two vectors of observed and predicted values 
#' in order to generate a dataframe necessary for producing a calibration plot.
#' 
#' @param observed      Vector of observed values (Numeric type).
#' @param predicted     Vector of predicted values (Numeric type).
#' @param n             Number of splits for observed and predicted.
#' Default is 10.
#' 
#' @return A data frame used to generate calibration plots (using \code{"ggplot"}) for training and testing data.
#' 
#' @export
calplotDFgen <- function(observed,
                         predicted,
                         n = 10
                         ) {

    # AGGREGATED DATAFRAME
    df <- bind_cols(observed = observed,
                             predicted = predicted) %>%
                    bind_cols(
                              dec = predicted %>% ntile(., 10)
                              ) %>%
                    left_join(
                              bind_cols(observed = observed,
                                        predicted = predicted) %>%
                              bind_cols(
                                        dec = predicted %>% ntile(., 10)
                              ) %>%
                              group_by(.data$dec) %>%
                              summarise(avg_val =  mean(.data$predicted)),
                          by = "dec"
                    )

    if(any(is.na(df$dec))){
      print("NA's present in group")
      df <- df %>%
        filter(!is.na(.data$dec))
    }

    # PREDICTED VALUE BY DECILE DATAFRAME
    dfpreddec <- bind_cols(observed = observed, 
                             predicted = predicted) %>%
                    bind_cols(
                              dec = predicted %>% ntile(., 10)
                              ) %>%
                    group_by(.data$dec) %>%
                    summarise(avg_val =  mean(.data$predicted))
        
    # Observed Data Portion
    dffitmed <- df %>%
        group_by(.data$dec) %>%
        do(fitmed = MedianCI(.[,"observed"][[1]], conf.level= 0.95,
                             method = "exact", R = 10000))
        #do(fitpois = glm(tug ~ 1, data= .))

    dfmedcoef <- tidy(dffitmed, fitmed)

    # GENERATE PLOT DATAFRAME
    plotdf <- dfpreddec %>%
        left_join(spread(dfmedcoef, .data$names, -.data$dec), by = "dec") %>%
        rename(ul = .data$upr.ci,
               ll = .data$lwr.ci)

    # set minimum and maximum based on the range of predicted and observed TUG values
    minc <- floor(min(plotdf$ll, plotdf$avg_val, na.rm=TRUE)) 
    maxc <- ceiling(max(plotdf$ul, plotdf$avg_val, na.rm=TRUE)) 
    return(
           list(
                #df = df,
                #dfpreddc=dfpreddc,
                plotdf = plotdf, 
                minc = minc, maxc = maxc
           )
    )
}
