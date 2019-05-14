# - - - - - - - - - - - - #
# Plotting Function: Precision Curve ----
# - - - - - - - - - - - - #

plot_prec <- function(plotobj, n_range){
  precplot <- rbindlist(lapply(plotobj$loocv_res, function(x) { 
    rbindlist(x$precisionvec, idcol="id") }), 
    idcol="n") %>%
    filter(grepl(n_range, .data$n)) %>%
    ggplot(., aes(x=.data$time, y=.data$prec, group=.data$id)) +
    geom_line() +
    geom_point() +
    xlim(0, 365) +
    ylim(0, 10) +
    facet_wrap( ~.data$n)
  
  return(precplot)
  
}
