#' Fit reference GAMLSS on the full data
#' 
#' The function fits an initial GAMLSS distribution to the full data
#' in order to be used as a reference distribution for the following
#' patient like me GAMLSS fits
#' @param ref  Reference model from the \code{\link{fitrefgamlss}}.
#' @param matchmodel - dataset generated from the \code{\link{matchTrainDataGen}} 
#' @return There are many possible return values 
#' @export
# - - - - - - - - - - - - - - - - - - - - - - #
# Fit GAMLSS model to nearest n matches
# - - - - - - - - - - - - - - - - - - - - - - #
plm <- function(
                ref=ref,
                matchmodel=matchmodel
                ) {

    mint <- ref$mint; maxt <- ref$maxt; gamlss_dist <- ref$gamlss_dist; 
    spl <- ref$spl; spls <- ref$spls; spln <- ref$spln; splt <- ref$splt; ref <- ref$ref;

    out <- tryCatch(
                    {
                        #message("TRY PART")
                        update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE))

                    },
                    error = function(cond)
                    {
                        message(paste("ERROR in GAMLSS"))
                        #message("ORIGINAL ERROR:")
                        message(cond)
                        return(NA)
                    },
                    warning=function(cond)
                    {
                        message(paste("GAMLSS WARNING"))
                        #message("ORIIGNAL WARNING:")
                        message(cond)
                        return(update(ref, data=matchmodel, control = gamlss.control(n.cyc=1000, trace=FALSE)))
                        #return(NULL)
                    },
                    finally={
                        #message("PROCESSING GAMLSS")
                        #message("GAMLSS EOL SUCCESS")
                    }
    )
    return(out)
}
