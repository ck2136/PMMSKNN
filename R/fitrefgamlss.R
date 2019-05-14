#' Fit reference GAMLSS on the full data
#' 
#' The function fits an initial GAMLSS distribution to the full data
#' in order to be used as a reference distribution for the following
#' patient like me GAMLSS fits
#' @param dist_fam  gamlss distribution specification
#' @param train_post - datasets, typically the \code{train_post} list component 
#'  of the object produced by \code{\link{preproc}}.
#' @param test_post Idem, component \code{test_post}
#' @param outcome    Name of the outcomes variable
#' @param time_elapsed - Name of the time variable. (type=string)
#' @param time_window - vector of numbers for `centiles.pred()`, `xvalues` argument 
#' @param cs Logical that specifies whether to use cubic spline. 
#'  The default \code{cs = FALSE} uses ...
#' @param dfspec Logical that specifies whether to the user sets  
#' degrees of freedom (...not clear to me what it does, and how it 
#' interacts with the next set of arguments)
#' @param d_f_m ... explain: arguments probably mean different things for different distributions. Might be preferable to package it with \code{dist_fam}
#' @param ptr_m -
#' @param d_f_s -
#' @param d_f_n - 
#' @param d_f_t -
#' @param \dots Passed down to \code{gamlss}
#' @return There are many possible return values 
#' @export
# - - - - - - - - - - - - - - - - - - - -#
# LOOCV Function ----
# - - - - - - - - - - - - - - - - - - - -#
# now using that we need to come up with the right number of matches that will give us the best bias and coverage.
fitrefgamlss <- function(
                         dist_fam = NULL, # for gamlss distribution
                         train_post, 
                         test_post, 
                         outcome, time_elapsed, 
                         time_window=NULL,
                         cs=FALSE,
                         dfspec=NULL,
                         d_f_m=1.8, ptr_m=1,
                         d_f_s=1.2,
                         d_f_n=1,
                         d_f_t=1,
                         ...) {

    # - - - - - - - - - - - - - - - - - - - - - # 
    # Setting spline option
    # - - - - - - - - - - - - - - - - - - - - - # 
    if(cs){
        if(is.null(dfspec)) {
            spl=paste0("cs(")
            spls=paste0("cs(")
            spln=paste0("cs(")
            splt=paste0("cs(")
        } else {
            spl=paste0("cs(df=",d_f_m,",")
            spls=paste0("cs(df=",d_f_s,",")
            spln=paste0("cs(df=",d_f_n,",")
            splt=paste0("cs(df=",d_f_t,",")
        }
    } else {
        spl="pb("
    }

    # - - - - - - - - - - - - - - - - - - - - - # 
    # Setting distribution for GAMLSS
    # If NULL then use NO as default
    # - - - - - - - - - - - - - - - - - - - - - # 
    if(is.null(dist_fam)) {

        # - - - - - - - - - - - - - - - - - - - - - # 
        # Reference gamlss fitting 
        # - - - - - - - - - - - - - - - - - - - - - # 
        message("Distribution of GAMLSS not chosen")
        message("Will determine based on Normal, Gamma, and Box Cox Cole and Green")

        ref1<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m,")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    data=na.omit(train_post), family=NO,
                    #gamlss control
                    ...)
        ref2<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m,")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    data=na.omit(train_post), family=GA,
                    #gamlss control
                    ...)
        ref3<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = as.formula(paste0(" ~ ",spln, time_elapsed, ")")),
                    data=na.omit(train_post), family=BCCGo,
                    #gamlss control
                    ...)
        ref4<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = as.formula(paste0(" ~ ",spln, time_elapsed, ")")),
                    data=na.omit(train_post), family=BCTo,
                    #gamlss control
                    ...)
        ref5<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = as.formula(paste0(" ~ ",spln, time_elapsed, ")")),
                    data=na.omit(train_post), family=BCPEo,
                    #gamlss control
                    ...)
        # here the BCCGo distribution with tau is intentionally left out bc the fit is good but not as when we jsut do the nu
        #ref5<-gamlss(eval(as.formula(paste0(outcome," ~ ", "pb(", time_elapsed, ")"))),
                    #sigma.formula = as.formula(paste0(" ~ ", "pb(", time_elapsed, ")")),
                    #nu.formula = as.formula(paste0(" ~ ", "pb(", time_elapsed, ")")),
                    #tau.formula = as.formula(paste0(" ~ ", "pb(", time_elapsed, ")")),
                    #data=na.omit(train_post), family=BCCGo)

        ref <- list(ref1,ref2,ref3,ref4,ref5)[which.min(lapply(list(ref1,ref2,ref3,ref4,ref5), function(x){x$aic}))][[1]]
    } else if(length(dist_fam()$parameters) == 2) {
        ref<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    data=data.frame(na.omit(train_post)), family=dist_fam,
                    #gamlss control
                    ...)
    } else if(length(dist_fam()$parameters) == 3) {
        ref<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = ~1,
                    data=na.omit(train_post), family=dist_fam,
                    #gamlss control
                    ...)
    } else if(length(dist_fam()$parameters) == 4) {
        ref<-gamlss(eval(as.formula(paste0(outcome," ~ ",spl, time_elapsed,"^",ptr_m, ")"))),
                    sigma.formula = as.formula(paste0(" ~ ",spls, time_elapsed, ")")),
                    nu.formula = ~1,
                    tau.formula = ~1,
                    data=na.omit(train_post), family=dist_fam,
                    ...
                    )
    }
    gamlss_dist <- ref$family[2]


    # - - - - - - - - - - - - - - - - - - - - - # 
    # Full training post operative data fitting (i.e reference data) with GAMLSS 
    # - - - - - - - - - - - - - - - - - - - - - # 

    # - - - - - - - - - - - - - - - - - - - - - # 
    # Specify time_window if NULL for GAMLSS centiles prediction
    # - - - - - - - - - - - - - - - - - - - - - # 
    if(is.null(time_window)){
        mint <- min(min(train_post[, time_elapsed]), min(test_post[, time_elapsed]))
        maxt <- max(max(train_post[, time_elapsed]), max(test_post[, time_elapsed]))
        iqrfull <- centiles.pred(ref, type='centiles', xname = time_elapsed, xvalues=c(mint:maxt),
                                data = na.omit(train_post),
                                cent=c(25,75), plot=FALSE)
    } else {
        iqrfull <- centiles.pred(ref, type="centiles", xname = time_elapsed, xvalues=c(time_window[1]:time_window[2]),
                                 data = na.omit(train_post),
                                 cent=c(25,75), plot=FALSE)
    }
    iqrfull$iqr<-iqrfull$C75-iqrfull$C25

    # - - - - - - - - - - - - - - - - - - - - - # 
    # Return the reference (ref) fit and iqr
    # - - - - - - - - - - - - - - - - - - - - - # 
    return(
           list(
                ref = ref,
                iqrfull = iqrfull,
                mint=mint,
                maxt=maxt,
                gamlss_dist=gamlss_dist,
                spl=spl, spls=spls, spln=spln, splt=splt
           )
    )

}
