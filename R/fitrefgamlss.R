#' Fit reference GAMLSS on the full data
#' 
#' The function fits an initial GAMLSS distribution to the full data
#' in order to be used as a reference distribution for the following
#' matched individuals' GAMLSS fits
#' 
#' @param dist_fam  Gamlss distribution specification using the \code{\link{gamlss.dist}} package. The specification for a normal distribution would be \code{gamlss.dist::NO}. For other distributions see \code{\link{gamlss.dist}}.
#' @param train_post Data frame that contains the post-baseline observations from the training dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param test_post Data frame that contains the post-baseline observations from the testing dataset. Typically this would be the \code{train_post} list component that was generated from the \code{\link{preproc}} function
#' @param outcome    Name of the outcomes variable (type=string)
#' @param time_elapsed Name of the time variable. (type=string)
#' @param time_window Vector of numbers for `centiles.pred()`, `xvalues` argument. For example, specify such as \code{c(10:30)}
#' @param cs Logical that specifies whether to use cubic spline. 
#'  The default \code{cs = FALSE} uses ...
#' @param dfspec Logical (\code{TRUE/FALSE}) that specifies whether to 
#' specify degrees of freedoms for the location, scale, and shape parameters
#' for the distribution specified with \code{dist_fam}.
#' Default value is \code{NULL}.
#' @param d_f_m  Numeric value that specifies the degrees of freedom for the cubic spline specified for the mean parameter of the distribution specified according to \code{dist_fam}
#' @param ptr_m  Numeric value that specifies the power transformation of time variable. Default value is 1.
#' @param d_f_s  Numeric value that specifies the degrees of freedom for the cubic spline specified for the scale parameter of the distribution specified according to \code{dist_fam}
#' @param d_f_n  Numeric value that specifies the degrees of freedom for the cubic spline specified for the shape parameter, specifically the \eqn{\nu} parameter, of the distribution specified according to \code{dist_fam}
#' @param d_f_t  Numeric value that specifies the degrees of freedom for the cubic spline specified for the shape parameter, specifically the \eqn{\tau} parameter, of the distribution specified according to \code{dist_fam}
#' @param \dots Passed down to \code{gamlss}
#' 
#' @return Returns a gamlss object as described in \code{\link{gamlss}}
#' 
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
        mint <- min(time_window)
        maxt <- max(time_window)
        iqrfull <- centiles.pred(ref, type="centiles", xname = time_elapsed, xvalues=time_window,
                                 data = na.omit(train_post),
                                 cent=c(25,75), plot=FALSE)
    }
    # iqrfull$iqr<-iqrfull$C75-iqrfull$C25
    iqrfull$iqr<-iqrfull[[3]]-iqrfull[[2]]

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
