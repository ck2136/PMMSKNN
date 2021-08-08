#' PreProcess Function: Split full data into train/test and match individuals using predictive mean matching
#' 
#' \code{preproc()} function essentially takes a full dataset along with a list of arguments
#' necessary to split the data into train/test datasets. After splitting the data 
#' into a training and testing split, The training data is used to fit a linear mixed model
#' with a b-s;line using the \code{brokenstick} package. The user specifies the knots 
#' and the distal outcome of interest (\eqn{\hat{y}}) which essentially is used to perform the 
#' predictive mean matching. The \eqn{\hat{y}} is used to fit a linear model based on the 
#' user specified list of matching characteristics provided through the `varlist` argument.
#' 
#' @param dff  Full dataset containing both training and testing dataset specified by 
#' `split_var` column.
#' @param split_var  A string representing the name of the splitting variable used to split
#' the full dataset (i.e. `dff`) into training and testing data. THe variable/column should be
#' numeric type.
#' @param trainval  A numeric value indicating training set observations. 
#' If \code{train_val = 1}, then rows in the \code{dff} should have 1 in each of the
#' `split_var` column.
#' @param testval  A numeric value indicating testing set observations. 
#' If \code{test_val = 0}, then rows in the \code{dff} should have 0 in each of the
#' `split_var` column.
#' Expression can be of form \code{"time > 3"}. Default is \code{NULL}
#' @param knots_exp   Numeric vector with break points for \code{brokenstick} model
#' @param out_time    Single numeric value specifying the distal time. Has to be one of 
#'                   \code{knots_exp}.
#' @param outcome     String representing the name of the outcomes variable
#' @param time_var    String representing the name of the time variable
#' @param idcol      String representing the name of the variable with patient identification
#' @param baseline_var  String representing the name of the pre-op/post-op 
#' indicator variable. For example, for pre-op (baseline = 1); 
#'                   for post-op: (baseline = 0), or otherwise
#' @param varlist    Names of additional variables for prediction. If not 
#'                   specified, all variables in \code{dff} are predictors. 
#'                   Categorical variables may be factor or character.
#' @param pmmform    - formula representing the model used for predictive mean matching. For example, to regress \code{outcome} on the variables in \code{varlist}, \code{outcome ~ var1 + var2 + var3} 
#' @param modelselect  - A logical (\code{TRUE/FALSE}) specifying whether to go through a stepwise selection of variables for the predictive mean matching algorithm
#' @param \dots        - Specification for linear model in the predictive mean matching algorithm.
#' @param m           - Numeric value representing the Number of repititions of obtaining \eqn{\dot{y}} (i.e. the predicted value from predictive mean matching)
#' 
#' @return  A list with six components. 
#'          1. Post-baseline training data 
#'          2. Dataframe with training set patient id and \eqn{\dot{y}} values ordered
#'          3. Regression dataframe used for the predictive mean matching. the `yhat` column here is the predicted mean values.
#'          4. Predictive mean matching model object
#'          5. Post-baseline testing data 
#'          6. Dataframe with testing set patient id and \eqn{\dot{y}} values ordered
#' 
#' @export
preproc <- function(dff,
                    split_var = 'train_test',
                    trainval = 1,
                    testval = 2,
                   knots_exp = c(0, 14, 50, 90),    # select knots that are clinically relevat
                   out_time = 90,                   # this variable has to be within the knots above
                   outcome = "tug",
                   time_var = "time",
                   idcol = "patient_id",           # Need the id column to be specified as "patient_id"
                   baseline_var = "baseline",       # pre-op/post-op indicator variable. preop: baseline = 1; postop: baseline = 0 or otherwise
                   varlist = NULL,                  # need user to fill in var list otherwise will use all vars, categorical variables need to be factor or character
                   pmmform = NULL,
                   modelselect = FALSE,
                   m=5,
                   ...) {

    # - - - - - - - - - - - - - - - - - - - - - - #
    # If baseline_var not  supplied, then stop
    # - - - - - - - - - - - - - - - - - - - - - - #
    if(is.null(dff[,baseline_var])){
        stop("baseline_var is NULL. 
             Specify  baseline var as string. (e.g.baseline_var = 'baseline').
             Utility function baselinemk() may be used to create baseline variable.
             ")
    }
    # - - - - - - - - - - - - - - - - - - - - - - #
    # If time variable isn't integer then matching won't be stable due to floating point error
    # - - - - - - - - - - - - - - - - - - - - - - #
    if(!is.integer(dff[,time_var])){
        message(paste0(time_var, " is not an integer! converting to integer! May need to check if this makes sense!"))
        dff[,time_var] <- as.integer(dff[[time_var]])
    }

    # - - - - - - - - - - - - - - - - - - - - - - #
    # If varlist supplied, then form dataframe that only contains 
    # ID, Outcome, Time, train/test split variable, and Listed variables 
    # else stop function
    # - - - - - - - - - - - - - - - - - - - - - - #
    if(!is.null(varlist)){
        dff <- dff %>%
            .[,c(idcol, outcome, time_var, split_var, baseline_var, varlist)] %>%
          as.data.frame
    } else {
        stop("varlist not populated: specify varlist = c('var1','var2',...)")
    }
    # - - - - - - - - - - - - - - - - - - - - - - #
    # Non complete covariate matrix! Stop
    # - - - - - - - - - - - - - - - - - - - - - - #
    if(
       any(
           !dff[,c(idcol, outcome, time_var, split_var, baseline_var, varlist)] %>% complete.cases
       )
       ){
        stop("missing data in the data!")
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - #
    # Check for duplicated values in baseline and postoperative
    # - - - - - - - - - - - - - - - - - - - - - - #
    if(any(dff %>% filter(.$baseline == 1) %>% dplyr::select(!!dplyr::sym(idcol)) %>% unlist %>% duplicated)){
        warning("Duplicate baseline values exist within training and testing set! remove them before running preproc")
    }
  
  if(
    any(c(
      # Training set duplicated?
      dff %>%
      dplyr::filter(.$train_test == 1 & .$baseline == 1) %>%
      dplyr::select(!!dplyr::sym(idcol)) %>%
      unlist %>% duplicated,
      # Test set duplicated?
      dff %>%
      dplyr::filter(.$train_test == 2 & .$baseline == 1) %>%
      dplyr::select(!!dplyr::sym(idcol)) %>%
      unlist %>% duplicated
    ))
  ){
    warning("Duplicate baseline values exist within either the training and testing set! remove them before running preproc")
    }

    # if(
    #    dff %>% 
    #        filter(.data$baseline == 0) %>% 
    #        dplyr::select_(idcol, time_var) %>% nrow != 
    #        dff %>% 
    #        filter(.data$baseline == 0) %>% 
    #        dplyr::select_(idcol, time_var) %>% 
    #        distinct_(idcol, time_var)  %>% nrow
    # ){
    #     stop("Duplicate post operative values exist! remove them before running preproc")
    # }

    # - - - - - - - - - - - - - - - - - - - - - - #
    # Check whether train and test cases have both pre and post op values
    # - - - - - - - - - - - - - - - - - - - - - - #
    exclude <- dff %>% 
        filter(.data$baseline == 1) %>%
        dplyr::select(!!dplyr::sym(idcol), !!dplyr::sym(outcome)) %>% 
        full_join(
                  dff %>% 
                      filter(.data$baseline ==1) %>%
                      distinct(!!dplyr::sym(idcol), .keep_all =TRUE) %>%
                      dplyr::select(!!dplyr::sym(idcol), !!dplyr::sym(outcome)) %>% 
                      dplyr::rename("p_outcome" = !!dplyr::sym(outcome)),
                  by = c(idcol)
                  ) %>%
        filter(is.na(.data$p_outcome) | is.na(.[outcome])) 

    if(nrow(exclude) != 0){
        message(paste0("patients ", exclude," don't have post-operative values" ))
        stop("Not all patients have Post-Op values for LOOCV!")
    }

    # - - - - - - - - - - - - - - - - - - - - - - #
    # Split test/train 
    # - - - - - - - - - - - - - - - - - - - - - - #

    # df_train <-dff  %>% filter_(paste0(split_var, "==", trainval)) 
    # df_test <- dff %>% filter_(paste0(split_var, "==",testval))
    df_train <-dff  %>% filter(!!dplyr::sym(split_var) == trainval)
    df_test <- dff %>% filter(!!dplyr::sym(split_var) == testval)

    # - - - - - - - - - - - - - - - - - - - - - - #
    # Split test/train by pre and post using baseline_var
    # - - - - - - - - - - - - - - - - - - - - - - #

    # pre_train_df <- df_train %>% filter_(paste0(baseline_var,"== 1"))
    # pre_test_df <- df_test %>% filter_(paste0(baseline_var, "== 1"))
    pre_train_df <- df_train %>% filter(!!dplyr::sym(baseline_var) == 1)
    pre_test_df <- df_test %>% filter(!!dplyr::sym(baseline_var) == 1)

    # - - - - - - - - - - - - - - - - - - - - - - #
    # Filter out the post-baseline observations
    # - - - - - - - - - - - - - - - - - - - - - - #
    
    # post_train_df <- df_train %>% filter_(paste0(baseline_var, "== 0"))
    # post_test_df <- df_test %>% filter_(paste0(baseline_var, "== 0"))
    post_train_df <- df_train %>% filter(!!dplyr::sym(baseline_var) == 0)
    post_test_df <- df_test %>% filter(!!dplyr::sym(baseline_var) == 0)

    # - - - - - - - - - - - - - - - - - - - - - - #
    # use brokenstick to predict values at knots_exp
    # - - - - - - - - - - - - - - - - - - - - - - #

    # fit <- brokenstick(y = unlist(post_train_df[,outcome]),
    #                    x = unlist(post_train_df[,time_var]),
    #                    subjid = unlist(post_train_df[,idcol]),
    #                    knots = knots_exp
    #                    )
    
    fit <- brokenstick(
      formula(paste0(outcome, " ~ ", time_var , " | ", idcol)),
      data = post_train_df,
      knots = knots_exp
    )

    # est1 <- predict(fit, post_train_df, at="knots")
    # est1 <- predict(fit, post_train_df, x="knots", shape="wide")
    est1 <- predict(fit, post_train_df, x="knots", shape="long") %>% dplyr::select(!!dplyr::sym(idcol), !!dplyr::sym(time_var), .pred) 
    
    alldf <- left_join(
                       est1[round(est1[[time_var]], 3) == round(out_time,3),] %>%
                       dplyr::select(!!dplyr::sym(idcol), .pred) %>%
                       setNames(c(idcol,"yhat")) %>%
                          # need to change id bc it is a factor when outputted through predict()
                          mutate(!!idcol := !!parse_quo(paste0("as.numeric(as.character(",idcol,"))"), env=rlang::caller_env()))
                          # mutate(!!idcol := !!parse_quosure(paste0("as.numeric(as.character(",idcol,"))")))
                      ,
                     pre_train_df %>%
                       .[,c(idcol, time_var, split_var, baseline_var, varlist)],
                      by = idcol
                      ) %>%
    # - - - - - - - - - - - - - - - - - - - - - - #
    # only keep complete cases 
    # - - - - - - - - - - - - - - - - - - - - - - #
    .[complete.cases(.),]
    # - - - - - - - - - - - - - - - - - - - - - - #
    # lm() using gamlss package to get predicted outcome at out_time
    # - - - - - - - - - - - - - - - - - - - - - - #
    #pmm <- gamlss(yhat ~ get(outcome) + age + gender + bmi,family=NO, data=alldf)
    #pmm <- lm(yhat ~ get(outcome) + age + gender + bmi, data=alldf)
    if(is.null(pmmform)){
        pmm <- lm(formula(paste0("yhat ~ ", paste0(varlist, collapse="+"))), data=alldf)
    } else {
        pmm <- lm(formula = pmmform, data=alldf )
    }

    # - - - - - - - - - - - - - - - - - - - - - - #
    # If modelselect = TRUE, use stepAIC to choose variables
    # - - - - - - - - - - - - - - - - - - - - - - #
    if(modelselect){
        pmm <- stepAIC(pmm)
    }

    # - - - - - - - - - - - - - - - - - - - - - - #
    # Create dataset with fitted outcome at out_time for training patients
    # - - - - - - - - - - - - - - - - - - - - - - #
    train_ordered <- alldf %>%
        dplyr::select(!!dplyr::sym(idcol)) %>%
        cbind(pmm$fitted.values) %>%
        dplyr::rename(id = !!dplyr::sym(idcol),
               Fitted=`pmm$fitted.values`) %>%
               #Fitted=`pmm$mu.fv`) %>%
        arrange(.data$Fitted)

    train_ordered <- train_ordered %>%
        bind_cols(
                  ydot = pmmydotgen(alldf, 
                                    formula = formula(paste0("yhat ~ ", paste0(varlist, collapse="+"))),
                                    m, 
                                    pat_id = idcol, 
                                    seed=1234, 
                                    dftest=NULL)
        )

    # - - - - - - - - - - - - - - - - - - - - - - #
    # Create dataset with fitted outcome at out_time for testing patients
    # Here we still use the linear model used to fit the training data (i.e. pmm)
    # - - - - - - - - - - - - - - - - - - - - - - #
    
    test_ordered <- pre_test_df %>% 
        dplyr::select(!!dplyr::sym(idcol)) %>%
        bind_cols(pred = predict(pmm, data=alldf, 
                                 newdata=pre_test_df %>% 
                                            .[,c(outcome, varlist)]
                             )
        ) %>%
        rename(id = !!dplyr::sym(idcol)) %>%
        arrange(.data$pred)

    test_ordered <- test_ordered %>%
        bind_cols(
                  ydot = pmmydotgen(alldf, 
                                    formula = formula(paste0("yhat ~ ", paste0(varlist, collapse="+"))),
                                    m, 
                                    pat_id = idcol, 
                                    seed=1234, 
                                    dftest=pre_test_df)
        )

    # - - - - - - - - - - - - - - - - - - - - - - #
    # Change patient_id column for LOOCV function use
    # - - - - - - - - - - - - - - - - - - - - - - #
    post_train_df <- post_train_df 

    post_test_df <- post_test_df 

    return(list(train_post = post_train_df, 
                train_o =  train_ordered,
                reg_df = alldf,
                reg_obj = pmm,
                test_post = post_test_df, 
                test_o =  test_ordered,
                bs_obj = fit,
                varname = c(outcome, time_var, idcol, baseline_var, split_var)
                )
    )
}
