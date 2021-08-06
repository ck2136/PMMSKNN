# FILEINFO {{{
# Filename      : test_loocv.R
# Purpose       : Test LOOCV Function
# Date created  : Fri 17 May 2019 07:42:49 AM MDT
# Last modified : Fri 17 May 2019 07:53:53 AM MDT
# Created by    : ck1
# Modified by   : ck1
# }}}

## Load packages
library('dplyr')
library('brokenstick')

# LOOCV: Non Parallel {{{
context("Running LOOCV using ChickWeight Data")

test_that(
  "LOOCV works using ChickWeight Data",
  {
    data("ChickWeight") 
    
    set.seed(1234)
    
    ## Select 50% as training set
    cid <- ChickWeight %>% 
      as.data.frame %>%
      distinct(Chick) %>%
      pull() %>%
      sample(.,size = {.[] %>% length}*0.5)
    
    ## Create new data to input into loocv_function*()
    full <- ChickWeight %>%
      mutate(
        train_test = ifelse(Chick %in% cid, 1, 2),
      ) %>%
      ## Specify baseline weight
      left_join(
        ChickWeight %>%
          group_by(Chick) %>%
          filter(row_number() == 1) %>%
          dplyr::select(weight, Chick) %>%
          rename(b_weight = weight)
      ) %>%
      # Set Time == 0 observations as baseline = 1; 0 otherwise
      mutate(
        baseline = ifelse(Time == 0, 1, 0),
        Chick = as.numeric(as.character(Chick))
      ) 
    
    test_proc <- preproc(
      dff=full,                 # specify full dataset name
      split_var = 'train_test', # train test split variable
      trainval = 1,             # training set value
      testval = 2,              # test set value
      knots_exp = c(0, 7, 14, 21), # Specify broken stick knots
      out_time = 21,            # specify which timepoint to use 
      outcome = "weight",          # specify outcome variable name
      time_var = "Time",        # specify time variable name
      pat_id = "Chick",    # specify patient id variable name
      baseline_var = "baseline",
      m = 20,
      varlist = c("b_weight") # specify list of covariates for pmm
    )
    
    ## loocv_function() {{{
    fin <- loocv_function(
      
      # specify number or vector of numbers from {1,...,total number of patients in training data} 
      nearest_n = c(13:15),
      # enter training and testing post operative and fitted y90 dataset
      preproc=test_proc,
      # interval = 10,
      mtype=0,
      # Specify use of cubic spline or not
      cs=TRUE,
      # specify degrees of freedom use or not
      dfspec=TRUE,
      # specify degree of freedom for location, scale and shape (d_f_* where * = {m, s} for location and scale default for shape is 1.
      # specify power transformation of location (ptr_m)
      d_f_m=3, ptr_m=0.5,
      d_f_s=1,
      # Specify distribution for location, scale and shape 
      #dist_fam = gamlss.dist::NO)
      dist_fam = gamlss.dist::NO,
      perfrank = "totscore"
    )
    
    expect_that(names(fin), equals(c("pred_res","test_score","loocv_res","loocv_score","nearest_n")))
    expect_false(is.null(fin$nearest_n))
  }
)

## }}}

# }}}

