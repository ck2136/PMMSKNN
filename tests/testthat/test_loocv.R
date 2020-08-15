# FILEINFO {{{
# Filename      : test_loocv.R
# Purpose       : Test LOOCV Function
# Date created  : Fri 17 May 2019 07:42:49 AM MDT
# Last modified : Fri 17 May 2019 07:53:53 AM MDT
# Created by    : ck1
# Modified by   : ck1
# }}}

# LOOCV: Non Parallel {{{
context("LOOCV Using ChickWeight Data")

data("ChickWeight")
full  <- ChickWeight
full %<>%
  mutate(
    Chick = as.numeric(as.character(Chick)),
    train_test = ifelse(Chick %in% c(1,2,20, 30, 40), 2, 1),
    Diet = as.numeric(as.character(Diet))
  ) %>% 
  rename(time = Time) %>%
  # Need to have distinct patient id's for the full data
  distinct(Chick, time, .keep_all=TRUE)
set.seed(1234)
full <- PMMSKNN:::baselinemk(full, "Chick", "time")
test_proc <- preproc(
                dff=full,                 # specify full dataset name
                split_var = 'train_test', # train test split variable
                trainval = 1,             # training set value
                testval = 2,              # test set value
                knots_exp = c(0, 4, 8, 16), # Specify broken stick knots
                out_time = 16,            # specify which timepoint to use 
                outcome = "weight",          # specify outcome variable name
                time_var = "time",        # specify time variable name
                pat_id = "Chick",    # specify patient id variable name
                varlist = c("Diet") # specify list of covariates for pmm
                # filter_exp = "time > 3"   # Filter observations that will be included
)

## loocv_function() {{{
fin <- loocv_function(
  
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  nearest_n = c(13:14),
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "weight",
  time_elapsed = "time",
  interval = 10,
  
  # Specify use of cubic spline or not
  cs=TRUE,
  
  # specify degrees of freedom use or not
  dfspec=TRUE,
  
  # specify degree of freedom for location, scale and shape (d_f_* where * = {m, s} for location and scale default for shape is 1.
  # specify power transformation of location (ptr_m)
  d_f_m=3, ptr_m=0.5,
  d_f_s=1, m=5,
  
  # Specify distribution for location, scale and shape 
  #dist_fam = gamlss.dist::NO)
  dist_fam = gamlss.dist::NO)

## }}}

## test_that() {{{

test_that("LOOCV performance list created" , {
             expect_that(names(fin), equals(c("pred_res","loocv_res","loocv_score","nearest_n")))
             expect_false(is.null(fin$nearest_n))
})

## }}}

# }}}

# LOOCV: Parallel {{{

## loocv_function() {{{
fin <- loocv_function(
  
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  nearest_n = c(15:17),
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "weight",
  #outcome = "knee_flex",
  time_elapsed = "time",
  interval = NULL,
  
  # Specify use of cubic spline or not
  cs=TRUE,
  
  # specify degrees of freedom use or not
  #dfspec=NULL,
  dfspec=TRUE,
  
  # specify degree of freedom for location, scale and shape (d_f_* where * = {m, s} for location and scale default for shape is 1.
  # specify power transformation of location (ptr_m)
  d_f_m=3, ptr_m=0.5,
  #d_f_m=3, ptr_m=1,
  d_f_s=1,
  # parallel
  parallel=3,m=5,
  
  # Specify distribution for location, scale and shape 
  #dist_fam = gamlss.dist::NO)
  dist_fam = gamlss.dist::NO)

## }}}

## test_that() {{{

test_that("LOOCV performance list created" , {
             expect_that(names(fin), equals(c("pred_res","loocv_res","loocv_score","nearest_n")))
             expect_false(is.null(fin$nearest_n))
})

## }}}

# }}}

