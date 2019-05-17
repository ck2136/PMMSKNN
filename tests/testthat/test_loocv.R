# FILEINFO {{{
# Filename      : test_loocv.R
# Purpose       : Test LOOCV Function
# Date created  : Fri 17 May 2019 07:42:49 AM MDT
# Last modified : Fri 17 May 2019 07:53:53 AM MDT
# Created by    : ck1
# Modified by   : ck1
# }}}

# LOOCV: Non Parallel {{{

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
  outcome = "tug",
  time_elapsed = "time",
  interval = 10,
  
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
  outcome = "tug",
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
  parallel=3,
  
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

