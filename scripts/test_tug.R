# FILEINFO {{{
# Filename      : C.K.
# Purpose       : Test out new functions in PMMSKNN package
# Date created  : somedate
# Last modified : Thu 16 May 2019 02:21:31 PM MDT
# Created by    : C.K.
# Modified by   : ck1
# }}}

# Load Libraries {{{
library("pacman")
p_load(PMMSKNN, readxl, dplyr, here)
# }}}

# Load Data and Wrangle {{{
data(tug_full)

## Wrangle Steup {{{
# load only the TUG dataset
full  <- tug_full

# need to exclude the above patients
# exclude also time > 200
full <- full %>%
    #filter(!patient_id %in% exclude$patient_id & time < 200)
    filter(time < 200) %>%
    mutate(gender = as.factor(gender))

# Select patient id's that have TUG < 2 or > 70 after time > 3 
exclude <- full %>% filter(tug < 2 | (tug > 70 & time > 3)) %>% dplyr::select(patient_id) %>%
    bind_rows(
              # need to exclude patients that have no post operative time beyond 2 from the train pre and possibly test pre because if people don't have post operative time in test it doesn't make sense
              full %>%
                  group_by(patient_id) %>%
                  filter(max(time) < 3) %>%
                  distinct(patient_id))

full <- full %>%
    filter(!patient_id %in% exclude$patient_id & time < 200)

# Train and Test split for all TKA outcomes: create 
set.seed(1234)
full <- PMMSKNN:::baselinemk(full, "patient_id", "time")

# make sure there aren't any missing data and etc
summary(full) ;  sapply(full, function(x) {
                          table(is.na(x))})

full <- full %>%
    distinct(patient_id, time, .keep_all=TRUE)
## }}}

# }}}

# PREPROCESS DATA {{{
## extract fitted values for the linear model (PREOP only here)
#-- preprocess train and test data involves:
###   1. Taking training post data and running broken stick to get y90 or yNU where NU is value chosen
###   2. Processes test data set so that yNU is matched with training data
###   3. The matched test and train data according to yNU is used to later match the personalized predicted gamlss

test_proc <- preproc(
                dff=full,                 # specify full dataset name
                split_var = 'train_test', # train test split variable
                trainval = 1,             # training set value
                testval = 2,              # test set value
                knots_exp = c(0, 14, 50, 90), # Specify broken stick knots
                out_time = 90,            # specify which timepoint to use 
                outcome = "tug",          # specify outcome variable name
                time_var = "time",        # specify time variable name
                pat_id = "patient_id",    # specify patient id variable name
                varlist = c("age","gender","bmi","b_tug"), # specify list of covariates for pmm
                filter_exp = "time > 3"   # Filter observations that will be included
)

# }}}

# RUN LOOCV on Training {{{

# LOOCV: Non Parallel {{{
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

# }}}

# LOOCV: Parallel {{{

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

# }}}

# }}}

# Plots {{{

# Bias, Cov, Pred {{{
plot_cal(plotobj = fin, test_proc = test_proc, obs_dist = "median")
# }}}

# Calibration {{{
plot_cal(plotobj = fin, test_proc = test_proc, obs_dist = "median", loocv = FALSE)
# }}}

# }}}

# Performance Measures {{{

# Internal Validation {{{
loocvperf(fin$loocv_res, 
          test_proc$train_o, 
          bias="raw",
          nearest_n=fin$loocv_score$nearest_n) 
# }}}

# External Validation {{{
extvalid(fin, test_proc)
# }}}

# Inspection of fin object {{{

fin %>% str(max.level = 2)
fin$pred_res$iqrcoverage
# there are many NA values in predicting the coverage for test set... NOT GOOD!

# }}}

# }}}
