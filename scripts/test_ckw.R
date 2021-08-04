# FILEINFO {{{
# Filename      : C.K.
# Purpose       : Test out chickweight
# Date created  : somedate
# Date created  : Mon 25 Feb 2020 06:52:50 PM MST
# Last modified : Mon 25 Feb 2020 06:52:50 PM MST
# Created by    : C.K.
# Modified by   : ck
# }}}

# Load Libraries {{{ ----------------
rm(list=ls())
library("pacman")
p_load(
  PMMSKNN,readxl,dplyr,here,testthat,gamlss, 
  doParallel,future.apply,janitor,fastDummies
  )
devtools::reload(pkgload::inst("PMMSKNN"))
# }}}

# Load Data and Wrangle {{{ --------------
data("ChickWeight")

## Select 50% as training set
set.seed(1234)
cid <- ChickWeight %>% 
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
    ) %>%
  dummy_cols("Diet")

# }}}

# PREPROCESS DATA {{{ ---------------
## extract fitted values for the linear model (PREOP only here)
#-- preprocess train and test data involves:
###   1. Taking training post data and running broken stick to get y90 or yNU where NU is value chosen
###   2. Processes test data set so that yNU is matched with training data
###   3. The matched test and train data according to yNU is used to later match the personalized predicted gamlss
#debug(preproc)
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
                varlist = c("Diet_2","Diet_3","Diet_4","b_weight") # specify list of covariates for pmm
)

# }}}
summary(test_proc$reg_obj)

# RUN LOOCV on Training {{{ ---------------------

# RUN LOOCV on Training with BrokenStick {{{ ----------------------------

res_bs <- loocv_function_bs(
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  nearest_n = seq(10, 20, 1),
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "weight",
  mtype=0, # Use straight up predicted value
  # m = 10, # if we're using 
  # interval = 10,
  # loocv = FALSE,
  # userchoose = 9,
  bs_obj = test_proc$bs_obj,
  parallel = 7L
)

res_bs$loocv_score
res_bs$test_score
res_bs$nearest_n

res_bs2 <- loocv_function_bs(
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  nearest_n = seq(10, 30, 1),
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "weight",
  mtype=0, # Use straight up predicted value
  # m = 10, # if we're using 
  # interval = 10,
  # loocv = FALSE,
  # userchoose = 9,
  bs_obj = test_proc$bs_obj,
  parallel = 7L,
  perfrank = "totscore",
  opt_cov = 0.5
)

res_bs2$test_score
res_bs2$loocv_score
res_bs2$nearest_n
## }}}

# LOOCV PMM: Non Parallel {{{ ------------------
res_ps1 <- loocv_function(
  
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  nearest_n = c(10:20),
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "weight",
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

# }}}
res_ps1$loocv_score
res_ps1$test_score

test_that("LOOCV performance list created" , {
             expect_that(names(res_ps1), equals(c("pred_res","test_score","loocv_res","loocv_score","nearest_n")))
             expect_false(is.null(res_ps1$nearest_n))
})


# LOOCV PMM: Parallel {{{ -----------------
res_ps2 <- loocv_function(
  
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  nearest_n = c(10:25),
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "weight",
  #outcome = "knee_flex",
  interval = NULL,
  mtype=0,
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
  parallel=7, m=5,
  # Specify distribution for location, scale and shape 
  #dist_fam = gamlss.dist::NO)
  dist_fam = gamlss.dist::BCCGo)

# }}}
res_ps2$loocv_score
res_ps2$nearest_n
res_ps2$test_score

# LOOCV_SKNN {{{ -------------------------
sknnlooObj<- loocv_function_sknn(
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  fulldata = full,
  nearest_n = c(10:14),
  patid = "Chick",
  formula =  ~ Diet_2 + Diet_3 + Diet_4 + b_weight,
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "weight",
  # interval = 10,
  
  # Specify use of cubic spline or not
  cs=TRUE,
  
  # specify degrees of freedom use or not
  dfspec=TRUE,
  
  # specify degree of freedom for location, scale and shape (d_f_* where * = {m, s} for location and scale default for shape is 1.
  # specify power transformation of location (ptr_m)
  d_f_m=3, ptr_m=0.5,
  d_f_s=1,
  # parallel
  parallel=5, m=5,
  
  # Specify distribution for location, scale and shape 
  #dist_fam = gamlss.dist::NO)
  dist_fam = gamlss.dist::NO
)

sknnlooObj$loocv_score
sknnlooObj$test_score
sknnlooObj$nearest_n

# Try ranking by totscore
sknnlooObj2<- loocv_function_sknn(
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  fulldata = full,
  nearest_n = c(15,30,45),
  patid = "Chick",
  formula =  ~ Diet_2 + Diet_3 + Diet_4 + b_weight,
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "weight",
  # interval = 10,
  
  # Specify use of cubic spline or not
  cs=TRUE,
  
  # specify degrees of freedom use or not
  dfspec=TRUE,
  
  # specify degree of freedom for location, scale and shape (d_f_* where * = {m, s} for location and scale default for shape is 1.
  # specify power transformation of location (ptr_m)
  d_f_m=3, ptr_m=0.5,
  d_f_s=1,
  # parallel
  parallel=5, m=5,
  
  # Specify distribution for location, scale and shape 
  #dist_fam = gamlss.dist::NO)
  dist_fam = gamlss.dist::NO,
  perfrank = "totscore"
)

sknnlooObj2$loocv_score
sknnlooObj2$test_score
sknnlooObj2$nearest_n

## }}}

## Plot Nth Person {{{
plot_NthP_plm(
              test_proc=test_proc,
              outcome = "weight",
              nvec=c(0.1,0.9),
              mtype = 1,
              n=10,
              dist=BCCGo,
              df_m=2,
              df_s=1,
              df_n=1,
              df_t=1,
              xvalues=0:21,
              chartname="Reference Growth Chart",
              name="Chick"
              )

## }}}

## Plot individual trajectory {{{

plotobj <- plot_ind(
  test_proc=test_proc,
  outcome = "weight",
  idnum = 1,
  mtype = 1,
  n=10,
  dist=NO,
  df_m=2,
  df_s=1,
  df_n=1,
  df_t=1,
  xvalues=0:21,
  xlab="Time", ylab="Chick Weight (grams)"
)
# Plot 
plotobj$plot
# Plot Data Frame
plotobj$plotdf
## }}}


## Extract patients based on id {{{
## in the below function, the parameter 'i' is the idnumber 
matches <- matchIdExtractTest(test_proc = test_proc, 
                              mtype = 1, 
                              i = 1,
                              n = 10) 

matchmodel <- test_proc$train_post %>% 
  filter(.data$patient_id %in% matches)

## }}}

## test_that() {{{------------------------

test_that("LOOCV performance list created" , {
             expect_that(names(res_bs), equals(c("pred_res","test_score","loocv_res","loocv_score","nearest_n")))
             expect_false(is.null(res_bs$nearest_n))
})

## }}}

# }}}

# Plots {{{ ----------------

## Bias, Cov, Pred {{{
plot_cal(plotobj = res_ps1, 
         test_proc = test_proc, 
         obs_dist = "median",
         outcome = "weight",
         filt=FALSE,
         pred_sum="mean",
         #plot_by=seq(10,150,5),
         loocv=TRUE,
         filter_exp = NULL,
         plot_cal_zscore=FALSE,
         wtotplot=FALSE,
         plotvals=FALSE,
         iqrfull=NULL,
         bs=FALSE
         )
## }}}

## Calibration {{{
#debug(plot_cal)
plot_cal(plotobj = res_ps1, 
         test_proc = test_proc, 
         obs_dist = "median",
         outcome = "weight",
         filt=FALSE,
         pred_sum="mean",
         #plot_by=seq(10,150,5),
         loocv=FALSE,
         filter_exp = NULL,
         plot_cal_zscore=FALSE,
         wtotplot=FALSE,
         plotvals=FALSE,
         iqrfull=NULL,
         bs=FALSE
         )
## }}}

# }}}

# Profiling {{{

## LOOCV BS {{{

profvis({
   res_bs <- loocv_function_bs(
                               # specify number or vector of numbers from {1,...,total number of patients in training data} 
                               nearest_n = seq(10, 30, 10),
                               # enter training and testing post operative and fitted y90 dataset
                               train_post = test_proc$train_post,
                               ord_data = test_proc$train_o,
                               test_post = test_proc$test_post,
                               test_o = test_proc$test_o,
                               # Specify outcome variable and time variable name
                               outcome = "tug",
                               time_elapsed = "time",
                               mtype=0, # Use straight up predicted value
                               # m = 10, # if we're using 
                               # interval = 10,
                               # loocv = FALSE,
                               # userchoose = 9,
                               bs_obj = test_proc$bs_obj,
                               parallel = 7L
   )
})
## }}}

# LOOCV PMM: Non Parallel {{{ ------------------

profvis({
   res_ps1 <- loocv_function(

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
                             dist_fam = gamlss.dist::NO)
})

# }}}

# LOOCV PMM: Parallel {{{ -----------------

res_ps2 <- loocv_function(
  
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  nearest_n = c(35:37),
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
  mtype=0,
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
  parallel=3, m=5,
  # Specify distribution for location, scale and shape 
  #dist_fam = gamlss.dist::NO)
  dist_fam = gamlss.dist::BCCGo)

# }}}
res_ps2$loocv_score
res_ps2$test_score

# LOOCV_SKNN {{{ -------------------------

sknnlooObj<- loocv_function_sknn(
  # specify number or vector of numbers from {1,...,total number of patients in training data} 
  fulldata = full,
  nearest_n = c(15,30,45),
  formula =  ~ age + bmi + b_tug,
  # enter training and testing post operative and fitted y90 dataset
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  # Specify outcome variable and time variable name
  outcome = "tug",
  time_elapsed = "time",
  # interval = 10,
  
  # Specify use of cubic spline or not
  cs=TRUE,
  
  # specify degrees of freedom use or not
  dfspec=TRUE,
  
  # specify degree of freedom for location, scale and shape (d_f_* where * = {m, s} for location and scale default for shape is 1.
  # specify power transformation of location (ptr_m)
  d_f_m=3, ptr_m=0.5,
  d_f_s=1,
  # parallel
  parallel=5, m=5,
  
  # Specify distribution for location, scale and shape 
  #dist_fam = gamlss.dist::NO)
  dist_fam = gamlss.dist::NO
)

sknnlooObj$test_score

## }}}

# }}}
