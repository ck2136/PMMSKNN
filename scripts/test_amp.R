# - - - - - - - - - - - - - - - - - - - - - #
# File Name: 	'test_amp.R'
# Purpose:
# Created on: 	22-02-2019
# Modified on 	17-12-2018
# Created by:
# Contact Info: ck
# - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - #
# Environment configuration
# - - - - - - - - - - - - - - - - - - - - #

## Remove objects in workspace
rm(list=ls())
## Set working directory
library(pacman)
p_load(here, readxl, dplyr, rio,pcr)
# if you are in the code folder get out!
#library("pcr")
setwd(here())

# - - - - - - - - - - - - - - - - - - - - #
# Data Import and testing
# - - - - - - - - - - - - - - - - - - - - #

amp  <- pcr::amp

#-- add baseline indicator
amp <- baselinemk(amp, "patient_id", "time")

#-- add baseline value of plusm as predictor
amp <- amp %>% 
    left_join(
              amp %>% 
                  arrange(patient_id, time) %>% 
                  distinct(patient_id, .keep_all = TRUE) %>%
                  select(patient_id, plusm) %>%
                  rename(b_plusm = plusm)
    ) 

#-- filter out duplicated values

amp <- amp %>%
    distinct(patient_id, time, .keep_all=TRUE) 

#-- Need to have patients that have both the pre op AND post op values
amp %>% 
    filter(baseline == 1) %>%
    dplyr::select(patient_id, plusm) %>% 
    left_join(
              amp %>% 
                  filter(baseline ==0) %>%
                  distinct(patient_id, .keep_all =TRUE) %>%
                  dplyr::select(patient_id, plusm) %>%
                  rename(p_plusm = plusm)
    ) %>%
    filter(is.na(p_plusm)) %>% nrow

exclude <- amp %>% 
    filter(baseline == 1) %>%
    dplyr::select(patient_id, plusm) %>% 
    left_join(
              amp %>% 
                  filter(baseline ==0) %>%
                  distinct(patient_id, .keep_all =TRUE) %>%
                  dplyr::select(patient_id, plusm) %>%
                  rename(p_plusm = plusm)
    ) %>%
    filter(is.na(p_plusm)) %>%
    dplyr::select(patient_id) %>% unlist %>% as.vector

amp <- amp %>%
    filter(!patient_id %in% exclude)

# quick data check
# make sure there aren't any missing data and etc
summary(amp) ;  sapply(amp, function(x) {
                          table(is.na(x))})



#- - - - - - - - - - - - - - - - - - - - - - - - - - #
## Run the Long LOOCV Training Set Matchem Up!
#- - - - - - - - - - - - - - - - - - - - - - - - - - #

##extract fitted values for the linear model (PREOP only here)
#-- preprocess train and test data: involves
#1. Taking training post data and running broken stick to get y90 or yNU where NU is value chosen
#2. Processes test data set so that yNU is matched with training data
#3. The matched test and train data according to yNU is used to later match the personalized predicted gamlss


# Process the full data to get
# 1. Post operative test data 
# 2. Fitted 90-Day TUG for test patients
# 3. Post operative train data 
# 4. Fitted 90-Day TUG for train patients
test_proc <- preproc(
                dff=amp, 
                split_var = 'train_test', # train test split variable
                #split_var = 'source', # for all_tka dataset
                trainval = 1, # training set value
                testval = 2, # test set value
                knots_exp = c(0, 50, 181), # Specify broken stick knots
                out_time = 181,  # specify which timepoint to use 
                #outcome = "tug", # specify outcome variable name
                outcome = "plusm",
                time_var = "time", # specify time variable name
                pat_id = "patient_id", # specify patient id variable name
                varlist = c("age","gender","bmi", "b_plusm"), # specify list of covariates for pmm
                filter_exp = NULL               
)


# - - - - - - - - - - - - - - -#
# TRAINING DATA: LOOCV Result
# - - - - - - - - - - - - - - -#

# need to run everything in parallel by 3 and store results for later use
# let's say we do 5 at a time for the 3 cores then we go up by 15

# run from 5 to 100, 5 to 14 run with just a sequential forlopp 
fin <- loocv_function(nearest_n = seq(50,50,10),
                      train_post = test_proc$train_post,
                      ord_data = test_proc$train_o,
                      test_post = test_proc$test_post,
                      test_o = test_proc$test_o,
                      outcome = "plusm",time_elapsed = "time",
                      cs=TRUE,
                      dfspec=TRUE,
                      matchprobweight=FALSE,
                      d_f_m=3, ptr_m=0.5,
                      d_f_s=1,
                      userchoose=NULL,                    # User gets to choose number 
                      seed=1234,
                      parallel=3,
                      dist_fam = gamlss.dist::NO
                      )
# run things in parallel but be cognizant about memory

#saveRDS(fin,paste0("./data/fin_BCCGo_50_150_10_cs_dfspec_plusm.RDS"))# saved as the Gamma distribution
#rm(fin)


# - - - - - - - - - - - - - - - - - - - - #
# Plots: Calibration, Bias/RMSE and Coverage
# - - - - - - - - - - - - - - - - - - - - #

#-- Load Final Results

#fin <- readRDS(paste0(here(),"/data/fin_BCCGo_50_150_10_cs_dfspec_plusm.RDS"))

#-- Bias, Cov, Pred
#png("./report/figure/LOOCV_BCCGo_plusm_test.png")
plot_cal(
         # Specify plotobj which is the object saved from the loocv_func()
         # if multiple distributes were used then specify the distribution (e.g. myfiles$BCCGo) otherwise specify just the object (e.g. fin)
         #plotobj=finfilt,
         #plotobj=finnofilt,
         plotobj=fin,
         # specify the processed file to match test data and train data for calibration
          test_proc=test_proc,
         # specify calibration plot x 
          pred_sum='mean',
          # specify observed distribution to use (default = "median")
          obs_dist="median",
          wtotplot=FALSE
          )
#dev.off()

#-- Calibration
#png("./report/figure/CAL_BCCGo_plusm_test.png")
plot_cal(
         #plotobj=myfiles$BCCGo,
         plotobj=fin,
         #plotobj=finfilt,
         #plotobj=finnofilt,
          test_proc=test_proc,
          outcome = 'plusm',
          #outcome = 'tug',
          #outcome = 'knee_flex',
          pred_sum='mean',
          #obs_dist="gaussian",
          #obs_dist="poisson",
          obs_dist="median",
          # specify loocv=FALSE to check calibration
          loocv = FALSE
          )
#dev.off()

# - - - - - - - - - - - - - - - - - - - - #
# Internal Validation Performance Measures
# - - - - - - - - - - - - - - - - - - - - #

loocvperf(fin$loocv_res, test_proc$train_o)

# - - - - - - - - - - - - - - - - - - - - #
# External Validation Performance Measures
# - - - - - - - - - - - - - - - - - - - - #

extvalid(fin, test_proc)
