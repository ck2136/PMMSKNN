# File:             test_tug_stef.R
# Objective:        One pass through all functions for TUG outcome

library(pcr)
library(dplyr) # ez data manipulation

# load only the TUG dataset
full  <- pcr::tug_full

# need to exclude the above patients
# exclude also time > 200
full <- full %>%
  #filter(!patient_id %in% exclude$patient_id & time < 200)
  filter(time < 200) %>%
  mutate(gender = as.factor(gender))

# Select patient id's that have TUG < 2 or > 70 after time > 3 
# 
# Need to exclude patients that have no post operative time beyond 
# 2 from the train pre and possibly test pre because if people 
# don't have post operative time in test it doesn't make sense
exclude <- full %>% 
  filter(tug < 2 | (tug > 70 & time > 3)) %>% 
  dplyr::select(patient_id) %>%
  bind_rows(full %>%
              group_by(patient_id) %>%
              filter(max(time) < 3) %>%
              distinct(patient_id))
full <- full %>%
  filter(!patient_id %in% exclude$patient_id & time < 200)

# Train and Test split for all TKA outcomes: create 
set.seed(1234)
full <- baselinemk(full, "patient_id", "time")

# make sure there aren't any missing data and etc
summary(full)
sapply(full, function(x) table(is.na(x)))

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
  dff = full,
  split_var = 'train_test',
  trainval = 1,
  testval = 2, 
  knots_exp = c(0, 14, 50, 90),
  out_time = 90, 
  outcome = "tug", 
  time_var = "time",
  pat_id = "patient_id", 
  varlist = c("age", "gender", "bmi"),
  filter_exp = NULL)

# TRAINING DATA: LOOCV Result
fin <- loocv_function(
  nearest_n = 10,
  train_post = test_proc$train_post,
  ord_data = test_proc$train_o,
  test_post = test_proc$test_post,
  test_o = test_proc$test_o,
  outcome = "tug",
  time_elapsed = "time",
  interval = 40,
  cs = TRUE,
  dfspec = TRUE,
  d_f_m = 3, ptr_m = 0.5, d_f_s = 1,
  dist_fam = gamlss.dist::BCCGo)

# Model performance plot with bias, coverage, width 50% pred in
plot_cal(plotobj = fin, 
         test_proc = test_proc, 
         obs_dist = "median")

# Calibration plot (does not yet work)
plot_cal(plotobj = fin, 
         test_proc = test_proc, 
         obs_dist = "median", 
         loocv = FALSE)

# Internal Validation Performance Measures
loocvperf(fin$loocv_res, test_proc)

# External Validation Performance Measures
extvalid(fin, test_proc)
