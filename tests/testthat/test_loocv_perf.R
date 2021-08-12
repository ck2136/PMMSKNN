# FILEINFO {{{
# Filename      : test_loocv_perf.R
# Purpose       : testing the output of loocv_perf()
# Date created  : Tue 14 May 2019 07:08:51 PM MDT
# Last modified : Fri 17 May 2019 07:55:05 AM MDT
# Created by    : ck1
# Modified by   : ck1
# }}}
context("Preprocessing Data")

test_that(
  "loocv_perf() works to generate performance statistics",
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
    
    # Final default case
    test_proc <- preproc(
      dff=full,                 # specify full dataset name
      split_var = 'train_test', # train test split variable
      trainval = 1,             # training set value
      testval = 2,              # test set value
      knots_exp = c(0, 7, 14, 21), # Specify broken stick knots
      out_time = 21,            # specify which timepoint to use 
      outcome = "weight",          # specify outcome variable name
      time_var = "Time",        # specify time variable name
      id = "Chick",    # specify patient id variable name
      baseline_var = "baseline",
      m = 20,
      varlist = c("b_weight") # specify list of covariates for pmm
    )
    
    ## loocv_function() {{{
    fin <- loocv_function(
      
      # specify number or vector of numbers from {1,...,total number of patients in training data} 
      nearest_n = c(10,30,50),
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
    
    # create null pred_train
    fin_temp <- fin
    fin_temp$loocv_res$nearest_10$pred_train <- NULL
    
    # expect NA from the 1st row and 1st column
    expect_true(is.na(loocv_perf(
      fin_temp$loocv_res,
      outcome="weight",
      nearest_n=c(10,30,50),
      perf_round_by=4
    )[1,1]))
    
    # expect 
    
    fin_notune <- loocv_function(
      
      # specify number or vector of numbers from {1,...,total number of patients in training data} 
      nearest_n = c(10),
      userchoose = 10,
      loocv = FALSE,
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
    
    loocv_perf(
      fin_notune$pred_res,
      outcome="weight",
      nearest_n=c(10),
      perf_round_by=4,
      train = FALSE
    )
    
  }
)

# }}}
