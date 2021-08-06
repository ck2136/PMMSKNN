# FILEINFO {{{
# Filename      : test_preproc.R
# Purpose       : testing the output of preproc()
# Date created  : Tue 14 May 2019 07:08:51 PM MDT
# Last modified : Fri 17 May 2019 07:55:05 AM MDT
# Created by    : ck1
# Modified by   : ck1
# }}}
context("Preprocessing Data")

test_that(
  "preproc() works using ChickWeight Data",
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
    
    test_proc %>%
      expect_type("list") %>%
      expect_length(8) 
      
    expect_that(test_proc, is.list) 
    expect_is(test_proc$reg_obj, "lm")
    expect_true(any(grepl("test_post", names(test_proc))))
    expect_true(any(grepl("test_o", names(test_proc))))
    expect_true(any(grepl("train_o", names(test_proc))))
    expect_true(any(grepl("train_post", names(test_proc))))
    expect_true(any(grepl("reg_df", names(test_proc))))
    expect_true(any(grepl("reg_obj", names(test_proc))))
    expect_true(any(grepl("bs_obj", names(test_proc))))
    expect_true(any(grepl("varname", names(test_proc))))
  }
)

# }}}
