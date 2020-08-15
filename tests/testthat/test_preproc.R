# FILEINFO {{{
# Filename      : test_preproc.R
# Purpose       : testing the output of preproc()
# Date created  : Tue 14 May 2019 07:08:51 PM MDT
# Last modified : Fri 17 May 2019 07:55:05 AM MDT
# Created by    : ck1
# Modified by   : ck1
# }}}
context("Preprocessing Data")

# Load data and wrangel {{{
data("ChickWeight")

## Wrangle Steup {{{
# load only the TUG dataset
full  <- ChickWeight

# need to exclude the above patients
# exclude also time > 200
full <- full %>%
    #filter(!patient_id %in% exclude$patient_id & time < 200)
    filter(Time < 200) 

set.seed(1234)
full %<>%
  mutate(
    Chick = as.numeric(as.character(Chick)),
    train_test = ifelse(Chick %in% c(1,2,20, 30, 40), 2, 1),
    Diet = as.numeric(as.character(Diet))
  ) %>% 
  # Need to have distinct patient id's for the full data
  distinct(Chick, Time, .keep_all=TRUE)
    

## }}}

# }}}

# baselinemk() {{{
set.seed(1234)
full <- PMMSKNN:::baselinemk(full, "Chick", "Time")
test_that("baseline column created" , {
             expect_true(any(grepl("baseline", colnames(full))))
             expect_identical(max(full$baseline),1)
             expect_lte(min(full$baseline),0)
})
# }}}

# preproc() {{{
test_proc <- preproc(
                dff=full,                 # specify full dataset name
                split_var = 'train_test', # train test split variable
                trainval = 1,             # training set value
                testval = 2,              # test set value
                knots_exp = c(0, 4, 8, 16), # Specify broken stick knots
                out_time = 16,            # specify which timepoint to use 
                outcome = "weight",          # specify outcome variable name
                time_var = "Time",        # specify time variable name
                pat_id = "Chick",    # specify patient id variable name
                varlist = c("Diet") # specify list of covariates for pmm
                # filter_exp = "time > 3"   # Filter observations that will be included
)

test_that("preproc() contains appropriate data" , {
             expect_that(test_proc, is.list)
             expect_is(test_proc$reg_obj, "lm")
             expect_true(any(grepl("test_post", names(test_proc))))
             expect_true(any(grepl("test_o", names(test_proc))))
             expect_true(any(grepl("train_o", names(test_proc))))
             expect_true(any(grepl("train_post", names(test_proc))))
})
# }}}
