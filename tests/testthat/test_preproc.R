# FILEINFO {{{
# Filename      : test_preproc.R
# Purpose       : testing the output of preproc()
# Date created  : Tue 14 May 2019 07:08:51 PM MDT
# Last modified : Tue 14 May 2019 09:19:46 PM MDT
# Created by    : ck1
# Modified by   : ck1
# }}}
context("Preprocessing Data")

# Load Libs and data  {{{
library("pacman")
p_load(PMMSKNN)
data(tug_full)
# }}}

# baselinemk() {{{
set.seed(1234)
tug_full <- PMMSKNN:::baselinemk(tug_full, "patient_id", "time")
test_that("baseline column created" , {
             expect_true(any(grepl("baseline", colnames(tug_full))))
             expect_identical(max(tug_full$baseline),1)
             expect_lte(min(tug_full$baseline),0)
})
# }}}

# preproc() {{{
tug_full <- tug_full %>%
    distinct(patient_id, time, .keep_all=TRUE)
test_proc <- preproc(
                dff=tug_full,                 # specify full dataset name
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

test_proc%>% str(max.level=1)
   
test_that("preproc() contains appropriate data" , {
             expect_that(test_proc, is.list)
             expect_is(test_proc$reg_obj, "lm")
             expect_true(any(grepl("test_post", names(test_proc))))
             expect_true(any(grepl("test_o", names(test_proc))))
             expect_true(any(grepl("train_o", names(test_proc))))
             expect_true(any(grepl("train_post", names(test_proc))))
})
# }}}
