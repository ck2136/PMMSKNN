# FILEINFO {{{
# Filename      : test_baselinemk.R
# Purpose       : testing the output of baselinemk()
# Date created  : Thu 5 Aug 2021 07:08:51 PM MDT
# Last modified : Thu 5 Aug 2021 07:55:05 AM MDT
# Created by    : ck1
# Modified by   : ck1
# }}}
context("Creating baseline variable Data")

test_that(
  "preproc() works using ChickWeight Data",
  {
    data("ChickWeight") 
    
    ## Create new data to input into loocv_function*()
    CWn <- ChickWeight %>%
      baselinemk("Chick","Time") 
    
    CWn %>%
      expect_type("list") %>%
      expect_length(5) 
      
    expect_true(any(grepl("baseline", names(CWn))))
  }
)

# }}}