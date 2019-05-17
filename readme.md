
Predictive Mean Matched Sequential K-Nearest Neighbor
=====================================================

Introduction
------------

The purpose of the this repository is to provide a method for determining the trajectory of Knee Surgery Outcomes for patients based on obtaining predictions for patients using a 'patient-like-me' algorithm (a.k.a. [sequential k-nearest neighbor](https://www.ncbi.nlm.nih.gov/pubmed/20676226) (SKNN)). We extend the SKNN approach by matching similar patients using the [predictive mean matching](https://amstat.tandfonline.com/doi/abs/10.1080/07350015.1988.10509663).

Data
----

The data is coming from patient information recorded in the UC Health system (REDCap).

Code for Analysis
-----------------

[code folder](https://github.com/ck2136/PMMSKNN/tree/master/code) contains the analysis code and there a scripts in the [script folder](https://github.com/ck2136/PMMSKNN/tree/master/scripts). Testing on [TUG data set](https://github.com/ck2136/PMMSKNN/blob/master/data/tug_full.rda) is done with [`test_tug.R`](https://github.com/ck2136/PMMSKNN/tree/master/tests/testthat). The main prediction method is using the R package [brokenstick](https://github.com/stefvanbuuren/brokenstick), along with [predictive mean matching](https://books.google.com/books?hl=en&lr=&id=rM8eSRUYYHYC&oi=fnd&pg=PA442&dq=%22predictive+mean+matching%22++rubin&ots=OM-74mXZoX&sig=H-tIcTl7xqIfbgumXuHBktBTfkQ#v=onepage&q=%22predictive%20mean%20matching%22%20%20rubin&f=false) and [gamlss](https://www.gamlss.com/). Currently the code is under development to work within the [caret](https://github.com/topepo/caret) and [mlr](https://github.com/mlr-org/mlr) packages.

Installation/Compilation Tip
----------------------------

-   Download the github folder through

        devtools::install_github('ck2136/PMMSKNN')

-   If not available then `git clone` then `R CMD Install`

        git clone https://github.com/ck2136/PMMSKNN.git
        R CMD Install PMMSKNN

-   There will be dependencies that should be resolved if installation isn't done through the standard R method (in R):

        devtools::install_github("stefvanbuuren/brokenstick")
        devtools::install_deps('.')
        devtools::install_local('.')

Example workflow
----------------

### Load Libraries and TUG data

``` r
library("pacman")
p_load(PMMSKNN, readxl, dplyr, here, magrittr)
data(tug_full) ## example tug data
```

### Wrangle TUG data

``` r
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

# Need to have distinct patient id's for the full data
full %<>%
    distinct(patient_id, time, .keep_all=TRUE)
```

### preproc() creates matched test/train based on PMM

``` r
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
## time is not an integer! converting to integer! May need to check if this makes sense!
## Warning: number of observations (=1325) <= number of random effects
## (=1990) for term (0 + x1 + x2 + x3 + x4 + x5 | subjid); the random-effects
## parameters and the residual variance (or scale parameter) are probably
## unidentifiable
## boundary (singular) fit: see ?isSingular
test_proc %>% str(max.level=1)
## List of 5
##  $ train_post:Classes 'tbl_df', 'tbl' and 'data.frame':  1325 obs. of  9 variables:
##  $ train_o   :'data.frame':  398 obs. of  8 variables:
##  $ reg_obj   :List of 13
##   ..- attr(*, "class")= chr "lm"
##  $ test_post :Classes 'tbl_df', 'tbl' and 'data.frame':  602 obs. of  9 variables:
##  $ test_o    :Classes 'tbl_df', 'tbl' and 'data.frame':  201 obs. of  8 variables:
```

-   Depending on the `knots_exp` specified and the `out_time` time chosen, there will likely be warnings about parameter estimation within the `brokenstick()` algorithm. Here is where the researcher needs to consider the appropriate values for the variable in terms of clinical relevance and the data at hand.

### LOOCV: loocv\_function() calculates performance measure

``` r
res <- loocv_function(
  
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

  # Specify number of cores for parallel processing
  parallel=3,
  
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
## PREDICTING FOR TEST = 408
## PREDICTING FOR TEST = 443
## PREDICTING FOR TEST = 604
## PREDICTING FOR TEST = 450
## PREDICTING FOR TEST = 492
## PREDICTING FOR TEST = 536
## Current count is: 1
## PREDICTING FOR TEST = 516
## PREDICTING FOR TEST = 595
## PREDICTING FOR TEST = 564
## Current count is: 11
## Current count is: 21
## Current count is: 31
## PREDICTING FOR TEST = 565
## Current count is: 41
## Current count is: 51
## PREDICTING FOR TEST = 411
## Current count is: 61
## PREDICTING FOR TEST = 574
## Current count is: 71
## PREDICTING FOR TEST = 493
## PREDICTING FOR TEST = 541
## Current count is: 81
## Warning in predict.gamlss(obj, what = "sigma", newdata = newx, type = "response", : There is a discrepancy  between the original and the re-fit 
##  used to achieve 'safe' predictions 
## 
## Current count is: 91
## PREDICTING FOR TEST = 531
## Current count is: 101
## Current count is: 111
## PREDICTING FOR TEST = 520
## Current count is: 121
## Current count is: 131
## Current count is: 141
## Current count is: 151
## Current count is: 161
## Current count is: 171
## Current count is: 181
## Current count is: 191
## Current count is: 201
## PREDICTING FOR TEST = 434
## Current count is: 211
## Current count is: 221
## Current count is: 231
## Current count is: 241
## PREDICTING FOR TEST = 423
## Current count is: 251
## Current count is: 261
## Current count is: 271
## Current count is: 281
## Current count is: 291
## PREDICTING FOR TEST = 478
## Current count is: 301
## Current count is: 311
## PREDICTING FOR TEST = 476
## Current count is: 321
## Current count is: 331
## PREDICTING FOR TEST = 563
## Current count is: 341
## PREDICTING FOR TEST = 468
## PREDICTING FOR TEST = 433
## Current count is: 351
## PREDICTING FOR TEST = 526
## Current count is: 361
## Current count is: 371
## Current count is: 381
## PREDICTING FOR TEST = 473
## PREDICTING FOR TEST = 420
## PREDICTING FOR TEST = 461
## Current count is: 391
```

### Plots: plot\_cal() returns a plot of the performance measures from the LOOCV

``` r
plot_cal(plotobj = res, test_proc = test_proc, 
         obs_dist = "median")
```

<img src="README_figs/README-plotloocv-1.png" width="672" />

### Plots: plot\_cal() also returns plot of the calibration

``` r
plot_cal(plotobj = res, test_proc = test_proc,
         obs_dist = "median", loocv = FALSE)  
## [1] "creating training calibration plot"
## [1] "creating testing calibration plot"
## [1] "creating temp and test matching data"
## [1] "binding with decile"
## [1] "filt df made in test"
## [1] "Predicting TUG values"
```

<img src="README_figs/README-plotcal-1.png" width="672" />

``` r
      ## Note: specify loocv = FALSE 
```

Authors
-------

-   [Chong Hoon Kim](mailto:chong.kim@ucdenver.edu)
-   [Dr. Kathryn Colborn](mailto:KATHRYN.COLBORN@UCDENVER.EDU)
-   [Dr. Timothy Loar](mailto:TIMOTHY.LOAR@UCDENVER.EDU)
-   [Dr. Andrew Kittelson](mailto:andrew.kittelson@ucdenver.edu)
-   [Dr. Stef van Buuren](mailto:S.vanBuuren@uu.nl)
