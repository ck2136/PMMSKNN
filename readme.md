
Predictive Mean Matched Sequential K-Nearest Neighbor
=====================================================

Introduction
------------

The purpose of the this repository is to provide a method for determining the trajectory of Knee Surgery Outcomes for patients based on obtaining predictions for patients using a 'patient-like-me' algorithm (a.k.a. [sequential k-nearest neighbor](https://www.ncbi.nlm.nih.gov/pubmed/20676226)).

Data
----

The data is coming from patient information recorded in the UC Health system (REDCap).

Code for Analysis
-----------------

/code folder contains the analysis code and there a scripts in the /script folder. Testing on TUG data set is done with `test_tug.R`. The main prediction method is using the R package [brokenstick](https://github.com/stefvanbuuren/brokenstick), along with [predictive mean matching](https://books.google.com/books?hl=en&lr=&id=rM8eSRUYYHYC&oi=fnd&pg=PA442&dq=%22predictive+mean+matching%22++rubin&ots=OM-74mXZoX&sig=H-tIcTl7xqIfbgumXuHBktBTfkQ#v=onepage&q=%22predictive%20mean%20matching%22%20%20rubin&f=false) and [gamlss](https://www.gamlss.com/). Currently the code is under development to work within the [caret](https://github.com/topepo/caret) and [mlr](https://github.com/mlr-org/mlr) packages.

Installation/Compilation Tip
----------------------------

-   Download the github folder through

        devtools::install_github('ck2136/PMMSKNN')

-   If not available then `git clone` then `R CMD Install`

        git clone https://github.com/ck2136/PMMSKNN.git
        R CMD Install PMMSKNN

-   There will be dependencies that should be resolved if installation isn't done through the standard R method (in R):

        ddvtools::install_github("stefvanbuuren/brokenstick")
        devtools::install_deps('.')
        devtools::install_local('.')

Example workflow
----------------

``` r
library("pacman")
p_load(PMMSKNN, readxl, dplyr, here, magrittr)
data(tug_full) ## example tug data
```

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

# Need 
full %<>%
    distinct(patient_id, time, .keep_all=TRUE)
```

Authors
-------

-   [Chong Hoon Kim](mailto:chong.kim@ucdenver.edu)
-   [Dr. Kathryn Colborn](mailto:KATHRYN.COLBORN@UCDENVER.EDU)
-   [Dr. Timothy Loar](mailto:TIMOTHY.LOAR@UCDENVER.EDU)
-   [Dr. Andrew Kittelson](mailto:andrew.kittelson@ucdenver.edu)
-   [Dr. Stef van Buuren](mailto:S.vanBuuren@uu.nl)
