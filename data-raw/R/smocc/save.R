# smocc.R
# Fits the broken stick model to the smocc data
# using birth + standard visits as break points (0-2 years)
# Temporary storage in .store
# This script updates: 
#   store/fit_hgt, store/fit_wgt, store/fit_hdc (if refresh == TRUE)
#   data/smocc.rda and data/smocc_bs.rda
# SvB, 11July2016/13oct2016

library("brokenstick")  # >= 0.44
library("clopus")
library("mice")

# Set the refresh flag to TRUE to recalculate the time-consuming 
# brokenstick models for hgt, wgt and hdc
# if refresh == TRUE, then we overwrite store/fit_hgt, store/fit_wgt and
# store/fit_hdc
refresh <- TRUE

# define project
src <- "smocc"
project <- path.expand("~/Package/donordata/donordata")
wd <- file.path(project, "data-raw/smocc")
source(file.path(wd, "R/fetchandstore.R"))
.store <- file.path(wd, "store/")

# load the smocc data files files
fetch(child, time)

# remove records with unknown ages
# reduces number of records from 18586 to 18537
# measured on 1933 children
data <- time
data <- data[!is.na(data$age), ]

## calculate Z-scores
data$hgt.z <- groeivoorspeller::transform_z(data, "hgt")
data$wgt.z <- groeivoorspeller::transform_z(data, "wgt")
data$hdc.z <- groeivoorspeller::transform_z(data, "hdc")
data$bmi <- data$wgt / (data$hgt/100) ^ 2
data$bmi.z <- groeivoorspeller::transform_z(data, "bmi")

# perform broken stick analyses per outcome
# separately for preterm and a term references
knots <- round(c(0, 28/365.25, 56/365.25, 1/4, 1/3, 1/2,
                 7.5/12, 9/12, 11/12, 14/12, 18/12, 2), 4)
Boundary.knots <- c(0, 3)
# NOTE: The visits for SMOCC occured at
# knots <- round(c(0, 1, 2, 3, 6, 9, 12, 15, 18, 24)/12, 4)
# Boundary.knots <- c(0, 3)
# there are no data in smocc for 4m (x5) and 7.5m (x7)
# 
# Match basistakenpakket BTP and SMOCC
# BTP   SMOCC random effect
# 0w    0w  x1
# 4w    4w  x2
# 8w    8w  x3
# 3m    3m  x4
# 4m    --  x5
# 6m    6m  x6
# 7.5m  --  x7
# 9m    9m  x8
# 11m   --  x9
# --    12m 
# 14m   --  x10
# --    15m 
# 18m   18m x11
# 24m   24m x12
# 36m   36m x13  

if (refresh) {
    # time consuming refresh of broken stick models
    yvar <- "hgt.z"
    fit_hgt <- brokenstick(y = data[, yvar], x = data$age, subjid = data$id,
                           knots = knots, boundary = Boundary.knots)
    store(fit_hgt)
    
    yvar <- "wgt.z"
    fit_wgt <- brokenstick(y = data[, yvar], x = data$age, subjid = data$id,
                           knots = knots, boundary = Boundary.knots)
    store(fit_wgt)
    
    #yvar <- "hdc.z"
    #fit_hdc <- brokenstick(y = data[, yvar], x = data$age, subjid = data$id,
    #                       knots = knots, boundary = Boundary.knots)
    #store(fit_hdc)
    
    yvar <- "bmi.z"
    fit_bmi <- brokenstick(y = data[, yvar], x = data$age, subjid = data$id,
                           knots = knots, boundary = Boundary.knots)
    store(fit_bmi)
}

# read models back
fetch(fit_hgt, fit_wgt, fit_hdc, fit_bmi)

hgt <- donordata:::convert2y("hgt", fit_hgt, child, 1)
wgt <- donordata:::convert2y("wgt", fit_wgt, child, 2)
hdc <- donordata:::convert2y("hdc", fit_hdc, child, 1)
bmi <- donordata:::convert2y("bmi", fit_bmi, child, 2)

# determine id (1933)
id <- unique(data$id)
child <- child[child$id %in% id, ]

# combine into donorset
don <- dplyr::bind_cols(child, data.frame(hgt), data.frame(wgt), 
                        data.frame(hdc), data.frame(bmi))

# initial imputation model
don2 <- don
don2$sex <- as.factor(don2$sex)
don2$etn <- as.factor(don2$etn)
don2$edu <- as.factor(don2$edu)
ini <- mice(don2, m = 1, maxit = 0)

# single impute missing values
imp <- mice(don2, m = 1, method = "pmm", seed = 72981)
child2 <- complete(imp)
child <- child2
child$sex <- as.character(child$sex)
child$etn <- as.character(child$etn)
child$edu <- as.character(child$edu)

# overwrite any imputed bs estimates by NA
#idx <- c(grep("hgt_", names(child)), grep("wgt_", names(child)),
#         grep("hdc_", names(child)), grep("bmi_", names(child)))
#bs_names <- names(child)[idx]
#child[, bs_names] <- don[, bs_names]

# add child names of 10 demo children
id.smocc <- c(34071, 34072, 34073, 34075, 34076, 34077, 34078, 34079, 34080, 34081)
idx <- child$id %in% id.smocc
child$name <- ""
child$name[idx] <- c("Laura S", "Thomas S", "Anne S", "Jeroen S",
                     "Mark S", "Kevin S", "Linda S", "Iris S",
                     "Tim S", "Rick S")

# make id consistent with child and time
id <- intersect(unique(child$id), unique(time$id))

# create smocc donordata
smocc <- list(
    child = tibble::as_tibble(child[child$id %in% id, ]),
    time = tibble::as_tibble(data[data$id %in% id, ]))

# number of donor cases: n = 1933

# save 
usethis::use_data(smocc, overwrite = TRUE)


# export the broken stick models
export_hgt <- export(fit_hgt)
export_wgt <- export(fit_wgt)
export_hdc <- export(fit_hdc)
export_bmi <- export(fit_bmi)

# store the combined broken stick model in main data
# directory for lazy loading when attached
smocc_bs <- list(
    hgt = export_hgt,
    wgt = export_wgt,
    hdc = export_hdc,
    bmi = export_bmi)
usethis::use_data(smocc_bs, overwrite = TRUE)

