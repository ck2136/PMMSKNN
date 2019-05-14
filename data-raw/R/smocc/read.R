# import.data.R
#
# imports the SMOCC data for curve matching
# this version includes date of birth and data of measurement

library(haven)
library(car)
library(clopus)

# define paths
project <- path.expand("~/Package/donordata/donordata")
wd <- file.path(project, "data-raw/smocc")
setwd(wd)

# Input files
longname   <- file.path(project, "data-raw", "SMOCC", "data", "SMOCCtime.sav")
broadname  <- file.path(project, "data-raw", "SMOCC", "data", "SMOCCchild.sav")

# read the data
smocc.broad <- read_spss(file = broadname)
smocc.long  <- read_spss(file = longname)

# select all visits from the long data
smocc.long <- smocc.long[smocc.long$visit == 1,]

# calculate time varying variables from long
id <- smocc.long$pnr

# data "long" contains 2040 unique children and 16953 rows
length(unique(id))
length(id)

# calculate dates of measurement and dates of birth
dom <- with(smocc.long, ISOdate(1900 + year, month, day))
dob <- with(smocc.long, ISOdate(1900 + k101617, k101415, k101213))

# calculate decimal age
day <- as.numeric(dom - dob)
age <- round(day / 365.25, 4)

# simple data transforms and edits
hgt <- smocc.long$k103336 / 10
hgt[hgt < 35] <- NA
wgt <- smocc.long$k102832/1000
wgt[wgt == 0] <- NA
hdc <- smocc.long$k103739 / 10
hdc[hdc == 0] <- NA

# nutrution indicator define three types of milkfeeding
bf1 <- as.logical(recode(smocc.long$k1041, "1 = TRUE; 0 = FALSE; else = NA"))  # breast
bf2 <- as.logical(recode(smocc.long$k1042, "1 = TRUE; 0 = FALSE; else = NA"))  # humanized
bf3 <- as.logical(recode(smocc.long$k1043, "1 = TRUE; 0 = FALSE; else = NA"))  # cow
bf4 <- as.logical(recode(smocc.long$k1044, "1 = TRUE; 0 = FALSE; else = NA"))  # other
bf5 <- as.logical(recode(smocc.long$k1045, "1 = TRUE; 0 = FALSE; else = NA"))  # remarks

# define three types of milkfeeding
breast      <- bf1 & !(bf2 | bf3 | bf4)
mixed       <- bf1 &  (bf2 | bf3 | bf4)
bottle      <- (!bf1) & (bf2 | bf3 | bf4)

# save dates in Dutch format
dob <- format(as.Date(dob), "%d-%m-%y")
dom <- format(as.Date(dom), "%d-%m-%y")

long <- data.frame(src = "smocc", id = id,
                   dob, dom, age, hgt, wgt, hdc,
				   breast, mixed, bottle, 
				   stringsAsFactors = FALSE)

# aggregate to find duration of exc and mixed breastfeeding
# in: id, age, breast, mixed, bottle
# out: bfexcl06, durbrst 

## FIXME: still need to do, calculate bfexcl06 and durbrst

agg1 <- aggregate(long[, c("id", "age", "breast", "mixed")], 
                 by = list(long$id, long$breast), 
                 FUN = max, na.rm = TRUE)
agg1 <- agg1[agg1$Group.2, c("id", "age")]
agg1$age <- round(365.25 * agg1$age)
agg1$age[agg1$age > 182] <- 182

agg2 <- aggregate(long[, c("id", "age", "breast", "mixed")], 
                 by = list(long$id, long$mixed), 
                 FUN = max, na.rm = TRUE)
agg2 <- agg2[agg2$Group.2, c("id", "age")]
agg2$age <- round(365.25 * agg2$age)
agg2$age[agg2$age > 182] <- 182

bf <- smocc.broad[, "pnr"]
bf$id <- bf$pnr
bf1 <- merge(bf, agg1, by = "id", all.x = TRUE)
bf2 <- merge(bf1, agg2, by = "id", all.x = TRUE)

## find birth wgt, birth hgt and birth hdc in broad data
bw <- smocc.broad$k061417
wgt.0 <- bw / 1000
hgt.0 <- smocc.broad$k062223
hgt.0[hgt.0 == 0] <- NA
hdc.0 <- smocc.broad$k062425
hdc.0[hdc.0 == 0] <- NA

# add birth measurements to long data
dob <- with(smocc.broad, ISOdate(1900 + k011617, k011415, k011213))
dob <- format(as.Date(dob), "%d-%m-%y")

birth <- data.frame(src = "smocc", id = smocc.broad$pnr,
                    dob = dob, dom = dob,
                    age = 0, hgt = hgt.0, wgt = wgt.0, hdc = hdc.0,
					breast = NA, mixed = NA, bottle = NA,
					stringsAsFactors = FALSE)
long <- rbind(birth, long)


# calculate record number, sorted according to id and age
idx <- order(long$id, long$age)
mrec <- table(long[idx, "id"])
rec <- unlist(sapply(mrec, seq)) # record number
nrec <- rep(mrec, times = mrec)  # number of records

# add sex, etnicity and gestational age for calculating Z-scores
sex <- as.character(recode(as_factor(smocc.broad$k0118), "'Jongen' = 'male'; 'Meisje' = 'female'"))
sexrep <- rep(sex, times = mrec)
etn <- as.character(recode(as_factor(smocc.broad$k026971), "'TURKIJE' = 'TU'; 'MAROKKO' = 'MA'; else = 'NL'"))
etnrep <- rep(etn, times = mrec)
ga <- smocc.broad$k051213
ga[ga == 0] <- NA
garep <- rep(ga, times = mrec)
bw <- smocc.broad$k061417
bwrep <- rep(bw, times = mrec)

# reorder long data set
long <- long[idx, ]
long <- data.frame(src = long$src, id = long$id,
                   rec = rec, nrec = nrec,
                   dob = long$dob, dom = long$dom,
                   age = long$age, sex = sexrep,
                   etn = etnrep, ga = garep, bw = bwrep,
                   hgt = long$hgt, wgt = long$wgt, hdc = long$hdc,
                   stringsAsFactors = FALSE)

# compare timing of visits to current schedule
brk <- round(c(0, 28/365.25, 56/365.25, 1/4, 1/3, 1/2, 7.5/12,
               9/12, 11/12, 14/12, 18/12, 2, 3, 3+9/12), 4)
plot(long$age, long$wgt)
abline(v = brk, lwd = 3, col = "red", lty = 2)


# calculate Z-scores w.r.t national references
data <- long

# eliminate some values that were missed during data cleaning
data[data$id == 15123 & data$rec == 1 , c("hgt")] <- NA 
data[data$id == 15194 & data$rec == 5 , c("hgt")] <- NA 
data[data$id == 35027 & data$rec == 6 , c("hgt","wgt","hdc")] <- NA 

# hgt.z <- y2z(y = data$hgt, x = data$age, ref = nl1997, yname = "hgt",
#              sex = data$sex, sub = data$etn, drop = TRUE)
# wgt.z <- y2z(y = data$wgt, x = data$age, ref = nl1997, yname = "wgt",
#              sex = data$sex, sub = data$etn, drop = TRUE)
# hdc.z <- y2z(y = data$hdc, x = data$age, ref = nl1997, yname = "hdc",
#              sex = data$sex, sub = data$etn, drop = TRUE)
# calculate Z-scores w.r.t preterm references (ga < 37)
# hgt.z.pt <- y2z(y = data$hgt, x = data$age, ref = preterm, yname = "hgt",
#                 sex = data$sex, sub = data$ga, drop = TRUE)
# wgt.z.pt <- y2z(y = data$wgt, x = data$age, ref = preterm, yname = "wgt",
#                 sex = data$sex, sub = data$ga, drop = TRUE)
# hdc.z.pt <- y2z(y = data$hdc, x = data$age, ref = preterm, yname = "hdc",
#                 sex = data$sex, sub = data$ga, drop = TRUE)

# use preterm z-score for those ga < 37
# idx <- !is.na(data$ga) & data$ga < 36 
# hgt.z[idx] <- hgt.z.pt[idx]
# wgt.z[idx] <- wgt.z.pt[idx]
# hdc.z[idx] <- hdc.z.pt[idx]

long2 <- data.frame(data)

# only use children with at least 6 out of 11 visits
longitudinal <- long2[long2$nrec >= 6, ]

# extract relevant background characteristics
# health of the child,
# smoking of the mother during pregnancy, education of
# father and mother (low versus middle/high),
# age of the mother, and height of father and mother

id <- smocc.broad$pnr
twin <- smocc.broad$k0612 - 1
goodhealth <- recode(smocc.broad$k3072, "1=1; c(2,3,4)=0; else=NA")
smo <- recode(smocc.broad$k0540, "1=0; c(2,3,4,5)=1; else=NA")
eduf <- as.character(recode(as_factor(smocc.broad$k0246), "c('Blo','Lo','Lbo') = 'low';
               c('Mavo','Mbo','Havo','Vwo') = 'middle';
               c('Hbo','Univ.') = 'high';
               else = NA"))
edum <- as.character(recode(as_factor(smocc.broad$k0273), "c('Blo','Lo','Lbo') = 'low';
               c('Mavo','Mbo','Havo','Vwo') = 'middle';
               c('Hbo','Univ.') = 'high';
               else = NA"))
edu  <- edum


# age of mother
bdm <- smocc.broad$k026364
bdm[bdm == 0 | bdm == 99 | is.na(bdm)] <- 1
bmm <- smocc.broad$k026566
bmm[bmm == 0 | bmm == 99 | is.na(bmm)] <- 1
bym <- smocc.broad$k026768
bym[bym == 0] <- NA
daym <- as.numeric(with(smocc.broad,
                        ISOdate(1900 + k011617, k011415, k011213) -
                            ISOdate(1900 + bym, bmm, bdm)))
agem <- trunc(daym / 365.25)

hgtm <- smocc.broad$k032123
hgtm[hgtm == 0] <- NA
hgtf <- smocc.broad$k026062
hgtf[hgtf == 0] <- NA


background <- data.frame(src = "smocc", id,
                         dob = dob,
                         sex, etn, edu, 
                         ga, bw, twin, 
                         agem, smo, hgtm, hgtf,
                         stringsAsFactors = FALSE)

child <- background
time <- longitudinal

# save
save(child, file = file.path(wd, "store","child"))
save(time, file = file.path(wd, "store", "time"))

## now run save.R