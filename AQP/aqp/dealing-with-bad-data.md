<!--
	manually process like this:
	library(knitr); library(markdown)
knit('aqp-intro.Rmd')
markdownToHTML('aqp-intro.md', 'aqp-intro.html')
-->
	


## Dealing With Troublesome Data
A quick tutorial on how to search for and accommodate messy data. Inconsistent horizon depths, missing horizons, overlapping horizons, and other such mayhem can lead to unexpected results. It is best to search for and filter (or fix) these kind of errors before proceeding to analyze soil profile data. While classes and methods within the `aqp` package are fairly tolerant of messy data, it is recommended that you apply these tests to your data before feeding into `aqp` functions. This tutorial requires the `plyr` package; if you don't already have it, you can get it like this: `install.packages('plyr', dep=TRUE)`.

### Sample Data
Copy and paste this code into an R session to familiarize yourself with the sample data set used in this tutorial.

```r
library(aqp)

# load sample data set, a simple data.frame object with horizon-level data from 10 profiles
data(sp1)
str(sp1)

# optionally read about it...
# ?sp1

# upgrade to SoilProfileCollection
# 'id' is the name of the column containing the profile ID
# 'top' is the name of the column containing horizon upper boundaries
# 'bottom' is the name of the column containing horizon lower boundaries
depths(sp1) <- id ~ top + bottom

# check it out:
class(sp1)
print(sp1)
plot(sp1)
```


### Data Cleaning
Soil data often contain records where the lower depth of the deepest horizon is missing (`NA`). These data can either be filtered out or "cleaned" by replacing the missing lower depth with the corresponding upper depth. Older soils data with may have O horizons that start above the soil surface (e.g. "Oe 3 to 0 cm). Data containing these type of horizons must either be filtered or fixed before use with classes or methods defined in the `aqp` package. 

```r
# replace missing bottom depths
sp1$bottom[!is.na(sp1$top) & is.na(sp1$bottom)] <- sp1$top[!is.na(sp1$top) & is.na(sp1$bottom)]

# remove O horizons where top > bottom
bad.O.hz.idx <- which(sp1$top > sp1$top)
if(length(bad.O.hz.idx) > 0)
	sp1 <- sp1[-bad.O.hz.idx, ]
```


### Checking for Bad Data
Lets break some good data to make a point.

```r
# load required libraries
library(aqp)
library(plyr)

# load sample data (you already tried this right?)
data(sp1)

# make a copy of the example data set, we will insert errors later
bad <- sp1

# insert a missing horizon boundary in profile P001
bad$top[4] <- NA

# create an overlapping horizon in profile P002
bad$top[9] <- 15

# check:
head(sp1[, c('id', 'top', 'bottom')], 10) # good data
```

```
##      id top bottom
## 1  P001   0      2
## 2  P001   2     14
## 3  P001  14     49
## 4  P001  49     57
## 5  P001  57     89
## 6  P001  89     89
## 7  P002   0      5
## 8  P002   5     19
## 9  P002  19     33
## 10 P002  33     59
```

```r
head(bad[, c('id', 'top', 'bottom')], 10) # bad data
```

```
##      id top bottom
## 1  P001   0      2
## 2  P001   2     14
## 3  P001  14     49
## 4  P001  NA     57
## 5  P001  57     89
## 6  P001  89     89
## 7  P002   0      5
## 8  P002   5     19
## 9  P002  15     33
## 10 P002  33     59
```

```r
# test for missing horizon depths OR overlapping horizons
ddply(sp1, 'id', test_hz_logic, topcol='top', bottomcol='bottom') # TRUE is good
```

```
##     id hz_logic_pass
## 1 P001          TRUE
## 2 P002          TRUE
## 3 P003          TRUE
## 4 P004          TRUE
## 5 P005          TRUE
## 6 P006          TRUE
## 7 P007          TRUE
## 8 P008          TRUE
## 9 P009          TRUE
```

```r
ddply(bad, 'id', test_hz_logic, topcol='top', bottomcol='bottom') # FALSE is bad
```

```
##     id hz_logic_pass
## 1 P001         FALSE
## 2 P002         FALSE
## 3 P003          TRUE
## 4 P004          TRUE
## 5 P005          TRUE
## 6 P006          TRUE
## 7 P007          TRUE
## 8 P008          TRUE
## 9 P009          TRUE
```


### Filtering out Bad Data
The `test_hz_logic()` function is applied to each profile using the `ddply()` function, with results returned as a simple table. We can use this table to determine which profiles should either be filtered or fixed. Note that we previously corrupted the horizons in profiles `P001` and `P002`.

```r
# apply test to our 'bad' data and save result
bad.test <- ddply(bad, 'id', test_hz_logic, topcol='top', bottomcol='bottom')

# which are the good (valid) profiles?
good.profiles <- as.character(bad.test$id[which(bad.test$hz_logic_pass)])

# keep the good ones
good <- subset(bad, id %in% good.profiles)

# promote to SoilProfileCollection
depths(good) <- id ~ top + bottom

# set figure margins to 0
par(mar=c(0,0,0,0))
# plot
plot(good)
```

<img src="figure/filtering.png" title="plot of chunk filtering" alt="plot of chunk filtering" style="display: block; margin: auto;" />


### Concluding Notes
Classes and methods within the `aqp` package can accommodate missing horizons (e.g. a horizon that was not submitted for lab characterization), however the following conditions are not allowed:
  * missing horizon boundaries (bottom of an R horizon)
  * top depth > bottom depth (old-style O horizon)
  * overlapping horizons (data-entry error or multiple samples from overlapping depths)


----------------------------
This document is based on `aqp` version 1.7-4.
