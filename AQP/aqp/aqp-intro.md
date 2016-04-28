---
output:
  html_document:
    mathjax: null
    jquery: null
    smart: no
---



## Introduction to SoilProfileCollection Objects
This is a very basic introduction to the `SoilProfileCollection` class object defined in the `aqp` package for **R**. The `SoilProfileCollection` class was designed to simplifiy the process of working with the collection of data associated with soil profiles: site-level data, horizon-level data, spatial data, diagnostic horizon data, metadata, etc. Examples listed below are meant to be copied/pasted from this document and interactively run within **R**. Comments (green text) briefly describe what the code in each line does. This document assumes a basic level of proficiency with **R** which can be gained by reviewing some of the material in [tutorials like this](http://cran.r-project.org/doc/contrib/Verzani-SimpleR.pdf). Further documentation on objects and functions from the `aqp` package can be accessed by typing `help(aqp)` (or more generally, `?function_name`) at the **R** console.


### Object Creation
`SoilProfileCollection` objects are typically created by "promoting" `data.frame` objects (rectangular tables of data) that contain at least three essential columns: 
 * 1) an ID column uniquely identifying groups of horizons (e.g. pedons)
 * 2) horizon top boundaries
 * 3) horizon bottom boundaries

The `data.frame` should be pre-sorted according to the profile ID and horizon top boundary. Formula notation is used to define the columns used to promote a `data.frame` object:


```r
idcolumn ~ hz_top_column + hz_bottom_column
```

In this tutorial we will use some sample data included with the `aqp` package, based on characterization data from 10 soils sampled on serpentinitic parent material as described in [McGahan et al, 2009](https://www.soils.org/publications/sssaj/articles/73/6/2087).

```r
# load required packages, you may have to install these if missing:
# install.packages('aqp', dep=TRUE)
library(aqp)
library(Hmisc)
library(lattice)
library(MASS)

# load sample data set, a data.frame object with horizon-level data from 10 profiles
data(sp4)
str(sp4)
```

```
## 'data.frame':	30 obs. of  13 variables:
##  $ id         : chr  "colusa" "colusa" "colusa" "colusa" ...
##  $ name       : chr  "A" "ABt" "Bt1" "Bt2" ...
##  $ top        : int  0 3 8 30 0 9 0 4 13 0 ...
##  $ bottom     : int  3 8 30 42 9 34 4 13 40 6 ...
##  $ K          : num  0.3 0.2 0.1 0.1 0.2 0.3 0.2 0.6 0.8 0.4 ...
##  $ Mg         : num  25.7 23.7 23.2 44.3 21.9 18.9 12.1 12.1 17.7 16.4 ...
##  $ Ca         : num  9 5.6 1.9 0.3 4.4 4.5 1.4 7 4.4 24.1 ...
##  $ CEC_7      : num  23 21.4 23.7 43 18.8 27.5 23.7 18 20 31.1 ...
##  $ ex_Ca_to_Mg: num  0.35 0.23 0.08 0.01 0.2 0.2 0.58 0.51 0.25 1.47 ...
##  $ sand       : int  46 42 40 27 54 49 43 36 27 43 ...
##  $ silt       : int  33 31 28 18 20 18 55 49 45 42 ...
##  $ clay       : int  21 27 32 55 25 34 3 15 27 15 ...
##  $ CF         : num  0.12 0.27 0.27 0.16 0.55 0.84 0.5 0.75 0.67 0.02 ...
```

```r
# optionally read about it...
# ?sp4

# upgrade to SoilProfileCollection
# 'id' is the name of the column containing the profile ID
# 'top' is the name of the column containing horizon upper boundaries
# 'bottom' is the name of the column containing horizon lower boundaries
depths(sp4) <- id ~ top + bottom

# check it out:
class(sp4) # class name
```

```
## [1] "SoilProfileCollection"
## attr(,"package")
## [1] "aqp"
```

```r
print(sp4)
```

```
## Object of class SoilProfileCollection
## Number of profiles: 10
## Depth range: 16-49 cm
## 
## Horizon attributes:
##       id name top bottom   K   Mg  Ca CEC_7 ex_Ca_to_Mg sand silt clay   CF
## 1 colusa    A   0      3 0.3 25.7 9.0  23.0        0.35   46   33   21 0.12
## 2 colusa  ABt   3      8 0.2 23.7 5.6  21.4        0.23   42   31   27 0.27
## 3 colusa  Bt1   8     30 0.1 23.2 1.9  23.7        0.08   40   28   32 0.27
## 4 colusa  Bt2  30     42 0.1 44.3 0.3  43.0        0.01   27   18   55 0.16
## 5  glenn    A   0      9 0.2 21.9 4.4  18.8        0.20   54   20   25 0.55
## 6  glenn   Bt   9     34 0.3 18.9 4.5  27.5        0.20   49   18   34 0.84
## 
## Sampling site attributes:
##          id
## 1    colusa
## 2     glenn
## 3     kings
## 4  mariposa
## 5 mendocino
## 6      napa
```


### Accessing, Setting, and Replacing Data
"Accessor" functions are used to extract specfic components from within `SoilProfileCollection` objects.

```r
# methods for object inspection
idname(sp4) # self-explanitory
horizonDepths(sp4) # self-explanitory
depth_units(sp4) # defaults to 'cm'
metadata(sp4) # not much to start with
profile_id(sp4) # vector of profile IDs
horizonNames(sp4) # column names from horizon data
siteNames(sp4) # column names from site data

# further inspection with common function overloads
length(sp4) # number of profiles in the collection
nrow(sp4) # number of horizons in the collection
names(sp4) # column names from site and horizon data, concatenated into a single vector
min(sp4) # shallowest profile depth in collection
max(sp4) # deepest profile depth in collection
```


Several accessor functions can also be used to define, modify, or replace components within `SoilProfileCollection` objects. For example, 

```r
horizons(sp4)
```
will extract horizon data and 


```r
# careful!
horizons(sp4) <- x
```

will replace horizon data. Further examples include:

```r
# extraction of soil profile components
site(sp4)  # get or set site data
diagnostic_hz(sp4) # get or set diagnostic horizons
proj4string(sp4)
coordinates(sp4)
```


### Horizon and Site Data
Typically, horizon and site data are the most important components of `SoilProfileCollection` objects. Both are internally stored as `data.frame` objects; with one or more rows (per profile ID) in the horizon table and one row (per profile ID) in the site table. Columns from either table can be accessed with the `$` style notation of `data.frames`. New data can be assigned to either table in the same manner, as long as the length of the new data is either: 
 * 1) same length as the number of profiles in the collection (target is the site table)
 * 2) same length as the number of horizons in the collection (target is the horizon table)
 

```r
# assignment of new data to existing or new attributes
sp4$elevation <- rnorm(n=length(sp4), mean=1000, sd=150) # site-level, based on length of assigned data
sp4$thickness <- sp4$bottom - sp4$top # horizon-level

# extraction of specific attributes by name
sp4$clay # vector of clay content (horizon data)
```

```
##  [1] 21 27 32 55 25 34  3 15 27 32 25 31 33 13 21 23 15 17 12 19 14 14 22 25 40 51 67 24 25 32
```

```r
sp4$elevation # vector of simulated elevation (site data)
```

```
##  [1] 1171.8250 1020.4769 1142.4566 1470.1071  952.8598  905.9533 1070.3719  926.5812 1067.3850
## [10] 1092.5689
```

```r
# assign a single single value into horizon-level attributes
sp4$constant <- rep(1, times=nrow(sp4))

# promote horizon-level data to site-level data (when it makes sense to do so)
# note that this _moves_ the named column from horizon to site
site(sp4) <- ~ constant 
```


Horizon and site data can also be modified via extraction to `data.frame` followed by replacement (horizon data) or join (site data). Note that while this approach gives the most flexibility, it is also the most dangerous-- replacement of horizon data with new data that don't exactly conform to the original sorting may corrupt your `SoilProfileCollection`.

```r
# extract horizon data to data.frame
h <- horizons(sp4)

# add a new column and save back to original object
h$random.numbers <- rnorm(n=nrow(h), mean=0, sd=1)

# _replace_ original horizon data with modified version
# ! row-order should not be altered !
horizons(sp4) <- h

# extract site data to data.frame
s <- site(sp4)

# add a fake group to the site data
s$group <- rep(c('A', 'B'), length.out=nrow(s))

# join new site data with previous data: old data are _not_ replaced
site(sp4) <- s

# check:
sp4
```

```
## Object of class SoilProfileCollection
## Number of profiles: 10
## Depth range: 16-49 cm
## 
## Horizon attributes:
##       id name top bottom   K   Mg  Ca CEC_7 ex_Ca_to_Mg sand silt clay   CF thickness
## 1 colusa    A   0      3 0.3 25.7 9.0  23.0        0.35   46   33   21 0.12         3
## 2 colusa  ABt   3      8 0.2 23.7 5.6  21.4        0.23   42   31   27 0.27         5
## 3 colusa  Bt1   8     30 0.1 23.2 1.9  23.7        0.08   40   28   32 0.27        22
## 4 colusa  Bt2  30     42 0.1 44.3 0.3  43.0        0.01   27   18   55 0.16        12
## 5  glenn    A   0      9 0.2 21.9 4.4  18.8        0.20   54   20   25 0.55         9
## 6  glenn   Bt   9     34 0.3 18.9 4.5  27.5        0.20   49   18   34 0.84        25
##   random.numbers
## 1    0.374283715
## 2    0.349707479
## 3   -1.117226780
## 4    0.000831656
## 5   -0.360970422
## 6    0.692707609
## 
## Sampling site attributes:
##          id elevation constant group
## 1    colusa 1171.8250        1     A
## 2     glenn 1020.4769        1     B
## 3     kings 1142.4566        1     A
## 4  mariposa 1470.1071        1     B
## 5 mendocino  952.8598        1     A
## 6      napa  905.9533        1     B
```


### Diagnostic Horizons
Diagnostic horizons typically span several genetic horizons and may or may not be present in all profiles. To accomodate the wide range of possibilities, diagnostic horizon data are stored as a `data.frame` in "long format": each row corresponds to a diagnostic horizon, identified with a column matching the ID column used to initialize the `SoilProfileCollection` object.

```r
# manually create some diagnostic horizon data
# there is no restrictions on data format, as long as each row has an ID that exists within the collection
# be sure to use the ID column name that was used to initialize the SoilProfileCollection object
# check via: idname(sp4)
dh <- data.frame(id='colusa', kind='argillic', top=8, bottom=42, stringsAsFactors=FALSE)

# overwrite any existing diagnostic horizon data
diagnostic_hz(sp4) <- dh

# append to diagnostic horizon data
dh <- diagnostic_hz(sp4)
dh.new <- data.frame(id='napa', kind='argillic', top=6, bottom=20, stringsAsFactors=FALSE)

# overwrite existing diagnostic horizon data with appended data
diagnostic_hz(sp4) <- rbind(dh, dh.new)
```


### Spatial Data
Spatial data can be explicitly stored within a `SoilProfileCollection` object and accessed with methods imported from the [sp](http://cran.at.r-project.org/web/packages/sp/index.html) package. The use of `sp` classes (in this case `SpatialPoints` and `SpatialPointsDataFrame` objects) simplifies operations such as [plotting spatial data, coordinate system transformations, and spatial queries](http://cran.at.r-project.org/web/packages/sp/vignettes/intro_sp.pdf).

```r
# generate some fake coordinates as site level attributes
sp4$x <- rnorm(n=length(sp4), mean=354000, sd=100)
sp4$y <- rnorm(n=length(sp4), mean=4109533, sd=100)

# initialize spatial coordinates
coordinates(sp4) <- ~ x + y

# extract coordinates as matrix
coordinates(sp4)
```



|        x|       y|
|--------:|-------:|
| 353944.9| 4109535|
| 353930.9| 4109615|
| 353905.0| 4109617|
| 354078.0| 4109628|
| 353908.9| 4109632|
| 353872.6| 4109552|
| 354004.5| 4109608|
| 354024.0| 4109634|
| 353896.5| 4109535|
| 353908.6| 4109614|

```r
# get/set spatial reference system using PROJ4 syntax
proj4string(sp4) <- '+proj=utm +zone=11 +datum=NAD83'
proj4string(sp4)
```

```
## [1] "+proj=utm +zone=11 +datum=NAD83"
```

```r
# extract spatial data + site level attribtutes
# see ?SpatialPointsDataFrame for details
sp4.sp <- as(sp4, 'SpatialPointsDataFrame')

# plot the fake coordinates
plot(sp4.sp)
box()
```

<img src="figure/spatial-1.png" title="plot of chunk spatial" alt="plot of chunk spatial" style="display: block; margin: auto;" />


### Subsetting SoilProfileCollection Objects
`SoilProfileCollection` objects can be subset using the familiar `[`-style notiation used by matrix and `data.frame` objects, such that: `spc[i, j]` will return profiles identified by the integer vector `i`, and horizons identified by the integer vector `j`. Ommitting either index will result in all profiles (`i` ommitted) or all horizons (`j` ommitted). Typically, site-level attributes will be used as the subsetting criteria. Functions that return an index to matches (such as `grep()` or `which()`) provide the link between attributes and an index to matching profiles. Some examples:


```r
# explicit string matching
idx <- which(sp4$group == 'A')

# numerical expressions
idx <- which(sp4$elevation < 1000)

# regular expression, matches any profile ID containing 'shasta'
idx <- grep('shasta', profile_id(sp4), ignore.case=TRUE)

# perform subset based on index
sp4[idx, ]
```


The `subsetProfiles()` function can be used to subset collections based on criteria from both horizon and site level attributes. Several examples are given below.

```r
# matrix-style subsetting
sp4[1:2, ] # profiles 1--2 from the collection
sp4[, 1:2] # horizons 1--2 from each profile

# in the presence of spatial data, subsetting can either result in 
# 1. a new SoilProfileCollection (> 1 horizon / profile)
# 2. a new SpatialPointsDataFrame (1 horizon / profile)
class(sp4[1, ]) # profile 1 from the collection
class(sp4[, 1]) # horizon 1 from each profile

# filtering criteria from horizon and site level data
subsetProfiles(sp4, h='clay < 20', s='group == "A"')
```


### Concatenation of SoilProfileCollection Objects
If the internal structures of several `SoilProfileCollection` objects are identical (e.g. same ID name, horizon, site, diagnostic horizon, and spatial attributes) they can be safely combined into a new object. This is typically only possible when re-combining subsets of the same original object. Note that concatenation will re-order all input `SoilProfileCollection` objects such that the resulting IDs are in alpha-numeric order.

```r
# subset data into chunks
s1 <- sp4[1:2, ]
s2 <- sp4[4, ]
s3 <- sp4[c(6, 8, 9), ]

# re-combine: think "row-bind"
s <- rbind(s1, s2, s3)
```


### Plotting SoilProfileCollection Objects
The `plot()` method for `SoilProfileCollection` objects generates sketches of profiles within the collection based on horizon boundaries, vertically aligned to an integer sequence from 1 to the number of profiles. Horizon names are automatically extracted from a horizon-level attribute `name` (if present), or via an alternate attributed given as an argument: `name='column.name'`. Horizon colors are automatically generated from the horizon-level attribute `soil_color`, or any other attribute of **R**-compatible color description given as an argument: `color='column.name'`. This function is highly customizable, therefore, it is prudent to consult `help(plotSPC)` from time to time. Soil colors in Munsell notation can be converted to **R**-compatible colors via `munsell2rgb()`.

```r
# plot uses base graphics, so it is simple to add low-level annotation layers
plot(sp4, name='name')

# inspect plotting area, very simple to overlay graphical elements
abline(v=1:length(sp4), lty=3, col='blue')

# profiles are centered at integers, from 1 to length(obj)
axis(1, line=0, at=1:10, cex.axis=0.75, font=4, col='blue', col.axis='blue', lwd=2)
mtext('profile axis', side=1, line=2, font=4, col='blue')

# y-axis is based on profile depths
axis(2, line=-1, at=pretty(1:max(sp4)), cex.axis=0.75, font=4, las=1, col='blue', col.axis='blue', lwd=2)
mtext('depth axis', side=2, line=1, font=4, col='blue')
```

<img src="figure/plotting-1.png" title="plot of chunk plotting" alt="plot of chunk plotting" style="display: block; margin: auto;" />


Horizon-level attributes can be symbolized with color, in this case using the horizon-level attribute "clay":

```r
# plot again, this time using new colors
par(mar=c(0,0,3,0)) # tighter figure margins
plot(sp4, name='name', color='clay')
```

<img src="figure/plotting-with-color-1.png" title="plot of chunk plotting-with-color" alt="plot of chunk plotting-with-color" style="display: block; margin: auto;" />

Horizon-level attributes that represent a volume fraction (e.g. coarse-fragment percentage) can be added to an existing figure:

```r
par(mar=c(0,0,3,0)) # tighter figure margins
sp4$frag_pct <- sp4$CF * 100 # convert fraction -> percent
plot(sp4, name='name', color='frag_pct')
addVolumeFraction(sp4, 'frag_pct') # symbolize volume fraction data
```

<img src="figure/plotting-vol-fraction-1.png" title="plot of chunk plotting-vol-fraction" alt="plot of chunk plotting-vol-fraction" style="display: block; margin: auto;" />

Annotation of depth-intervals can be accomplished using:

```r
# generate some fake depth brackets:
fake.tops <- rep(0, times=length(sp4))
fake.bottoms <- runif(n=length(sp4), min=0, max=20)
par(mar=c(0,0,3,0)) # tighter figure margins
plot(sp4, name='name')
addBracket(fake.tops, fake.bottoms, col='red') # add depth brackets
```

<img src="figure/plotting-depth-interval-1.png" title="plot of chunk plotting-depth-interval" alt="plot of chunk plotting-depth-interval" style="display: block; margin: auto;" />


### Iterating Over Profiles in a Collection
The `profileApply()` function is an extension of the familiar *apply()* family of functions that operate on vectors (`sapply` and `tapply`), matrices (`apply`), and lists (`lapply`)-- extended to `SoilProfileCollection` objects. The function named in the `FUN` argument is evaluated once for each profile in the collection, typicially returning a single value per profile. In this case, the ordering of the results would match the ordering of values in the site level attribute table.

```r
# max() returns the depth of a soil profile
sp4$soil.depth <- profileApply(sp4, FUN=max)
# max() with additional argument give max depth to non-missing 'clay'
sp4$soil.depth.clay <- profileApply(sp4, FUN=max, v='clay')
# nrow() returns the number of horizons
sp4$n.hz <- profileApply(sp4, FUN=nrow)
# compute the mean clay content by profile using an inline function
sp4$mean.clay <- profileApply(sp4, FUN=function(i) mean(i$clay))
# estimate soil depth based on horizon designation
sp4$soil.depth <- profileApply(sp4, estimateSoilDepth, name='name', top='top', bottom='bottom')
```

When `FUN` returns a vector of the same length as the number of horizons in a profile, `profileApply()` can be used to create new horizon-level attributes. For example, the change in clay content with depth could be calculated via:

```r
# save as horizon-level attribute
sp4$delta.clay <- profileApply(sp4, FUN=function(i) c(NA, diff(i$clay)))

# check results:
horizons(sp4)[1:6, c('id', 'top', 'bottom', 'clay', 'delta.clay')]
```



|id     | top| bottom| clay| delta.clay|
|:------|---:|------:|----:|----------:|
|colusa |   0|      3|   21|         NA|
|colusa |   3|      8|   27|          6|
|colusa |   8|     30|   32|          5|
|colusa |  30|     42|   55|         23|
|glenn  |   0|      9|   25|         NA|
|glenn  |   9|     34|   34|          9|


More complex summaries can be generated by writing a custom function that is then called by `profileApply()`. Note that each profile is passed into this function and accessed via a temporary variable (`i`), which is a `SoilProfileCollection` object containing a single profile. A list of `SoilProfileCollection` objects returned from a custom function can be re-constituted into a single `SoilProfileCollection` object via `rbind`. See `help(profileApply)` for details.

```r
# compute hz-thickness weighted mean exchangeable-Ca:Mg
wt.mean.ca.mg <- function(i) {
	# use horizon thickness as a weight
	thick <- i$bottom - i$top
	# compute the weighted mean, accounting for the possibility of missing data
	m <- wtd.mean(i$ex_Ca_to_Mg, weights=thick, na.rm=TRUE)
	return(m)
	}

# apply our custom function and save results as a site-level attribute
sp4$wt.mean.ca.to.mg <- profileApply(sp4, wt.mean.ca.mg)
```

We can now use our some of our new site-level attributes to order the profiles when plotting. In this case profiles are ordered based on the horizon-thickness weighted mean, exchageable Ca:Mg values. Horizons are colored by exchangeable Ca:Mg values.

```r
# the result is an index of rank
new.order <- order(sp4$wt.mean.ca.to.mg)
# plot the data using our new order based on Ca:Mg
par(mar=c(4,0,3,0)) # tighten figure margins
plot(sp4, name='name', color='ex_Ca_to_Mg', plot.order=new.order)
# add an axis labeled with the sorting criteria
axis(1, at=1:length(sp4), labels=round(sp4$wt.mean.ca.to.mg, 3), cex.axis=0.75)
mtext(1, line=2.25, text='Horizon Thickness Weighted Mean Ex. Ca:Mg', cex=0.75)
```

<img src="figure/profileApply-4-1.png" title="plot of chunk profileApply-4" alt="plot of chunk profileApply-4" style="display: block; margin: auto;" />


### Slicing Soil Profile Collections
Collections of soil profiles can be sliced (or re-sampled) into regular depth-intervals with the `slice()` function. The slicing structure and variables of interest are defined via formula notation:


```r
# slice select horizon-level attributes
seq ~ var.1 + var.2 + var.3 + ...
# slice all horizon-level attributes
seq ~ .
```

where `seq` is a sequence of integers (e.g. `0:15`, `c(5,10,15,20)`, etc.) and `var.1 + var.2 + var.3 + ...` are horizon-level attributes to slice. Both continuous and categorical variables can be named on the right-hand-side of the formula. The results returned by `slice()` is either a `SoilProfileCollection`, or `data.frame` when called with the optional argument `just.the.data=TRUE`. For example, to slice our sample data set into 1-cm intervals, from 0--15 cm:



```r
# slice data
sliced <- slice(sp4, fm= 0:15 ~ sand + silt + clay + name + ex_Ca_to_Mg)
# check the result
class(sliced)
```

```
## [1] "SoilProfileCollection"
## attr(,"package")
## [1] "aqp"
```

```r
# plot sliced data
par(mar=c(0,0,3,0)) # tighten figure margins
plot(sliced, name='name', color='ex_Ca_to_Mg')
```

<img src="figure/slice-1.png" title="plot of chunk slice" alt="plot of chunk slice" style="display: block; margin: auto;" />


Once soil profile data have been sliced, it is simple to extract "chunks" of data by depth interval via matrix-style subsetting:


```r
# slice from 0 to max depth in the collection
sliced <- slice(sp4, fm= 0:max(sp4) ~ sand + silt + clay + name + ex_Ca_to_Mg)
# extract all data over the range of 5--10 cm:
plot(sliced[, 5:10])
# extract all data over the range of 25--50 cm:
plot(sliced[, 25:50])
# extract all data over the range of 10--20 and 40--50 cm:
plot(sliced[, c(10:20, 40:50)])
```


### Aggregating Soil Profile Collections Along Regular "Slabs"
Depth-wise summary of horizon-level attributes is performed with the `slab()` function. Profile grouping criterial and horizon attribute selection is parametrized via formula: either `groups ~ var1 + var2 + var3` where named variables are aggregated within `groups` OR where named variables are aggregated across the entire collection ` ~ var1 + var2 + var3`. The default summary function (`slab.fun`) computes the 5th, 25th, 50th, 75th, and 95th percentiles via Harrell-Davis quantile estimator.

The depth structure ("slabs") over which summaries are computed is defined with the `slab.structure` argument using:

  * a single integer (e.g. `10`): data are aggregated over a regular sequence of 10-unit thickness slabs
  * a vector of 2 integers (e.g. `c(50, 60)`): data are aggregated over depths spanning 50--60 units
  * a vector of 3 or more integers (e.g. `c(0, 5, 10, 50, 100)`): data are aggregated over the depths spanning 0--5, 5--10, 10--50, 50--100 units


```r
# aggregate a couple of the horizon-level attributes, 
# across the entire collection, 
# from 0--10 cm
# computing the mean value ignoring missing data
slab(sp4, fm= ~ sand + silt + clay, slab.structure=c(0,10), slab.fun=mean, na.rm=TRUE)
```



|variable | all.profiles| value| top| bottom| contributing_fraction|
|:--------|------------:|-----:|---:|------:|---------------------:|
|sand     |            1| 47.63|   0|     10|                     1|
|silt     |            1| 31.15|   0|     10|                     1|
|clay     |            1| 21.11|   0|     10|                     1|

```r
# again, this time within groups defined by a site-level attribute:
slab(sp4, fm= group ~ sand + silt + clay, slab.structure=c(0,10), slab.fun=mean, na.rm=TRUE)
```



|variable |group | value| top| bottom| contributing_fraction|
|:--------|:-----|-----:|---:|------:|---------------------:|
|sand     |A     | 48.26|   0|     10|                     1|
|silt     |A     | 31.52|   0|     10|                     1|
|clay     |A     | 20.30|   0|     10|                     1|
|sand     |B     | 47.00|   0|     10|                     1|
|silt     |B     | 30.78|   0|     10|                     1|
|clay     |B     | 21.92|   0|     10|                     1|

```r
# again, this time over several depth ranges
slab(sp4, fm= ~ sand + silt + clay, slab.structure=c(0,10,25,40), slab.fun=mean, na.rm=TRUE)
```



|variable | all.profiles|    value| top| bottom| contributing_fraction|
|:--------|------------:|--------:|---:|------:|---------------------:|
|sand     |            1| 47.63000|   0|     10|             1.0000000|
|sand     |            1| 42.38931|  10|     25|             0.8733333|
|sand     |            1| 32.14607|  25|     40|             0.5933333|
|silt     |            1| 31.15000|   0|     10|             1.0000000|
|silt     |            1| 29.41221|  10|     25|             0.8733333|
|silt     |            1| 31.34831|  25|     40|             0.5933333|
|clay     |            1| 21.11000|   0|     10|             1.0000000|
|clay     |            1| 28.10687|  10|     25|             0.8733333|
|clay     |            1| 36.26966|  25|     40|             0.5933333|

```r
# again, this time along 1-cm slices, computing quantiles
agg <- slab(sp4, fm= ~ Mg + Ca + ex_Ca_to_Mg + CEC_7 + clay)

# see ?slab for details on the default aggregate function
head(agg)
```



|variable | all.profiles|     p.q5|    p.q25|    p.q50|    p.q75|    p.q95| top| bottom| contributing_fraction|
|:--------|------------:|--------:|--------:|--------:|--------:|--------:|---:|------:|---------------------:|
|Mg       |            1| 4.168501| 10.78693| 15.37732| 22.22149| 27.71002|   0|      1|                     1|
|Mg       |            1| 4.168501| 10.78693| 15.37732| 22.22149| 27.71002|   1|      2|                     1|
|Mg       |            1| 4.174176| 11.32581| 19.02084| 25.60437| 28.04662|   2|      3|                     1|
|Mg       |            1| 4.254672| 12.52545| 20.31176| 25.89652| 32.56386|   3|      4|                     1|
|Mg       |            1| 4.254672| 12.52545| 20.31176| 25.89652| 32.56386|   4|      5|                     1|
|Mg       |            1| 4.254673| 12.52694| 20.46020| 27.00139| 32.88077|   5|      6|                     1|

```r
# plot median +/i bounds defined by the 25th and 75th percentiles
# this is lattice graphics, syntax is a little rough
xyplot(top ~ p.q50 | variable, data=agg, ylab='Depth',
			 xlab='median bounded by 25th and 75th percentiles',
			 lower=agg$p.q25, upper=agg$p.q75, ylim=c(42,-2),
			 panel=panel.depth_function,
			 alpha=0.25, sync.colors=TRUE,
			 par.settings=list(superpose.line=list(col='RoyalBlue', lwd=2)),
			 prepanel=prepanel.depth_function,
			 cf=agg$contributing_fraction, cf.col='black', cf.interval=5, 
			 layout=c(5,1), strip=strip.custom(bg=grey(0.8)),
			 scales=list(x=list(tick.number=4, alternating=3, relation='free'))
			 )
```

<img src="figure/slab-1.png" title="plot of chunk slab" alt="plot of chunk slab" style="display: block; margin: auto;" />


Depth-wise aggregation can be useful for visual evaluation of multivariate similarity among groups of profiles.

```r
# processing the "napa" and tehama profiles
idx <- which(profile_id(sp4) %in% c('napa', 'tehama'))
napa.and.tehama <- slab(sp4[idx, ], fm= ~ Mg + Ca + ex_Ca_to_Mg + CEC_7 + clay)

# combine with the collection-wide aggregate data
g <- make.groups(collection=agg, napa.and.tehama=napa.and.tehama)

# compare graphically:
xyplot(top ~ p.q50 | variable, groups=which, data=g, ylab='Depth',
			 xlab='median bounded by 25th and 75th percentiles',
			 lower=g$p.q25, upper=g$p.q75, ylim=c(42,-2),
			 panel=panel.depth_function,
			 alpha=0.25, sync.colors=TRUE, cf=g$contributing_fraction, cf.interval=10, 
			 par.settings=list(superpose.line=list(col=c('RoyalBlue', 'Red4'), lwd=2, lty=c(1,2))),
			 prepanel=prepanel.depth_function,
			 layout=c(5,1), strip=strip.custom(bg=grey(0.8)),
			 scales=list(x=list(tick.number=4, alternating=3, relation='free')),
			 auto.key=list(columns=2, lines=TRUE, points=FALSE)
			 )
```

<img src="figure/slab-2-1.png" title="plot of chunk slab-2" alt="plot of chunk slab-2" style="display: block; margin: auto;" />


### Pair-Wise Dissimilarity
Calculation of between-profile dissimilarity is performed using the `profile_compare()` function. Dissimilarity values depend on attributes selection (e.g. clay, CEC, pH , etc.), optional depth-weighting parameter (`k`), and a maximum depth of evaluation (`max_d`). This topic is further documented in a [related tutorial](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/aqp-profile-dissimilarity.html?root=aqp), the function manual page (type `help(profile_compare)` in the **R** console), and [this paper](http://dx.doi.org/10.1016/j.cageo.2012.10.020).

```r
# eval dissimilarity:
# using Ex-Ca:Mg and CEC at pH 7
# with no depth-weighting (k=0)
# to a maximum depth of 40 cm
d <- profile_compare(sp4, vars=c('ex_Ca_to_Mg', 'CEC_7'), k=0, max_d=40)

# check distance matrix:
round(d, 1)
```

```
## Dissimilarities :
##                colusa glenn kings mariposa mendocino napa san benito shasta shasta-trinity
## glenn            13.5                                                                     
## kings            16.0  12.7                                                               
## mariposa          8.4  11.3  16.5                                                         
## mendocino        11.5   8.0  16.4     15.0                                                
## napa             30.4  24.1  29.4     29.2      21.6                                      
## san benito       25.7  20.6  26.3     28.2      15.8 18.0                                 
## shasta           17.2  13.3   8.7     17.6      17.1 33.7       22.2                      
## shasta-trinity    6.4  16.6  22.3      9.6      16.5 29.8       27.2   23.3               
## tehama           28.7  22.9  27.9     27.3      20.0  8.8       15.1   31.4           27.9
## 
## Metric :  mixed ;  Types = I, I 
## Number of objects : 10
```

```r
# vizualize dissimilarity matrix via hierarchical clustering
h <- hclust(d)
plot(as.dendrogram(h))
```

<img src="figure/dissimilarity-1.png" title="plot of chunk dissimilarity" alt="plot of chunk dissimilarity" style="display: block; margin: auto;" />



### Object Metadata

```r
# print metadata:
metadata(sp4)

# alter the depth unit metadata attribute
depth_units(sp4) <- 'inches' # units are really 'cm'

# more generic interface for adjusting metadata
md <- metadata(sp4) # save original metadata

# add columns
md$describer <- 'DGM'
md$date <- as.Date('2009-01-01')
md$citation <- 'McGahan, D.G., Southard, R.J, Claassen, V.P. 2009. Plant-Available Calcium Varies Widely in Soils on Serpentinite Landscapes. Soil Sci. Soc. Am. J. 73: 2087-2095.'

# re-assign user defined metadata to original object
metadata(sp4) <- md

# check:
metadata(sp4)

# fix depth units, set back to 'cm'
depth_units(sp4) <- 'cm'
```


### Object Structure and Coercion

```r
# check our work by viewing the internal structure
str(sp4)

# deconstruct SoilProfileCollection into a data.frame, with horizon+site data
as(sp4, 'data.frame')

# extraction of site + spatial data as SpatialPointsDataFrame
as(sp4, 'SpatialPointsDataFrame')
```


----------------------------
This document is based on `aqp` version 1.8-6.
