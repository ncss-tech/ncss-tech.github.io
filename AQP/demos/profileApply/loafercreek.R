## ----setup, echo=F, warning=F, message=F---------------------------------
 library(knitr, quietly = TRUE)
opts_chunk$set(message = FALSE,
               warning = FALSE,
               background = '#F7F7F7',
               dpi = 100,
               fig.align = 'center',
               dev = 'png',
               dev.args = list(pointsize = 10),
               cache = FALSE,
               tidy = FALSE)

options(width=100)

## ---- echo=FALSE, fig.width=10, fig.height=6-----------------------------
library(soilDB)
x <- fetchOSD(c("Argonaut", "Auburn", "Bonanza",
                  "Exchequer", "Gardellones", "Gopheridge",
                  "Jasperpeak", "Loafercreek", "Motherlode",
                  "Sobrante"))

# terrible hack i am so evil >:)
x@horizons[which(horizons(x)$hzname == "R"),'bottom'] <- 200
x@horizons[which(horizons(x)$hzname == "R"),'soil_color'] <- NA

sharpshootR::SoilTaxonomyDendrogram(x, 
                                    y.offset = 0.35, 
                                    cex.taxon.labels = 1.2,
                                    cex.names=0.8,
                                    cex.id=1.1,
                                    max.depth=200)

## ----eval=FALSE----------------------------------------------------------
## # install remotes if needed
## # install.packages("remotes")
## remotes::install_github('ncss-tech/aqp', dependencies = FALSE, build = FALSE)
## remotes::install_github('ncss-tech/soilDB', dependencies = FALSE, build = FALSE)
## remotes::install_github('ncss-tech/sharpshootR', dependencies = FALSE, build = FALSE)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(soilDB)

## ---- message=FALSE, warning=FALSE---------------------------------------
data("loafercreek")

## ---- eval=FALSE---------------------------------------------------------
## # access the clay attribute from the horizons data frame
## horizons(spc)$clay
## 
## #add new site data by LEFT JOIN on UNIQUE site ID (assumed to be present in both spc and new.site.data)
## site(spc) <- new.site.data

## ------------------------------------------------------------------------
my.sub.set <- loafercreek[3:6, ]

## ------------------------------------------------------------------------
# number of rows in my.sub.set@site (# of sites or profiles)
nrow(site(my.sub.set))

## ------------------------------------------------------------------------
# adjust margins, units are inches, bottom, left, top, right; adjust as needed
par(mar=c(0, 0, 0, 1))

# make a SoilProfileCollection plot
plotSPC(my.sub.set, label = 'pedon_id', 
        id.style = "side", cex.names = 0.75,
        x.idx.offset = 0.1)

## ------------------------------------------------------------------------
length(loafercreek)

## ------------------------------------------------------------------------
nrow(loafercreek)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(soilDB)
series.names <- c("Argonaut", "Auburn", "Bonanza",
                  "Exchequer", "Gardellones", "Gopheridge", 
                  "Jasperpeak", "Loafercreek", "Motherlode", 
                  "Sobrante")

osds <- fetchOSD(soils = series.names, extended = TRUE)

# adjust margins, units are inches, bottom, left, top, right; modify as needed
par(mar=c(0, 0, 0, 1))

plotSPC(osds$SPC, max.depth=200)

## ---- message=FALSE, echo=FALSE------------------------------------------
k <- fetchKSSL("loafercreek")

hz.match <- '' # match all horizons

## ---- eval=FALSE---------------------------------------------------------
## library(soilDB)
## # use the `series` argument to specify taxonname. or use `mlra` or `bbox`
## # ?fetchKSSL for details
## k <- fetchKSSL(series = "loafercreek")
## 
## # count the number of rows (records) in the `loafercreek@site` data.frame
## n.pedons <- length(k)
## 
## #calculate some basic univariate summary statistics on _all_ horizons in the SPC
## median(k$fe_dith, na.rm = TRUE)
## min(k$fe_dith, na.rm = TRUE)
## max(k$fe_dith, na.rm = TRUE)
## 
## # here you would inspect the data further ...
## # other attributes or by horizon designation, perhaps?
## 

## ------------------------------------------------------------------------
# set spatial coordinates to create a Spatial object
coordinates(loafercreek) <- ~ x_std + y_std

## ------------------------------------------------------------------------
slot(loafercreek, 'sp')

## ------------------------------------------------------------------------
# when you set the proj4string, be sure it matches the formula 
# and the system/format of the data you sent to coordinates() above
proj4string(loafercreek) <- '+proj=longlat +datum=WGS84'

## ------------------------------------------------------------------------
slot(loafercreek, 'sp')

## ------------------------------------------------------------------------
# loafercreek is a SoilProfileCollection
class(loafercreek)

# coerce the loafercreek object to SpatialPointsDataFrame
loafercreek.spdf <- as(loafercreek, 'SpatialPointsDataFrame')

# loafercreek.spdf is a SpatialPointsDataFrame WITH JUST THE SITE DATA 
# the S4 method defined in the SPC does this by design
class(loafercreek.spdf)
nrow(loafercreek.spdf)

## ------------------------------------------------------------------------
# plot pedon locations as points
plot(loafercreek.spdf, sub = "Loafercreek Pedon Locations", pch = 19, cex = 0.5)

# "add" county level maps from maps package
maps::map(database = 'county', regions = 'CA', add = TRUE)

# add a neatline
box()

## ------------------------------------------------------------------------
## make 10 site ids `id` each witin site-level attributes `di` and `mlra`
new.site.data <- data.frame(id = 1:10, di = 10 - 0:9, mlra = "18") 
head(new.site.data)

## or read your data from file
# your.site.data <- read.csv("your.site.data.csv")

## make 10 random horizon datasets, with site id, top and bottom 
## horizon designation and 5 random horizon level attributes (p1, p2...)
your.spc <- do.call('rbind', lapply(new.site.data[['id']], random_profile)) 
head(your.spc)

## or read your data from file. 
# your.spc <- read.csv("your.horizon.data.csv")

#promote horizon data frame (in this case 10 random profiles) to SPC
depths(your.spc) <- id ~ top + bottom

# merge site data into the site slot of the SPC based on common site id `id`
site(your.spc) <- new.site.data

# merge horizon data into the horizon slot of the SPC based on common site id `hzHZ`
# makes some new fake horizon data, take the autogenerated unique horizon ID from the SPC 
new.horizon.data <- data.frame(hzID=your.spc$hzID, newvalue=2)
horizons(your.spc) <- merge(horizons(your.spc), new.horizon.data, by="hzID", all.x = TRUE)

# check to see if the horizon data merge worked
head(your.spc$newvalue)

# attribute names storing critical SPC information
idname(your.spc)
hzidname(your.spc)
horizonDepths(your.spc)

head(profile_id(your.spc)) #unique site/profile ids
head(hzID(your.spc)) #unique horizon ids (assigned at time of SPC creation)

## ---- eval=FALSE---------------------------------------------------------
## profileApply(object, FUN, simplify=TRUE, ...)

## ------------------------------------------------------------------------
depth.to.contact <- profileApply(loafercreek, estimateSoilDepth)

## ------------------------------------------------------------------------
#look at a density (continuous frequency) plot; depth on x axis
plot(density(depth.to.contact, na.rm = TRUE))

## ------------------------------------------------------------------------
quantile(depth.to.contact, 
         probs = c(0,0.01,0.05,0.25,0.5,0.75,0.95,0.99,1), 
         na.rm = TRUE)

## ------------------------------------------------------------------------
bad.peiid <- c("542129") 

#SPC with just the "bad" pedon (this one isn't that bad)
deep.one <- loafercreek[site(loafercreek)$peiid %in% bad.peiid, ]
length(deep.one)

#the inverse, loafercreek without the "bad" one
loafernew <-loafercreek[!(site(loafercreek)$peiid %in% bad.peiid), ]
length(loafernew)

## ------------------------------------------------------------------------
estimateSoilDepth(deep.one)

estimateSoilDepth(deep.one, no.contact.depth = 100, no.contact.assigned = NA)

## ------------------------------------------------------------------------
#plot difference "stored v.s. calculated"
plot(loafercreek$bedrckdepth ~ depth.to.contact, xlim = c(0,200), ylim = c(0,200))

#add a 1:1 line with intercept 0 and slope 1
abline(0, 1)

## ------------------------------------------------------------------------
loafercreek$bedrckdepth <- depth.to.contact

## ------------------------------------------------------------------------
# create new variable with clay proportion, calculated from %
loafercreek$new.hz.level.variable <- loafercreek$clay / 100 

## ------------------------------------------------------------------------
your.function.name <- function(...) { return("output") } 

## ------------------------------------------------------------------------
# your input `tempF` is a numeric vector in degrees fahrenheit
fahrenheit_to_celsius <- function(temp_F) {
  
  # perform the conversion from input to output
  temp_C <- ((temp_F - 32) * (5 / 9))
  
  # return output, a numeric vector in degrees C
  return(temp_C)
}

## ------------------------------------------------------------------------
# a numeric vector in degrees fahrenheit; e.g. Soil Temperature Regime breaks
temp.regime.F <- c(46.4, 59, 71.6)

# call the function we defined, and store the output
temp.regime.C <- fahrenheit_to_celsius(temp.regime.F)

temp.regime.C

## ------------------------------------------------------------------------
profileMaxClay <- function(p) {
    # access the numeric vector `clay` from the horizons data.frame
    d <- horizons(p)$clay
    
    # if all `d` are NA and na.rm = TRUE, max() makes a zero-length variable
    if(all(is.na(d))) 
      return(NA) # so if a pedon `p` has no clay data... return NA pre-emptively
    
    # calculate the maximum value in `d`, omitting NA
    return(max(d, na.rm = TRUE))
}

## ------------------------------------------------------------------------
# calculate the max clay content in all profiles
loafercreek$maxclay <- profileApply(loafercreek, profileMaxClay)

## ------------------------------------------------------------------------
# look at the density plot (estimated probability density function) of maximum clay contents
plot(density(loafercreek$maxclay, na.rm = TRUE, kernel = "rectangular"))

# calculate quantiles
quantile(loafercreek$maxclay, probs = c(0,0.01,0.05,0.25,0.5,0.75,0.95,0.99,1), na.rm = TRUE)

## ------------------------------------------------------------------------
profileMaxClayAttr <- function(p, attr = "clay") {
    # maybe you could calculate something more interesting
    # than the maximum in your version of this function?
    h <- horizons(p)
    
    # pre-emptively return NA if all clay are NA
    if(all(is.na(h$clay))) 
      return(NA)
    
    # here we calculate the horizon indices that equal THE PROFILE max clay
    d.max.idx <- which(h$clay == max(h$clay, na.rm = TRUE))
    
    # there may be multiple horizons with the max clay content... 
    # we just care about the first one (closest to the surface)
    # flattening a (possible) many:one relationship
    d.max.idx.first <- d.max.idx[1]
    
    # return an arbitrary attribute `name` 
    # default returns clay from horizon of max clay
    # you set attr to return attr from horizon of max clay
    return(h[d.max.idx.first, attr])
}

## ------------------------------------------------------------------------
# overwrite with identical values as obtained with profileMaxClay()
# from new generalized function! using default attr='clay'
loafercreek$maxclay <- profileApply(loafercreek, profileMaxClayAttr)

# calculate the DEPTH TO first horizon with max clay content in all profiles
# note we change default attr to return hzdept instead of clay
loafercreek$maxclaydepth <- profileApply(loafercreek, 
                                         profileMaxClayAttr, attr="hzdept")

## ------------------------------------------------------------------------
# look at the density plot (estimated probability density function)
# of minimum depth to maximum clay content
plot(density(loafercreek$maxclaydepth, na.rm = TRUE))

# calculate quantiles
quantile(loafercreek$maxclaydepth, 
         probs = c(0,0.01,0.05,0.25,0.5,0.75,0.95,0.99,1), 
         na.rm = TRUE)

## ------------------------------------------------------------------------
# create a named numeric vector as a "lookup table"
hue.lookup.table <- seq(5, 22.5, 2.5)
names(hue.lookup.table) <- c('5R','7.5R','10R','2.5YR',
                             '5YR','7.5YR','10YR','2.5Y')

## ----echo=F--------------------------------------------------------------
df.lut <- data.frame(names(hue.lookup.table), hue.lookup.table)
names(df.lut) <- c("Dry Hue (H)", "H*")
knitr::kable(df.lut, row.names = FALSE, digits=1,
             caption="Lookup table: Equivalent H and H* Values (after Hurst, 1977)")

## ------------------------------------------------------------------------
# determine H* using the lookup table
hstar <- hue.lookup.table[loafercreek$d_hue]

# calculate Hurst (1977) "Redness Index" H*(L/C)
loafercreek$hri <- hstar * loafercreek$d_value / loafercreek$d_chroma

## ------------------------------------------------------------------------
loafercreek$horizon.is.red <- loafercreek$hri <= 20

summary(loafercreek$horizon.is.red)

## ------------------------------------------------------------------------
loafercreek$depth.to.red <- profileApply(loafercreek, function(p) {
  # access horizon slot of a single profile `p` from loafercreek
  h <- horizons(p)
  
  # calculate indices where horizon.is.red == TRUE, take the first 
  shallowest.match.idx <- which(h$horizon.is.red)[1]
  
  # return top depth of first horizon matching 
  return(h[shallowest.match.idx, 'hzdept'])
})

## ------------------------------------------------------------------------
# calculate top depth of deepest horizon that met our redness criteria
density.cutoff <- max(loafercreek$depth.to.red, na.rm=T)+1
density.cutoff

## ------------------------------------------------------------------------
plot(density(loafercreek$depth.to.red, 
             from = 0, to = density.cutoff,  
             kernel="rectangular", na.rm=T))

## ------------------------------------------------------------------------
# subset of profiles with depth.to.red <= 20
sub1 <- subsetProfiles(loafercreek, s = 'depth.to.red <= 20')

## ------------------------------------------------------------------------
# subset of profiles with any horizon.is.red == TRUE
sub1.alternate <- subsetProfiles(loafercreek, h = 'horizon.is.red == TRUE')

## ------------------------------------------------------------------------
# has red color within 20 cm
length(sub1)

# has red color at any depth
length(sub1.alternate)

## ------------------------------------------------------------------------
# adjust margins, units are inches, bottom, left, top, right
# be sure to leave room for legend
par(mar=c(0, 0, 3, 2))

# convert legend variable to factor for coloring
sub1$horizon.is.red <- factor(sub1$horizon.is.red)

# calculate the plotting order
plotting.order <- order(sub1$depth.to.red)

# make a plot, used some exaggerated red and brown colors to display
# horizon class membership for our threshold
plotSPC(sub1, max.depth = 100, print.id = FALSE, 
        plot.order = plotting.order,
        axis.line.offset = -0.5, name = '',
        col.palette = parseMunsell(c('10YR 5/3', '5YR 3/8')),
        width = 0.4, color = 'horizon.is.red')

# we sorted the plot on "depth to red" to add a line to guide our eye
# plot a white line, with black dots along it
lines(1:length(sub1), sub1$depth.to.red[plotting.order], col = "white", lwd = 2)
lines(1:length(sub1), sub1$depth.to.red[plotting.order], lwd = 3, lty = 3)

# add an axis with depth to redness (corresponds to dotted line)
axis(side = 1, at = 1:length(sub1), cex.axis = 0.5,
     labels = sub1$depth.to.red[plotting.order])
mtext(text = 'Depth to Red (cm)', side = 1, line = 2.5)

## ------------------------------------------------------------------------
# create a logical vector of `peiids` in `loafercreek` that are IN sub1, 
# then invert it with NOT (!)
not.in.sub1 <- !(site(loafercreek)$peiid %in% site(sub1)$peiid)

# square bracket SPC subsetting
sub2 <- loafercreek[not.in.sub1,]

# how many in the "remainder" set?
length(sub2)

## ----fig.width=10, fig.height=4.5----------------------------------------
par(mar=c(0, 0, 3, 2))

# convert legend variable to factor for coloring
sub2$horizon.is.red <- factor(sub2$horizon.is.red)

# depth cut off at 100cm, hide IDs and make it clear which have data
plotSPC(sub2, max.depth=100, print.id = FALSE, 
        axis.line.offset=-0.5, name = '', 
        width = 0.4, color = 'horizon.is.red',
        col.palette = parseMunsell(c('10YR 5/3', '5YR 3/8')), 
        col.legend.cex=1.25, col.label='Horizon is red?')

## ------------------------------------------------------------------------
# plot estimate of the probability density function for redness index
plot(density(sub1$hri, na.rm=T))

# add a vertical line at HRI=20 (where we put our "break" in the data)
abline(v=20, lty=2, lwd=2, col="RED")

## ------------------------------------------------------------------------
# modify the  "genhz" (NASIS Comp. Layer ID) to combine BA and A
# we do this because there is only one BA horizon with dry color data
loafercreek$genhz[loafercreek$genhz == "BA"] <- "A"

# make a list of data frames split by "genhz" (NASIS Comp. Layer ID)
genhz.list <- split(horizons(sub1), f=sub1$genhz)

# calculate some quantiles of Hurst Redness Index for each genhz
# lapply applies an anoymous function to each data frame in genhz.list
qtiles.by.genhz <- do.call('rbind', lapply(genhz.list, function(d) {
  n.obs <- sum(!is.na(d$hri))
  names(n.obs) <- "n.obs"
  return(c(quantile(d$hri, 
                    probs=c(0,0.05,0.25,0.5,0.75,0.95,1), 
                    na.rm=TRUE), n.obs))
}))

#remove NA rows
qtiles.by.genhz  <- qtiles.by.genhz[complete.cases(qtiles.by.genhz),]

## ----echo=FALSE----------------------------------------------------------
# reorder soil horizons genetically and make a table
knitr::kable(qtiles.by.genhz[c("A", "Bt1", "Bt2", "Bt3", "BCt"),], digits = 0,
             caption="Quantiles of Hurst Redness Index - grouped by NASIS Component Layer ID")

## ------------------------------------------------------------------------
loafercreek$red.shallow <- loafercreek$depth.to.red <= 20

## ---- eval=FALSE---------------------------------------------------------
## sum(is.na(loafercreek$red.shallow))
## 
## length(loafercreek)

## ------------------------------------------------------------------------
all.na  <- profileApply(loafercreek, 
                        function(p) {
                          return(all(is.na(p$horizon.is.red))) 
                        })

## ------------------------------------------------------------------------
# change value of profiles that have NA red.shallow, without all NA color
loafercreek$red.shallow[is.na(loafercreek$red.shallow) & !all.na] <- 2

## ------------------------------------------------------------------------
summary(factor(loafercreek$red.shallow))

## ------------------------------------------------------------------------
# adjust margins, units are inches, bottom, left, top, right
par(mar=c(0, 0, 0, 2))

plotSPC(loafercreek[which(is.na(loafercreek$red.shallow)),], 
     color = 'd_hue', max.depth=100, 
     print.id = FALSE, name = '', axis.line.offset=-0.5)

## ----fig.width=12, fig.height=5.5----------------------------------------
# adjust margins, units are inches, bottom, left, top, right
# be sure to leave room for legend
par(mar=c(0, 0, 3, 2))
loafercreek$horizon.is.red <- factor(loafercreek$horizon.is.red)

# NOTE: GPPs DO allow NA in the grouping variable; plots as '<missing>'
# ...but a warning is generated
groupedProfilePlot(loafercreek,
                   group.name.cex = 1.5, groups = 'red.shallow', 
                   color = 'horizon.is.red', max.depth = 100, 
                   print.id = FALSE, name = '', width = 0.4, divide.hz = FALSE,
                   col.palette = parseMunsell(c('10YR 5/3', '5YR 3/8')), 
                   col.legend.cex = 1.25, col.label = 'Horizon is red?')

## ------------------------------------------------------------------------
# make some labels
#loafercreek$red.shallow :   1      0         2        NA
redness.levels <-        c("RED","DEEPRED", "NOTRED","NODATA")

# set the value of the NAs to 3
loafercreek$red.shallow[is.na(loafercreek$red.shallow)] <- 3

# create a site-level grouping factor
loafercreek$redness.class <- factor(loafercreek$red.shallow, 
                                    levels = c(1,0,2,3), 
                                    labels = redness.levels)

## ------------------------------------------------------------------------
summary(loafercreek$redness.class)

## ----fig.width=12, fig.height=5.5----------------------------------------
# adjust margins, units are inches, bottom, left, top, right
# be sure to leave room for legend
par(mar=c(0, 0, 3, 2))

groupedProfilePlot(loafercreek, max.depth=100, 
                   group.name.cex = 1.5, groups = 'redness.class', 
                   color = 'horizon.is.red', print.id = FALSE, 
                   name = '', axis.line.offset=-0.5, 
                   width = 0.4, divide.hz = FALSE,
                   col.palette = parseMunsell(c('10YR 5/3', '5YR 3/8')),
                   col.legend.cex=1.25, col.label='Horizon is red?')

## ----echo=FALSE----------------------------------------------------------
who.idx <- which(loafercreek$redness.class == "RED" | 
                   loafercreek$redness.class == "NOTRED") 
plot(density(loafercreek$maxclay[who.idx], na.rm = TRUE, ), 
     type="n", xlim = c(0, 60), ylim = c(0, 0.1), 
     main = "Probability density of profile maximum clay content by \"redness\" class", 
     sub = "Maximum Clay %; RED v.s NOTRED")

sub.idx <- c(2,4)
# set up plotting arguments
line.labelz <- c("ALL", levels(loafercreek$redness.class))
line.colorz <- c("BLACK","DARKRED","RED","BLUE","PURPLE")
plot.lty <- c(3,1,1,1,1)

# CREATE DATA FRAME (NO FACTORS TO PRESERVE ORDERING)
plot.params <- data.frame(labels=line.labelz, 
                          line.color=line.colorz, 
                          lty=plot.lty, stringsAsFactors=FALSE)

plot.params <- plot.params[sub.idx,]

# make a base plot with base R apply() :)
res <- apply(plot.params, MARGIN=1, FUN=function(i) {
  idx <- loafercreek$redness.class %in% i[['labels']]
  
  if(all(!idx)) # handle 'ALL' which is not a factor level; it is al levels
    idx <- !idx
  
  lines(density(loafercreek$maxclay[idx], 
                na.rm = TRUE, from=0, to=60, kernel="rectangular"), 
                lty=as.numeric(i[['lty']]), col=i[['line.color']], lwd = 2)
})

legend(x = 45, y = 0.1025, cex = 0.9,
       legend = plot.params$labels, col = plot.params$line.color, 
       lwd = 2, lty = plot.params$lty)

## ------------------------------------------------------------------------
# compare groups versus full set. Empty plot.
plot(density(loafercreek$maxclay, na.rm = TRUE, ), 
     type="n", xlim = c(0, 60), ylim = c(0, 0.1), 
     main = "Probability density of profile\n maxiumum clay content by \"redness\" class", 
     sub = "Maximum Clay %; Subsets versus full `loafercreek` group")

# set up plotting arguments
line.labelz <- c("ALL", levels(loafercreek$redness.class))
line.colorz <- c("BLACK","DARKRED","RED","BLUE","PURPLE")
plot.lty <- c(3,1,1,1,1)

# CREATE DATA FRAME (NO FACTORS TO PRESERVE ORDERING)
plot.params <- data.frame(labels=line.labelz, 
                          line.color=line.colorz, 
                          lty=plot.lty, stringsAsFactors=FALSE)

# make a base plot with base R apply() :)
res <- apply(plot.params, MARGIN=1, FUN=function(c) {
  idx <- loafercreek$redness.class %in% c[['labels']]
  
  if(all(!idx)) # handle 'ALL' which is not a factor level... it is all factor levels
    idx <- !idx
  
  lines(density(loafercreek$maxclay[idx], na.rm = TRUE, from=0, to=60,
                kernel="rectangular"), lty=as.numeric(c[['lty']]), 
                col=c[['line.color']], lwd = 2)
})

legend(x = 45, y = 0.1025, cex = 0.9,
       legend = plot.params$labels, col = plot.params$line.color, 
       lwd = 2, lty = plot.params$lty)

## ------------------------------------------------------------------------
table(toupper(loafercreek$taxonname))

table(toupper(loafercreek$taxonkind))

## ------------------------------------------------------------------------
# modify the  "genhz" (NASIS Comp. Layer ID) to combine BA and A
# we do this because there is only one BA horizon with dry color data
loafercreek$genhz[loafercreek$genhz == "BA"] <- "A"
loafercreek$genhz[is.na(loafercreek$genhz)] <- "not-assigned"

# create a horizon level redness grouping factor from site
# TODO: make an easier method for doing this in the SPC object
loafercreek$redness <- merge(horizons(loafercreek), 
                            site(loafercreek)[,c('peiid','redness.class')], 
                            all.x = TRUE)$redness.class

# make a list of data frames split by "redness" and "genhz" (NASIS Comp. Layer ID)
red.genhz.list <- split(horizons(loafercreek), 
                        f = list(loafercreek$redness, loafercreek$genhz))

## ------------------------------------------------------------------------
# calculate some quantiles of clay content for each redness.class*genhz
qtiles.redgenhz <- do.call('rbind', lapply(red.genhz.list, function(d) {
  
  # add number of obervations (named numeric vector) to output
  n.obs <- sum(!is.na(d$clay))
  names(n.obs) <- "n.obs"
  
  # calculate quantiles, concatenate n.obs & return
  return(data.frame(q = t(quantile(d$clay,
                                  probs = c(0,0.05,0.25,0.5,0.75,0.95,1), 
                                  na.rm = TRUE)), 
                    n.obs, 
                    redness = d$redness[1], 
                    genhz = d$genhz[1]))
}))

## ----results='asis'------------------------------------------------------
# print a list of data frames, split by redness class
library(knitr)

red.clayl <- split(qtiles.redgenhz, f = qtiles.redgenhz$redness)

res <- lapply(red.clayl, function(d) {
  rownames(d) <- d$genhz
  
  print(kable(d[c("A", "Bt1", "Bt2", "Bt3", "BCt", "Cr", "not-assigned"), ], 
              caption = "Selected Quantiles of Clay Content - 
                          grouped by NASIS Component Layer ID & Redness Group"))
})

## ------------------------------------------------------------------------
summary(loafercreek$redness.class)

## ---- eval=FALSE---------------------------------------------------------
## profileApply(object, FUN, simplify=TRUE, ...)

## ------------------------------------------------------------------------
# create an SPC with just one pedon (try different ones)
just.one <- loafercreek[1]

the.result <- estimateSoilDepth(just.one)
the.result

applied.result <- profileApply(just.one, estimateSoilDepth)
applied.result

## ------------------------------------------------------------------------
(the.result == applied.result)

## ------------------------------------------------------------------------
names(applied.result)

str(applied.result)

## ------------------------------------------------------------------------
# logical vector length of left hand side (1)
# "is left hand side character vector IN right hand side character vector?"
names(applied.result) %in% loafercreek$peiid

# get site index (we took the first one at the beginning)
# NB a vector of length equal to loafercreek number of site/profiles results
which(site(loafercreek)$peiid %in% names(applied.result))

# If we don't specify site(), we get `peiid` from @horizons. 
# So we get several horizon indexes, all matching peiid and belonging to the first site
# NB This may or may NOT be what you expect/want!
which(loafercreek$peiid %in% names(applied.result))

## ------------------------------------------------------------------------
numeric.vector <- profileApply(loafercreek, estimateSoilDepth)

head(numeric.vector, 3) # named numeric vector, names are peiid
class(numeric.vector)   # numeric
typeof(numeric.vector)  # double precision numeric

## ------------------------------------------------------------------------
a.list <- profileApply(loafercreek, estimateSoilDepth, simplify = FALSE)

head(a.list, 3)     # a named list, names are peiid

class(a.list)       # list can contain a mix of any data type

typeof(a.list[[1]]) # the first element of this list is numeric (integer)

str(unlist(a.list)) # create a named numeric vector from list (since all are numeric)

str(as.numeric(a.list)) # create an UNnamed numeric vector from list

## ------------------------------------------------------------------------
depth.weighted.average <- function(spc, tdepth, bdepth, attr, ...) {
  #expand `attr` in formula
  custom.formula <- formula(paste0(tdepth,":",bdepth," ~ ", paste0(attr, collapse=" + ")))
  
  # calculate a depth-weighted average using aqp::slice()
  return(mean(slice(spc, custom.formula, just.the.data=TRUE)[[attr]],
              na.rm = TRUE))
}

## ------------------------------------------------------------------------
#dwt.mean.wrap()
# x - a SoilProfileCollection; usually single-profile from profileApply()
# na.threshold - specifies the proportion (by number of records) of NA for warning 
# 
# Pedons with all NA for `attr` return NA (with a message containing peiid)
# Generate a warning if NA is over a threshold (with message containing peiid)

dwt.mean.wrap <- function(x, na.threshold, attr, ...) {
  # create local variable to limit SPC accessor calls
  my.attr <- horizons(x)[[attr]]
  
  # if all records are NA don't try to process the data, flag & return NA
  if(all(is.na(my.attr))) {
    print(paste0("All horizons are `NA` \'", attr, "\' -- PEIID:", site(x)$peiid))
    return(NA)
  }
  
  # calculate proportion of NA records and check against threshold
  prop.na <- sum(is.na(my.attr)) / length(my.attr)
  if(prop.na > na.threshold) {
    print(paste0("Greater than ", round(na.threshold * 100),
                 "% of horizons are `NA` \'", attr, "\' - PEIID:", site(x)$peiid))
    #we will not omit it, just flag it so we have the ID
  }
  
  # we used x and attr to max the above output, but depth.weighted.average()
  # also gets tdepth and bdepth from call to profileApply(), passed on with  `...`
  return(depth.weighted.average(x, attr, ...))
  # then we return the result
}

## ------------------------------------------------------------------------

wrapper.test <- profileApply(loafercreek, 
                             dwt.mean.wrap, 
                             na.threshold = 0.40,
                             tdepth = 25, bdepth = 75,
                             attr = 'clay')

plot(density(wrapper.test, na.rm=TRUE), 
     main = "25-75cm depth-weighted average clay content\nall `loafercreek` pedons")

sum(is.na(wrapper.test))

## ------------------------------------------------------------------------
#diagnostic horizon validation (dark surface, argillic)

## ------------------------------------------------------------------------
#particle size control section depths and weighted-average calculation

## ------------------------------------------------------------------------
#horizon-level validations and preparing SPC objects for analysis

## ------------------------------------------------------------------------
#do spatial example.. can we predict where the red/clayey ones are?

## ------------------------------------------------------------------------
#bring in the shallow and skeletal and deep data from the Loafercreek mapunits 

