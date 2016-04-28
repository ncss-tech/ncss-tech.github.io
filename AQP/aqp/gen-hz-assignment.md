---
output:
  html_vignette:
    mathjax: null
    jquery: null
    smart: no
---




Assigning Generalized Horizon Labels
===============================
D.E. Beaudette, J.M. Skovlin (edited for clarity by Jenny Sutherland)
<br>
2015-04-22
<br>
This document is based on `aqp` version 1.8-6 and `soilDB` version 1.5-3.


# Introduction

## An Example
Consider this situation: you have a collection of pedons that have been correlated to a named soil series (or component) and would like to *objectively* compute a range in characteristics ("low-rv-high" values) and horizon depths. As with most collections of pedon data, there may be considerable variation in description style and horizons used, horizon depths, and number of horizons described:

![alt text](genhz-sketch.png)

In this scenario, there are several obvious &ldquo;micro-correlation&rdquo; decisions that need to be made before horizons can be grouped for aggregation. For example, what horizonation prototype scheme (e.g., A-Bt1-Bt2-Bt3-Cr-R) best conveys the concept of this soil series or soil component? Does it make sense to group [Bt3, Bt4, BCt, CBt] horizons for aggregation? Or what about grouping [Bt3, 2Bt3] horizons? Do [BA, AB] horizons occur frequently enough to be included in the horizonation prototype?

Based on your knowledge of the area, pedon 2 might be a good "typical" pedon to use in developing a horizonation prototype. After careful review of the data and consultation with your crew, a new set of labels are assigned to each horizon (<span style="color: red;">red labels in figure above</span>) that define groups over which soil properties will be aggregated. These new labels define *functionally similar* groups that may span multiple genetic horizons.


## Formalized Approach 
This document describes an aggregation strategy based on the use of "generalized horizon labels" (GHL) through a combination of narrative, **R** code, and figures. You can follow along in an **R** session by copying code from the grey boxes and pasting it into the **R** console.

Here is a basic outline of the process:

  1. Select a set of GHL that best represents a group of pedons to be aggregated. This could be based on series descriptions, expert knowledge, or even inspection of the most frequently occurring horizon designations.
  
  2. Assign GHL to each horizon using whatever information  is available for grouping horizons. This micro-correlation of horizon designations will likely require slightly different "rules" in each instance. Careful inspection of horizon designation and observed properties is critical.
  
  3. Evaluate GHL assignments and manually refine as needed. 
  
  4. Keep track of final GHL assignments in NASIS or text file. 
  
  5. [Estimate a most likely horizonation](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/estimate-ML-horizonation.html?root=aqp) e.g., top and bottom depths for each generalized horizon label.
  
  6. [Compute range in characteristics](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/gen-hz-application.html?root=aqp), aka [low-rv-high values](http://casoilresource.lawr.ucdavis.edu/wiki/Low-rv-high), for clay, sand, pH, etc. over GHL. (next document in this series)
  



# Setup R Envionment
If you have never used the [aqp](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/aqp-intro.html?root=aqp) or [soildb](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/soilDB/soilDB-Intro.html?root=aqp) packages before, you will likely need to install them. This only needs to be done once. 

```r
# stable version from CRAN + dependencies
install.packages('ape', dep=TRUE) 
install.packages('latticeExtra', dep=TRUE)
install.packages('plyr', dep=TRUE) 
install.packages('aqp', dep=TRUE) 
install.packages('soilDB', dep=TRUE)
# latest versions from R-Forge:
install.packages('aqp', repos="http://R-Forge.R-project.org")
install.packages('soilDB', repos="http://R-Forge.R-project.org")
```

Once you have all of the R packages on which  this document depends, it is a good idea to load them. R packages must be **installed** anytime you change versions of R (e.g., after an upgrade) and **loaded** anytime you want to access functions from within those packages.


```r
library(aqp)
library(soilDB)
library(ape)
library(latticeExtra)
library(plyr)
library(lattice)
library(cluster)
library(MASS)
```

# Sample Data
While the methods outlined in this document can be applied to any collection of pedons, it is convenient to work with a standardized set of data. You can follow along with the analysis by copying code from the following blocks and running it in your **R** session. The sample data used in this document is based on 15 soil profiles that have been correlated to the [Loafercreek](https://soilseries.sc.egov.usda.gov/OSD_Docs/L/LOAFERCREEK.html) soil series from the Sierra Nevada Foothill Region of California. Note that the internal structure of the `loafercreek` data is identical to the structure returned by [`fetchNASIS()` from the soilDB package](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/soilDB/soilDB-Intro.html?root=aqp). All horizon-level values are pulled from the pedon horizon table of the pedons being analyzed.


```r
# load sample data from the soilDB package
data(loafercreek, package = 'soilDB')
# keep only the first 15 pedons
pedons <- loafercreek[1:15, ]
# plot profile sketches
par(mar=c(0,0,0,0))
plot(pedons, name='hzname', print.id=FALSE, cex.names=0.8, axis.line.offset=-4)
```

<img src="figure/load-data-1.png" title="plot of chunk load-data" alt="plot of chunk load-data" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">15 pedons correlated to the Loafercreek soil series.</p><hr>


## Optional: Follow Along with Your Data
The following code block demonstrates how to pull data in using the [`fetchNASIS()` function from the soilDB package](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/soilDB/soilDB-Intro.html?root=aqp).


```r
# first load the desired data set within NASIS into your NASIS selected set
# then load data from the NASIS selected set into R
# note that the `nullFragsAreZero` argument converts NULL rock fragment class data into 0s
pedons <- fetchNASIS(nullFragsAreZero=TRUE)
# optionally subset the data by taxon name - enter your taxon name
pedons <- pedons[grep(pattern='ENTER_YOUR_TAXON_NAME', pedons$taxonname, ignore.case=TRUE), ]
```

# Methods

## Selection of Generalized Horizon Labels

Generalized horizon labels represent an expert-guided selection of designations that were consistently observed in the field and are meaningful in terms of soil morphology and management. The very first step in this process is to tabulate the number of times each horizon designation occurs.


```r
sort(table(pedons$hzname), decreasing=TRUE)
```

```
## 
##    A  Bt1  Bt2  Bt3   Oi    R   Cr  Crt   BA  BCt  Bt4   Bw  CBt 2BCt 2Bt3  2CB  2Cr   2R 
##   15   15   15   10    9    9    8    5    3    3    2    2    2    1    1    1    1    1
```

In this case, [A, Bt1, Bt2, Bt3] horizon designations appear to be a good starting point. However, relying on field experience and expert knowledge of these soils, cross-checking against the [OSD](https://soilseries.sc.egov.usda.gov/OSD_Docs/L/LOAFERCREEK.html), and talking with other soil scientists may be more important than the previous tabulation. 

A quick summary of horizon depth mid-points (e.g., average depth of horizon) can help in organizing the various designations and possibly give some clues as to how they can be grouped. The following plot is called a [box and whisker plot](http://en.wikipedia.org/wiki/Box_plot).



```r
# compute horizon mid-points
pedons$mid <- with(horizons(pedons), (hzdept + hzdepb) / 2)

# sort horizon designation by group-wise median values
hz.designation.by.median.depths <- names(sort(tapply(pedons$mid, pedons$hzname, median)))

# plot the distribution of horizon mid-points by designation
bwplot(mid ~ factor(hzname, levels=hz.designation.by.median.depths), 
       data=horizons(pedons), 
       ylim=c(155, -5), ylab='Horizon Mid-Point Depth (cm)', 
       scales=list(y=list(tick.number=10)), 
       panel=function(...) {
  panel.abline(h=seq(0, 140, by=10), v=1:length(hz.designation.by.median.depths), col=grey(0.8), lty=3)
	panel.bwplot(...)
})
```

<img src="figure/horizonation-mid-point-1.png" title="plot of chunk horizonation-mid-point" alt="plot of chunk horizonation-mid-point" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Range in horizon depth mid-point for original horizon designations.</p><hr>

Next, a similar summary of soil properties (clay content and total rock fragment volume) is presented. The goal is to determine which horizons designations can be grouped and which generalized horizon labels to assign to each group.


```r
# box and wisker plot by clay content
bwplot(clay ~ factor(hzname, levels=hz.designation.by.median.depths), 
       data=horizons(pedons), 
       ylab='Clay Content (%)', 
       scales=list(y=list(tick.number=10)), 
       panel=function(...) {
  panel.abline(h=seq(0, 100, by=5), v=1:length(hz.designation.by.median.depths), col=grey(0.8), lty=3)
  panel.bwplot(...)
})
```

<img src="figure/univariate-eval-clay-1.png" title="plot of chunk univariate-eval-clay" alt="plot of chunk univariate-eval-clay" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Range in clay content for original horizon designations.</p><hr>


```r
# box and wisker plot by total rock fragment volume
bwplot(total_frags_pct ~ factor(hzname, levels=hz.designation.by.median.depths), 
       data=horizons(pedons), 
       ylab='Total Rock Fragment Volume (%)', 
       scales=list(y=list(tick.number=10)), 
       panel=function(...) {
  panel.abline(h=seq(0, 100, by=10), v=1:length(hz.designation.by.median.depths), col=grey(0.8), lty=3)
  panel.bwplot(...)
})
```

<img src="figure/univariate-eval-rf-1.png" title="plot of chunk univariate-eval-rf" alt="plot of chunk univariate-eval-rf" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Range in total rock fragment content for original horizon designations.</p><hr>


Sometimes looking at thematic soil profile sketches can be informative.

```r
# color horizons by clay content
par(mar=c(0,0,3,3))
plot(pedons, name='hzname', print.id=FALSE, cex.names=0.8, axis.line.offset=-4, color='clay')
```

<img src="figure/thematic-profile-sketch-1.png" title="plot of chunk thematic-profile-sketch" alt="plot of chunk thematic-profile-sketch" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Horizon colors are based on clay content.</p><hr>

## Assignment of Generalized Horizon Labels

Once a set of generalized horizon labels have been determined, a corresponding set of [regular expression](http://en.wikipedia.org/wiki/Regular_expression) (REGEX) rules are developed to convert field-described designations into GHL. Pattern matching with REGEX will typically assign useful GHL; however, there will always be cases where manual intervention is required. More on that will follow.

From the above analysis and the [OSD](https://soilseries.sc.egov.usda.gov/OSD_Docs/L/LOAFERCREEK.html), it seems like the following sequence of GHL is appropriate: (A, Bt1, Bt2, Bt3, Cr, R) -- an A horizon, followed by 3 Bt horizons, then Cr, and finally R. For each GHL, a corresponding REGEX rule is needed. For example, `'^A$|Ad|Ap'` will match 'A', 'Ad', and 'Ap'.


```r
# save our GHL
n <- c('A','Bt1','Bt2','Bt3','Cr','R')
# REGEX rules
p <- c('^A$|Ad|Ap',
       'Bt1$',
       '^Bt2$',
       '^Bt3|^Bt4|CBt$|BCt$|2Bt|2CB$|^C$',
       'Cr',
       'R')
```


Apply GHL pattern-matching rules, save to a new column called `genhz`, and cross-tabulate the occurrence of GHL and original designations.

```r
pedons$genhz <- generalize.hz(pedons$hzname, n, p)
# cross-tabulate original horizon designations and GHL
addmargins(table(pedons$genhz, pedons$hzname))
```



|/        | 2BCt| 2Bt3| 2CB| 2Cr| 2R|  A| BA| BCt| Bt1| Bt2| Bt3| Bt4| Bw| CBt| Cr| Crt| Oi|  R| Sum|
|:--------|----:|----:|---:|---:|--:|--:|--:|---:|---:|---:|---:|---:|--:|---:|--:|---:|--:|--:|---:|
|A        |    0|    0|   0|   0|  0| 15|  0|   0|   0|   0|   0|   0|  0|   0|  0|   0|  0|  0|  15|
|Bt1      |    0|    0|   0|   0|  0|  0|  0|   0|  15|   0|   0|   0|  0|   0|  0|   0|  0|  0|  15|
|Bt2      |    0|    0|   0|   0|  0|  0|  0|   0|   0|  15|   0|   0|  0|   0|  0|   0|  0|  0|  15|
|Bt3      |    1|    1|   1|   0|  0|  0|  0|   3|   0|   0|  10|   2|  0|   2|  0|   0|  0|  0|  20|
|Cr       |    0|    0|   0|   1|  0|  0|  0|   0|   0|   0|   0|   0|  0|   0|  8|   5|  0|  0|  14|
|R        |    0|    0|   0|   0|  1|  0|  0|   0|   0|   0|   0|   0|  0|   0|  0|   0|  0|  9|  10|
|not-used |    0|    0|   0|   0|  0|  0|  3|   0|   0|   0|   0|   0|  2|   0|  0|   0|  9|  0|  14|
|Sum      |    1|    1|   1|   1|  1| 15|  3|   3|  15|  15|  10|   2|  2|   2|  8|   5|  9|  9| 103|

In the above cross-tabulation, we can see that a couple of original designations were not matched (`not-used` in the table)  by our REGEX rules, namely, BA, Bw, and Oi horizons. For this example, we are going to  assume that those horizons are not common enough for inclusion in our set of GHL.

## Evaluation of Generalized Horizon Labels

An initial evaluation of our GHL assignment can be accomplished by plotting profile sketches with horizons colored by GHL. The result looks correct, but further investigation may be warranted.

```r
# make a palette of colors, last color is for not-used class
cols <- c(grey(0.33), 'orange', 'orangered', 'chocolate', 'green', 'blue', 'yellow')
# assign a color to each generalized horizon label
hz.names <- levels(pedons$genhz)
pedons$genhz.soil_color <- cols[match(pedons$genhz, hz.names)]
# plot generalized horizons via color and add a legend
par(mar=c(4,0,0,0))
plot(pedons, name='hzname', print.id=FALSE, cex.names=0.8, axis.line.offset=-4, color='genhz.soil_color')
legend('bottomleft', legend=hz.names, pt.bg=c(cols), pch=22, bty='n', cex=2)
```

<img src="figure/plot-ghl-1-1.png" title="plot of chunk plot-ghl-1" alt="plot of chunk plot-ghl-1" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Horizon colors are based on assigned GHL.</p><hr>


A box and whisker plot of the depth-ranges for each of the generalized horizon labels helps to visualize the degree of overlap. [Slicing](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/aqp-intro.html?root=aqp) the collection of profiles along 1-cm intervals generates a more precise figure.

```r
# slice profile collection from 0-150 cm
s <- slice(pedons, 0:150 ~ genhz)
# convert horizon name back to factor, using original levels
s$genhz <- factor(s$genhz, levels = hz.names)
# plot depth-ranges of generalized horizon slices
bwplot(hzdept ~ genhz, data=horizons(s), 
       ylim=c(155, -5), ylab='Generalized Horizon Depth (cm)', 
       scales=list(y=list(tick.number=10)), asp=1, 
       panel=function(...) {
          panel.abline(h=seq(0, 140, by=10), v=1:length(hz.names),col=grey(0.8), lty=3)
	        panel.bwplot(...)
          }
       )
```

<img src="figure/plot-ghl-2-1.png" title="plot of chunk plot-ghl-2" alt="plot of chunk plot-ghl-2" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Range in GHL horizon depth.</p><hr>




### Advanced: Multivariate Soil Property Summary

The evaluation of generalized horizon labels typically requires a review of more information than field-described horizon labels and depth. For example, clay content, sand content, pH, total rock fragment volume, and horizon mid-points are soil properties that can be used to differentiate horizons-- as long as the data are populated. In this document, clay content, total rock fragment volume, and horizon mid-points are used. Multivariate comparisons are commonly based on the concept of [pair-wise dissimilarity](http://hymenoptera.tamu.edu/courses/ento601/pdf/Sokal_1966.pdf) and subsequent reduction of the resulting [distance matrix](http://en.wikipedia.org/wiki/Distance_matrix) into a simpler representation. [Non-metric multidimensional scaling](http://en.wikipedia.org/wiki/Multidimensional_scaling) (nMDS) is one method for reducing the distance matrix into a new set of coordinates in two-dimensional space. A [related](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/aqp-profile-dissimilarity.html?root=aqp) document explores this idea further.


```r
# store the column names of our variables of interest
vars <- c('clay', 'mid', 'total_frags_pct')
# result is a list of several items
hz.eval <- evalGenHZ(pedons, 'genhz', vars)
# extract MDS coords
pedons$mds.1 <- hz.eval$horizons$mds.1
pedons$mds.2 <- hz.eval$horizons$mds.2
# extract silhouette widths and neighbor
pedons$sil.width <- hz.eval$horizons$sil.width
pedons$neighbor <- hz.eval$horizons$neighbor
```


In the following figure, generalized horizon labels are plotted at nMDS coordinates as colored circles, along with original horizon designations and pedon IDs. Ideally, groupings of GHL should form relatively homogeneous clusters in this figure; however, there will always be some overlap at the edges. Heterogeneous regions should be inspected in cases where the assigned GHL may not make sense. In this example, the original designation of "Bt3" for pedon ID 09CKS036 does not seem to fit into the "Bt3" generalized horizon. According to the soil properties used (clay content, total rock fragments, horizon mid-point), the horizon may  better fit the "Bt2" generalized horizon. In addition, it appears that "Bw" and "BA" horizons can be included into the "A" generalized horizon. Decisions on issues such as this can only be made based on field experience or at least some level of expert knowledge about the soils in question. Once a decision has been made, additions to the GLH rules (e.g. adding "Bw" and "BA" to the "A" GHL) and possibly manual adjustment can be used to accomodate the changes.

```r
# convert pedons to a data.frame
pedons.df <- as(pedons, 'data.frame')
# plot generalized horizon labels at MDS coordinates
mdsplot <- xyplot(mds.2 ~ mds.1, groups=genhz, data=pedons.df, 
                  xlab='', ylab='', aspect=1,
                  scales=list(draw=FALSE), 
                  auto.key=list(columns=length(levels(pedons.df$genhz))), 
                  par.settings=list(
                    superpose.symbol=list(pch=16, cex=3, alpha=0.5)
                  )
)

# annotate with original hzname and pedon ID
mdsplot +
  layer(panel.abline(h=0, v=0, col='grey', lty=3)) + 
  layer(panel.text(pedons.df$mds.1, pedons.df$mds.2, pedons.df$hzname, cex=0.85, font=2, pos=3)) +
  layer(panel.text(pedons.df$mds.1, pedons.df$mds.2, pedons.df$pedon_id, cex=0.55, font=1, pos=1))
```

<img src="figure/mds-plot-1.png" title="plot of chunk mds-plot" alt="plot of chunk mds-plot" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Horizon designations and GHL as plotted along nMDS axes.</p><hr>


The [silhouette width](http://www.sciencedirect.com/science/article/pii/0377042787901257) metric is a useful indicator of how consistently group labels partition observations in a dissimilarity matrix. Silhouette widths closer to 1 indicate excellent partitioning, while values closer to 0 indicate less effective partitioning. Silhouette widths less than 0 are often associated with inappropriate label assignment. A sketch of profiles with horizons colored by silhouette width can indicate possible horizons with inappropriate GHL assignment or a pedon that does not fit within the concept of Loafercreek. For example, many of the horizons associated with pedon ID 09CKS036 do not appear to fit within the collection of soil properties associated with generalized horizon labels "Bt1", "Bt2", and "Bt3". Keep in mind that this evaluation is only based on three soil properties-- expert judgement will be required for a final assignment of GHL.

```r
# plot silhouette width metric
par(mar=c(0,0,3,3))
plot(pedons, name='hzname', label='pedon_id', cex.names=0.8, axis.line.offset=-4, color='sil.width')
```

<img src="figure/silhouette-sketch-1.png" title="plot of chunk silhouette-sketch" alt="plot of chunk silhouette-sketch" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Blue and green colors warrant further investigation.</p><hr>


The following table contains information on those horizons with silhouette widths less than 0. The "neighbor" column contains the next closest GHL (in terms of the soil properties used in the pair-wise dissimilarity calculation), which may be more appropriate to use.

```r
# index those horizons with silhouette widths less than 0
check.idx <- which(pedons.df$sil.width < 0)
# sort this index based on min sil.width
check.idx.sorted <- check.idx[order(pedons.df$sil.width[check.idx])]
# list those pedons/horizons that may need some further investigation
pedons.df[check.idx.sorted, c('peiid', 'pedon_id', 'hzname', 'genhz', 'neighbor', 'sil.width', vars)]
```



|    |  peiid|pedon_id |hzname |genhz |neighbor |  sil.width| clay|  mid| total_frags_pct|
|:---|------:|:--------|:------|:-----|:--------|----------:|----:|----:|---------------:|
|46  | 374232|09CKS036 |Bt2    |Bt2   |Bt1      | -0.5085458|   21| 23.5|               5|
|47  | 374232|09CKS036 |Bt3    |Bt3   |Bt2      | -0.4101753|   23| 40.0|              10|
|45  | 374232|09CKS036 |Bt1    |Bt1   |A        | -0.3597560|   16|  9.5|               3|
|87  | 477042|10CKS015 |Bt2    |Bt2   |Bt1      | -0.2492451|   21| 36.5|               2|
|79  | 425293|10MJE005 |A      |A     |Bt1      | -0.2180060|   23|  1.5|               0|
|69  | 414936|10CKS009 |Bt3    |Bt3   |Bt2      | -0.1389465|   20| 53.5|              34|
|75  | 414938|10CKS011 |Bt2    |Bt2   |Bt1      | -0.0644906|   24| 30.5|              13|
|17  | 374201|09CKS006 |Bt3    |Bt3   |Bt2      | -0.0635606|   26| 40.0|              48|
|56  | 414895|10CKS001 |Bt3    |Bt3   |Bt2      | -0.0523289|   32| 43.5|              20|
|95  | 477055|10CKS030 |Bt2    |Bt2   |Bt3      | -0.0382704|   20| 40.5|              75|
|68  | 414936|10CKS009 |Bt2    |Bt2   |Bt1      | -0.0238213|   20| 31.0|              35|
|4   | 268820|07SKC003 |2Bt3   |Bt3   |Bt2      | -0.0222795|   38| 40.5|              10|
|102 | 477056|10CKS031 |CBt    |Bt3   |Bt2      | -0.0055023|   21| 43.0|              80|

The folowing table contains mean and standard deviations (values in parenthesis) of soil properties (clay content, total rock fragments, horizon mid-point) and silhouette widths summarized by GHL. Generalized horizons with the largest silhouette widths are the most internally consistent, while those with silhouette widths close to 0 are the least consistent.

```r
hz.eval$stats
```



|genhz    |clay         |mid           |total_frags_pct |sil.width   |
|:--------|:------------|:-------------|:---------------|:-----------|
|A        |14.27 (3.15) |3.27 (1.47)   |3.33 (3.39)     |0.51 (0.25) |
|Bt1      |19 (2.36)    |15.6 (5.19)   |11.87 (14.47)   |0.26 (0.2)  |
|Bt2      |24.13 (4.53) |32.93 (7.03)  |25.93 (24.65)   |0.04 (0.2)  |
|Bt3      |28.5 (5.17)  |57.8 (14.48)  |40.5 (24.71)    |0.12 (0.19) |
|Cr       |NaN (NA)     |76.75 (14.11) |0 (0)           |NaN (NA)    |
|R        |NaN (NA)     |142.15 (8.92) |0 (0)           |NaN (NA)    |
|not-used |15.8 (2.68)  |3.79 (5.07)   |1.36 (2.41)     |NaN (NA)    |

Plot profiles with only those horizons with silhouette widths less than 0 flagged.

```r
# add a column containing a color (red) that flags horizons with silhouette width less than 0
pedons$sil.flag <- ifelse(pedons$sil.width < 0, 'red', 'white')
par(mar=c(0,0,3,3))
plot(pedons, name='hzname', label='pedon_id', cex.names=0.8, axis.line.offset=-4, color='sil.flag')
```

<img src="figure/silhouette-sketch-flagged-1.png" title="plot of chunk silhouette-sketch-flagged" alt="plot of chunk silhouette-sketch-flagged" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Horizons with red colors warrant further investigation.</p><hr>

## Saving GHL to NASIS
GHL assignment as described thus far applies only to ***in-memory*** objects in the current R session. In order to use these GHL in further analysis or perform manual adjustments, the GHL need to be saved externally. If you are a NASIS user (working with data derived from NASIS) then the following code will create a text file that can be read by NASIS and stored in the `dspcomplayerid` field of the pedon horizon table. The "Update horizon group aggregations" calculation ("Calculations" -> "Pedon Horizon" -> "Update horizon group aggregations") expects a text file at the `C:/data/horizon_agg.txt` location on your system, containing phiid and generalized horizon labels.  

Note that some people prefer to adjust GHL assignments in R while others prefer to make adjustments after loading the data into NASIS.

```r
# clear-out any existing files
rules.file <- 'C:/data/horizon_agg.txt'
write.table(data.frame(), file=rules.file, row.names=FALSE, quote=FALSE, na='', col.names=FALSE, sep='|')
# extract horizon data
h <- horizons(pedons)
# strip-out 'not-used' genhz labels and retain horizon ID and genhz assignment
h <- h[which(h$genhz != 'not-used'), c('phiid', 'genhz')]
# append to NASIS import file
write.table(h, file=rules.file, row.names=FALSE, quote=FALSE, na='', col.names=FALSE, sep='|', append=TRUE)
```


# Concluding Remarks
In the next document, we will work with pedon data from NASIS and use our GHL to generate estimates of "low-rv-high" style range in characteristics. An alternative approach to assignment of GHL can be performed with the "Update comp layer id by hzname (horizon group Aggregation)" calculation. This calculation, however, performs only basic pattern matching and likely requires considerable manual intervention.


