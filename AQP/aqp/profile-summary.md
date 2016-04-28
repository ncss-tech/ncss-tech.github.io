---
output:
  html_document:
    theme: journal
    mathjax: null
    jquery: null
    smart: false
---





Profile Summaries
===============================
D.E. Beaudette
<br>
2015-04-21
<br>
This document is based on `aqp` version 1.8-6 and `soilDB` version 1.5-3.


# Introduction



```r
# load required libaries
library(soilDB)
library(lattice)
library(reshape2)
library(plyr)

# get all pedons from the selected set
x <- fetchNASIS(rmHzErrors = FALSE, nullFragsAreZero = FALSE)

# init vector of taxonnames to keep
soils <- c("Whiterock", "Copperopolis", "Shabarudo", "Dunstone", "Loafercreek", "Argonaut", "Crimeahouse", 
    "Hennekenot")

# convert vector of taxonnames into REGEX pattern for matching
pat <- paste0(soils, collapse = "|")

# subset pedons that match our REGEX pattern
idx <- grep(pat, x$taxonname, ignore.case = TRUE)
x <- x[idx, ]

# normalize taxonname via REGEX matching
for (i in soils) x$taxonname[grep(i, x$taxonname, ignore.case = TRUE)] <- i

# aggregate data by normalized taxonname, via slice-wise mean
a.colors <- slab(x, taxonname ~ m_r + m_g + m_b + clay + phfield + total_frags_pct, slab.fun = mean, 
    na.rm = TRUE)

# throw out aggregate data that are deeper than 150cm
a.colors <- subset(a.colors, subset = bottom < 150)

# convert long -> wide format
x.colors <- dcast(a.colors, taxonname + top + bottom ~ variable, value.var = "value")

# check
head(a.colors)
```



|variable |taxonname |     value| top| bottom| contributing_fraction|
|:--------|:---------|---------:|---:|------:|---------------------:|
|m_r      |Argonaut  | 0.3979233|   0|      1|             0.6666667|
|m_r      |Argonaut  | 0.3979233|   1|      2|             0.6666667|
|m_r      |Argonaut  | 0.4123208|   2|      3|             0.8333333|
|m_r      |Argonaut  | 0.4140862|   3|      4|             0.8333333|
|m_r      |Argonaut  | 0.4110800|   4|      5|             0.9166667|
|m_r      |Argonaut  | 0.4110800|   5|      6|             0.9166667|

```r
head(x.colors)
```



|taxonname | top| bottom|       m_r|       m_g|       m_b|     clay|  phfield| total_frags_pct|
|:---------|---:|------:|---------:|---------:|---------:|--------:|--------:|---------------:|
|Argonaut  |   0|      1| 0.3979233| 0.2492589| 0.1579567| 18.33333| 5.900000|        3.000000|
|Argonaut  |   1|      2| 0.3979233| 0.2492589| 0.1579567| 18.33333| 5.900000|        3.000000|
|Argonaut  |   2|      3| 0.4123208| 0.2575223| 0.1562377| 17.72727| 5.928571|        5.444444|
|Argonaut  |   3|      4| 0.4140862| 0.2569061| 0.1526347| 18.18182| 5.928571|        6.000000|
|Argonaut  |   4|      5| 0.4110800| 0.2567467| 0.1532129| 18.72727| 5.957143|        5.600000|
|Argonaut  |   5|      6| 0.4110800| 0.2567467| 0.1532129| 18.72727| 5.957143|        5.600000|

```r
# composite RGB triplets into an R-compatible color note that missing colors must be padded with NA
x.colors$soil_color <- NA
not.na <- which(complete.cases(x.colors[, c("m_r", "m_g", "m_b")]))
x.colors$soil_color[not.na] <- with(x.colors[not.na, ], rgb(m_r, m_g, m_b, maxColorValue = 1))

# aggregate bedrock depth probabilty by taxonname at 90% level of confidence
dp <- aggregateSoilDepth(x, "taxonname", crit.prob = 0.9)

# init a new SoilProfileCollection from aggregate data
depths(x.colors) <- taxonname ~ top + bottom
# join-in our depth to contact data
site(x.colors) <- dp

# generate index for new ordering
new.order <- match(c("Whiterock", "Copperopolis", "Shabarudo", "Dunstone", "Loafercreek", "Argonaut", 
    "Crimeahouse", "Hennekenot"), profile_id(x.colors))
```


```r
par(mar = c(1, 0, 3, 0))
plot(x.colors, divide.hz = FALSE, name = "", plot.order = new.order, col.label = "Soil Color", lwd = 1.25, 
    axis.line.offset = -6, cex.depth.axis = 1, cex.id = 1)
addBracket(x.colors$soil.top, x.colors$soil.bottom, col = "black", label = "P(soil >= 90%)", label.cex = 0.85)
title("Aggregate Soil Properties (mean)")
```

<img src="figure/plot-color-1.png" title="plot of chunk plot-color" alt="plot of chunk plot-color" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">sdfsdfsdf</p><hr>


```r
par(mar = c(1, 0, 3, 0))
plot(x.colors, divide.hz = FALSE, color = "clay", name = "", plot.order = new.order, col.label = "Clay Content (%)", 
    lwd = 1.25, axis.line.offset = -6, cex.depth.axis = 1, cex.id = 1)
addBracket(x.colors$soil.top, x.colors$soil.bottom, col = "black", label = "P(soil >= 90%)", label.cex = 0.85)
```

<img src="figure/plot-clay-1.png" title="plot of chunk plot-clay" alt="plot of chunk plot-clay" style="display: block; margin: auto;" />

```r
plot(x.colors, divide.hz = FALSE, color = "phfield", name = "", plot.order = new.order, col.label = "pH", 
    lwd = 1.25, axis.line.offset = -6, cex.depth.axis = 1, cex.id = 1)
addBracket(x.colors$soil.top, x.colors$soil.bottom, col = "black", label = "P(soil >= 90%)", label.cex = 0.85)
```

<img src="figure/plot-clay-2.png" title="plot of chunk plot-clay" alt="plot of chunk plot-clay" style="display: block; margin: auto;" />

```r
plot(x.colors, divide.hz = FALSE, color = "total_frags_pct", name = "", plot.order = new.order, col.label = "Total Rock Fragment Volume (%)", 
    lwd = 1.25, axis.line.offset = -6, cex.depth.axis = 1, cex.id = 1)
addBracket(x.colors$soil.top, x.colors$soil.bottom, col = "black", label = "P(soil >= 90%)", label.cex = 0.85)
```

<img src="figure/plot-clay-3.png" title="plot of chunk plot-clay" alt="plot of chunk plot-clay" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">sdfsdfsdf</p><hr>

