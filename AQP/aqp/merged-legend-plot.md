---
output:
  html_document:
    theme: journal
    mathjax: null
    jquery: null
    smart: false
---





Plotting multiple SoilProfileCollection objects with a merged legend
===============================
D.E. Beaudette
<br>
2015-05-13
<br>
This document is based on `aqp` version 1.8-7 and `soilDB` version 1.5-5.


```r
library(aqp)
library(soilDB)
library(scales)
library(RColorBrewer)

# get OSD data for select series
osd <- fetchOSD(c("sobrante", "auburn"))

# get KSSL data for these series
kssl.1 <- fetchKSSL("sobrante")
kssl.2 <- fetchKSSL("auburn")

# color ramp function
cr <- colorRamp(rev(brewer.pal(10, "Spectral")))

# get the full range of clay
combined.data <- c(kssl.1$clay, kssl.2$clay)
combined.range <- range(combined.data, na.rm = TRUE)

# NA-padded value -> color mapping for full range of some horizon attribute
mapColor <- function(x, r, col.ramp) {
    c.rgb <- cr(scales::rescale(x, from = r))
    cc <- which(complete.cases(c.rgb))
    cols <- rep(NA, times = nrow(c.rgb))
    cols[cc] <- rgb(c.rgb[cc, ], maxColorValue = 255)
    return(cols)
}

# convert non-NA values into colors
kssl.1$.color <- mapColor(kssl.1$clay, combined.range, cr)
kssl.2$.color <- mapColor(kssl.2$clay, combined.range, cr)

# generate combined range / colors for legend
pretty.vals <- pretty(combined.data, n = 8)
legend.data <- list(legend = pretty.vals, col = rgb(cr(scales::rescale(pretty.vals)), maxColorValue = 255))

# assemble into list
spc.list <- list(osd, kssl.1, kssl.2)
```


```r
# note the use of .color for pre-computed value->color
par(mar = c(1, 1, 3, 3))
plotMultipleSPC(spc.list, group.labels = c("OSD", "Sobrante", "Auburn"), args = list(list(name = "hzname", 
    color = "soil_color", id.style = "side"), list(name = "hzn_desgn", color = ".color", label = "pedon_id", 
    id.style = "side"), list(name = "hzn_desgn", color = ".color", label = "pedon_id", id.style = "side")))

# legend is the combined legend
mtext(side = 3, text = "Clay Content (%)", font = 2, line = 1.6)
legend("bottom", legend = legend.data$legend, col = legend.data$col, bty = "n", pch = 15, horiz = TRUE, 
    xpd = TRUE, inset = c(0, 0.99))
```

<img src="figure/plot-color-1.png" title="plot of chunk plot-color" alt="plot of chunk plot-color" style="display: block; margin: auto;" /><p class="caption" style="font-size:85%; font-style: italic; font-weight: bold;">Merged clay content legends.</p><hr>



