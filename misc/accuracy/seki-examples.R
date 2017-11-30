library(raster)
library(reshape2)
library(latticeExtra)
library(grid)


s <- stack('S:/NRCS/Archive_Dylan_Beaudette/SEKI/supervised_classification/2017-ideas/output/predictions/class-wise-probability-latest.tif')


## TODO: this stuff should be integrate into one of the AQP suite
source('local-functions.R')



x <- sampleRegular(s, 10000)

# NA in the results are areas without predictions
x <- na.omit(x)


#### the following required because of prior scaling of Pr to {0,100}
#### and randomForest generates 0sw

## 0-inflation only affects marginal Pr distribution investigation
## looking at the entire distribution is strongly influenced by 0-inflation

# rescale to probabilities
x <- x / 100.0
# add a very small fudge factor to remove 0
x[which(x==0)] <- 1e-15

####
####


class.labels <- toupper(letters[1:ncol(x)])

dimnames(x)[[2]] <- class.labels
d <- as.data.frame(x)
d$id <- rownames(d)

# uncertainty
H <- apply(x, 1, shannon.H, b=ncol(x))
CI <- apply(x, 1, confusion.index)

# prep for plotting
m <- melt(d, id.vars = c('id'), measure.vars = class.labels)
g <- make.groups(H, CI)


hist(apply(x, 1, max))
hist(x[x > 0.001])

plot(apply(x, 1, median), H, col=rgb(0,0,0,alpha = 0.125), pch=16)
plot(apply(x, 1, max), H, col=rgb(0,0,0,alpha = 0.125), pch=16)

plot(apply(x, 1, median), CI, col=rgb(0,0,0,alpha = 0.125), pch=16)
plot(apply(x, 1, max), CI, col=rgb(0,0,0,alpha = 0.125), pch=16)


bwplot(value ~ variable, data=m, subset= value > 0.001)


cols <- brewer.pal(9, 'Set1')
tps <- list(superpose.line=list(col=cols, lwd=1, alpha=0.85))

densityplot( ~ value , groups=variable, data=m, subset= value > 0.001,
                    pch=NA, xlim=c(-0.1, 1.1), scales=list(alternating=3, x=list(tick.number=5)), xlab='Class Probability',
                    strip=strip.custom(bg=grey(0.85)), auto.key=list(columns=5, lines=TRUE, points=FALSE),
                    par.settings=tps, panel=function(...) {
                      gs <- seq(0,1, by=0.1)
                      panel.abline(v=gs, lty=3, col='grey')
                      panel.densityplot(...)
                    })



# uncertainty metrics
cols <- wes_palette("Zissou")[c(1,5)]
tps <- list(superpose.line=list(col=cols, lwd=2, alpha=0.85))

densityplot( ~ data, groups=which, data=g, as.table=TRUE, pch=NA, auto.key=list(columns=2, lines=TRUE, points=FALSE), xlim=c(-0.1, 1.1), strip=strip.custom(bg=grey(0.85)), scales=list(alternating=3, y=list(rot=0), x=list(tick.number=5)), xlab='', par.settings=tps, panel=function(...) {
  gs <- seq(0,1, by=0.1)
  panel.abline(v=gs, lty=3, col='grey')
  panel.densityplot(...)
})


xyplot(H ~ CI, asp=1, scales=list(alternating=1, at=seq(0, 1, by=0.2)), xlim=c(-0.1, 1.1), ylim=c(-0.1, 1.1), xlab='Confusion Index', ylab='Shannon H',
       strip=strip.custom(bg=grey(0.85)), par.settings=list(plot.symbol=list(col='royalblue', pch=16, cex=0.85, alpha=0.25)), 
       panel=function(x, y, subscripts=subscripts, ...) {
         gs <- seq(0,1, by=0.1)
         med.H <- round(median(y), 2)
         med.CI <- round(median(x), 2)
         iqr.H <- round(IQR(y), 2)
         iqr.CI <- round(IQR(x), 2)
         cor.xy <- round(cor(x, y, method = 'spearman'), 2)
         ann <- paste0('H: ', med.H, ' (', iqr.H, ')\n', 
                       'CI: ', med.CI, ' (', iqr.CI, ')\n',
                       'cor: ', cor.xy)
         
         panel.abline(h=gs, v=gs, lty=3, col='grey')
         # panel.points(ss$CI, ss$Shannon.H, col=grey(0.85), alpha=0.1)
         panel.xyplot(x, y, subscripts=subscripts, ...)
         panel.abline(0, 1, lty=2)
         grid.text(ann, x = unit(0.05, 'npc'), unit(0.95, 'npc'), just = c('left', 'top'), gp = gpar(cex=0.75, font=2))
       })
