library(aqp)
library(lattice)
library(colorspace)
library(Cairo)

hues <- c('5YR', '7.5YR', '10YR')
d <- expand.grid(hue=hues, value=2:8, chroma=seq(2,8,by=2), stringsAsFactors=FALSE)
d$color <- with(d, munsell2rgb(hue, value, chroma))

d.rgb <- with(d, munsell2rgb(hue, value, chroma, return_triplets=TRUE))
d.lab <- as(with(d.rgb, RGB(r,g,b)), 'LAB')
d <- data.frame(d, d.lab@coords)

png(file='static-figures/munsell-soil_colors-LAB.png', width=1100, height=700)
xyplot(value ~ factor(chroma) | factor(hue, levels=hues),
       main=list("Common Soil Colors - Annotated with LAB Coordinates", cex=1.75), layout=c(3,1), 
       par.settings=list(layout.heights=list(strip=1.5)),
       scales=list(cex=1.25, alternating=1), strip=strip.custom(bg=grey(0.85), par.strip.text=list(cex=1.5)),
       data=d, as.table=TRUE, subscripts=TRUE, xlab=list('Chroma', cex=1.25), ylab=list('Value', cex=1.25),
       panel=function(x, y, subscripts, ...) {
         panel.xyplot(x, y, pch=15, cex=12, col=d$color[subscripts])
         lab.text <- with(d[subscripts, ], paste(round(L), round(A), round(B), sep='\n'))
         panel.text(x, y, labels=lab.text, cex=1.25, col='white', font=2)
       }
)

dev.off()


