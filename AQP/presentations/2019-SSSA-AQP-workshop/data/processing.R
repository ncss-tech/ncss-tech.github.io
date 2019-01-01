
library(aqp)
library(soilDB)
library(sharpshootR)
library(raster)

x.g <- read.csv('dahlgren-granitics.csv', stringsAsFactors=FALSE)
x.a <- read.csv(file='rasmussen-andisitic-lahar.csv', stringsAsFactors=FALSE)
x.as.site <- read.csv(file='rasmussen-andisitic-lahar-site.csv', stringsAsFactors=FALSE)

x.g$soil_color <- with(x.g, munsell2rgb(hue, value, chroma))
x.a$soil_color <- with(x.a, munsell2rgb(hue, value, chroma))

depths(x.g) <- id ~ top + bottom
site(x.g) <- ~ elev + MAAT + MAP + geo + x + y

 
depths(x.a) <- id ~ top + bottom
site(x.a) <- ~ elev + precip + MAP + MAT + veg + Fe_d_to_Fe_t
site(x.a) <- x.as.site


data(mineralKing)


coordinates(x.g) <- ~ x + y
proj4string(x.g) <- '+proj=longlat +datum=NAD83'

coordinates(x.a) <- ~ x + y
proj4string(x.a) <- '+proj=longlat +datum=NAD83'

coordinates(mineralKing) <- ~ x + y
proj4string(mineralKing) <- '+proj=longlat +datum=NAD83'


x.g$transect <- rep('Granite', times=length(x.g))
x.a$transect <- rep('Andesite', times=length(x.a))
mineralKing$transect <- rep('Mineral King', times=length(mineralKing))


# x.g$HzD <- hzDistinctnessCodeToOffset(substr(x.g$hz_boundary, 0, 1))

g <- aqp::union(list(x.g, x.a, mineralKing))

## prepare GIS data 
# source('prepare-GIS-data.R')
## 

gis.data <- read.csv('transect-GIS-data.csv', stringsAsFactors = FALSE)
site(g) <- gis.data


g$Fe_o_to_Fe_d <- g$Fe_o / g$Fe_d

intersect(horizonNames(x.a), horizonNames(x.g))

par(mar=c(0,0,3,0))
plot(g, plot.order=order(g$elev))
plot(g, plot.order=order(g$elev), color='clay')
plot(g, plot.order=order(g$elev), color='Al_p')
plot(g, plot.order=order(g$elev), color='Si_o')
plot(g, plot.order=order(g$elev), color='CEC')
plot(g, plot.order=order(g$elev), color='Fe_o_to_Fe_d')

groupedProfilePlot(g, groups = 'transect', group.name.offset = -15)

g.new.order <- order(x.g$elev)
a.new.order <- order(x.a$elev)

par(mfcol=c(1, 2))
plot(x.g, name='name', plot.order=g.new.order, hz.distinctness.offset='HzD')
axis(1, at=1:length(x.g), labels=x.g$elev[g.new.order], line=-2)

plot(x.a, name='name', plot.order=a.new.order)
axis(1, at=1:length(x.a), labels=x.a$elev[a.new.order], line=-2)

pdf(file='RAD-transect.pdf', width=11, height=7)
par(mar=c(3.5,3.5,3,1))
plotTransect(x.g, 'elev', crs=CRS('+proj=utm +zone=11 +datum=NAD83'), grad.axis.title='Elevation (m)')
plotTransect(x.g, 'elev', crs=CRS('+proj=utm +zone=11 +datum=NAD83'), grad.axis.title='Elevation (m)', color='clay', col.label='Clay (%)')
plotTransect(x.g, 'elev', crs=CRS('+proj=utm +zone=11 +datum=NAD83'), grad.axis.title='Elevation (m)', color='sand', col.label='Sand (%)')
plotTransect(x.g, 'elev', crs=CRS('+proj=utm +zone=11 +datum=NAD83'), grad.axis.title='Elevation (m)', color='BS', col.label='Base Saturation (%)')
plotTransect(x.g, 'elev', crs=CRS('+proj=utm +zone=11 +datum=NAD83'), grad.axis.title='Elevation (m)', color='Fe_o', col.label='Oxalate-Fe (g/kg)')
dev.off()





m.order <- order(mineralKing$elev_field)

par(mar=c(1,1,2,1))

groupedProfilePlot(mineralKing, groups='taxonname', print.id=FALSE)

par(mfcol=c(1, 3))
plot(mineralKing, name='hzname', plot.order=m.order, label='taxonname')
axis(1, at=1:length(mineralKing), labels=mineralKing$elev_field[m.order], line=-2)


plot(x.g, name='name', plot.order=g.new.order, hz.distinctness.offset='HzD')
axis(1, at=1:length(x.g), labels=x.g$elev[g.new.order], line=-2)

plot(x.a, name='name', plot.order=a.new.order)
axis(1, at=1:length(x.a), labels=x.a$elev[a.new.order], line=-2)


g <- union(list(g, mineralKing))

g$elev <- ifelse(is.na(g$elev), g$elev_field, g$elev)
g$taxonname <- ifelse(is.na(site(g)$taxonname), site(g)$id, site(g)$taxonname)

par(mar=c(1,1,2,1))
plot(g, label='taxonname', plot.order=order(g$elev))
