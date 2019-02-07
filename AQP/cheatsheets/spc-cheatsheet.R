## ---- eval=F, warning=F, message=F---------------------------------------
## install.packages(c("aqp","soilDB","sharpshootR"))

## ---- eval=F, warning=F, message=F---------------------------------------
## # install.packages("remotes") # if needed
## remotes::install_github(repo = 'ncss-tech/aqp', dependencies = FALSE, build = FALSE)
## remotes::install_github(repo = 'ncss-tech/soilDB', dependencies = FALSE, build = FALSE)
## remotes::install_github(repo = 'ncss-tech/sharpshootR', dependencies = FALSE, build = FALSE)

## ---- eval=F, warning=F, message=F---------------------------------------
## remotes::install_github(repo = 'ncss-tech/soilReports', dependencies = FALSE, build = FALSE)

## ---- eval=T, warning=F, message=F---------------------------------------
library(aqp)
library(soilDB)
library(sharpshootR)
library(soilReports)

## ---- eval=F, warning=F, message=F---------------------------------------
## data("loafercreek")

## ---- eval=F, warning=F, message=F---------------------------------------
## fetchNASIS() # default argument of fetchNASIS() gets pedons
## fetchNASIS_pedons()
## fetchNASIS_components()

## ---- eval=F, warning=F, message=F---------------------------------------
## fetchKSSL(series = "Holland")

## ---- eval=F, warning=F, message=F---------------------------------------
## fetchKSSL(mlra = "15")

## ---- eval=F, warning=F, message=F---------------------------------------
## fetchKSSL(bbox = c(-119, 37, -120, 38))

## ---- eval=F, warning=F, message=F---------------------------------------
## fetchSDA_component(WHERE = "compname = 'Mantree'")

## ---- eval=F, warning=F, message=F---------------------------------------
## SDA_query(q = "SELECT * FROM legend AS l
##                 INNER JOIN mapunit mu ON l.lkey = mu.lkey
##                 INNER JOIN component co ON mu.mukey = co.mukey
##                WHERE areasymbol = 'CA630' AND compname = 'Mantree';")

## ---- eval=F, warning=F, message=F---------------------------------------
## your.spc <- fetchOSD(soils = c("Holland", "Musick"))

## ---- eval=F, warning=F, message=F---------------------------------------
## your.extended.list <- fetchOSD(soils = "Holland", extended = TRUE)
## 
## #the SPC can be found in the SPC attribute when extended=TRUE
## your.extended.list$SPC

## ---- eval=F, warning=F, message=F---------------------------------------
## site(loafercreek)

## ----eval=F, warning=F, message=F----------------------------------------
## # overload of generic "length" function
## length(loafercreek)
## # this also works because each profile in the collection must have site level attributes
## nrow(site(loafercreek))

## ---- eval=F, warning=F, message=F---------------------------------------
## site(loafercreek)

## ---- eval=F, warning=F, message=F---------------------------------------
## horizons(loafercreek)

## ----eval=F, warning=F, message=F----------------------------------------
## # overload of generic "nrow" function
## nrow(loafercreek)
## nrow(horizons(loafercreek))

## ---- eval=F, warning=F, message=F---------------------------------------
## sub1 <- subsetProfiles(loafercreek, s='slope_field < 8')

## ---- eval=F, warning=F, message=F---------------------------------------
## sub2 <- subsetProfiles(loafercreek, h = 'clay > 35')

## ---- eval=F, warning=F, message=F---------------------------------------
## your.list <- c("07SKC016", "10MJE016", "10CKS034", "S2000CA007012", "06SMM013" )

## ---- eval=F, warning=F, message=F---------------------------------------
## matches.your.list <- site(loafercreek)$pedon_id %in% your.list

## ---- eval=F, warning=F, message=F---------------------------------------
## your.matching.indices <- which(matches.your.list)

## ---- eval=F, warning=F, message=F---------------------------------------
## loafercreek.subset <- loafercreek[your.matching.indices, ]

## ---- eval=F, warning=F, message=F---------------------------------------
## # be sure to specify the correct length
## loafercreek$foo <- 1:length(loafercreek)
## 
## # be sure to specify a vector that is the correct length
## loafercreek$foohz <- 1:nrow(loafercreek)
## 
## # note that this one causes an error. Test your understanding: why?
## horizons(loafercreek)$foohz2 <- 1:nrow(site(loafercreek))

## ---- eval=F, warning=F, message=F---------------------------------------
## # vectors shorter than length(SPC) are recycled
## site(loafercreek)$foo <- 2
## 
## # vectors shorter than nrow(SPC) are recycled
## horizons(loafercreek)$foo <- 2

## ---- eval=F, warning=F, message=F---------------------------------------
## site(your.spc) <- your.site.data

## ---- eval=F, warning=F, message=F---------------------------------------
## horizons(your.spc) <- merge(horizons(your.spc), your.horizon.data, by = "hzid", all.x=TRUE)

## ---- eval=T, warning=F, message=F---------------------------------------
data("loafercreek")
data("gopheridge")

## ---- eval=F, warning=F, message=F---------------------------------------
## loafergopher <- union(list(loafercreek, gopheridge))

## ---- eval=F, warning=F, message=F---------------------------------------
## plot(loafercreek[1:10, ], label = "pedon_id")

## ---- eval=T, warning=F, message=F, fig.width=9, fig.height=5------------
# adjust margins, units are inches, bottom, left, top, right
par(mar=c(1,1,0.5,1))

groupedProfilePlot(loafercreek[1:10, ], groups = "taxonname", 
                  label = "pedon_id", id.style="side",
                  axis.line.offset = -1.5, group.name.offset = -10)

## ---- eval=F, warning=F, message=F---------------------------------------
## coordinates(loafercreek) <- ~ x_std + y_std

## ---- eval=F, warning=F, message=F---------------------------------------
## proj4string(loafercreek) <- '+proj=longlat +datum=WGS84'

## ---- eval=F, warning=F, message=F---------------------------------------
## loafercreek.spdf <- as(loafercreek, 'SpatialPointsDataFrame')

## ---- eval=F, warning=F, message=F---------------------------------------
## plot(loafercreek.spdf)

## ---- eval=F, warning=F, message=F---------------------------------------
## rgdal::writeOGR(obj = loafercreek.spdf, layer = 'output.shapefile.name', driver = 'ESRI Shapefile')

