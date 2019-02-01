## install.packages(c("aqp","soilDB","sharpshootR"))
## # install.packages("devtools") # if needed
## devtools::install_github(repo = 'ncss-tech/aqp', dependencies = FALSE, build = FALSE)
## devtools::install_github(repo = 'ncss-tech/soilDB', dependencies = FALSE, build = FALSE)
## devtools::install_github(repo = 'ncss-tech/sharpshootR', dependencies = FALSE, build = FALSE)
## devtools::install_github(repo = 'ncss-tech/soilReports', dependencies = FALSE, build = FALSE)
library(aqp)
library(soilDB)
library(sharpshootR)
library(soilReports)
## data("loafercreek")
## fetchNASIS() # default argument of fetchNASIS() gets pedons
## fetchNASIS_pedons()
## fetchNASIS_components()
## fetchKSSL(series = "Holland")
## fetchKSSL(mlra = "15")
## fetchKSSL(bbox = c(-119, 37, -120, 38))
## fetchSDA_component(WHERE = "compname = 'Mantree'")
## SDA_query(q = "SELECT * FROM legend AS l
##                 INNER JOIN mapunit mu ON l.lkey = mu.lkey
##                 INNER JOIN component co ON mu.mukey = co.mukey
##                WHERE areasymbol = 'CA630' AND compname = 'Mantree';")
## your.spc <- fetchOSD(soils = c("Holland", "Musick"))
## your.extended.list <- fetchOSD(soils = "Holland", extended = TRUE)
## 
## #the SPC can be found in the SPC attribute when extended=TRUE
## your.extended.list$SPC
## site(loafercreek)
## nrow(site(loafercreek))
## length(loafercreek) # length(spc) is an alias for nrow(site(spc))
## horizons(loafercreek)
## nrow(horizons(loafercreek))
## nrow(loafercreek) # nrow(spc) is an alias for nrow(horizons(spc))
## sub1 <- subsetProfiles(loafercreek, s='slope_field < 8')
## sub2 <- subsetProfiles(loafercreek, h = 'clay > 35')
## your.list <- c("07SKC016","10MJE016","10CKS034","S2000CA007012","06SMM013" )
## matches.your.list <- site(loafercreek)$pedon_id %in% your.list
## your.matching.indices <- which(matches.your.list)
## loafercreek.subset <- loafercreek[your.matching.indices, ]
## loafercreek$foo <- 1:nrow(site(loafercreek))
## 
## loafercreek$foohz <- 1:nrow(horizons(loafercreek))
## 
## # note that this one causes an error. Test your understanding: why?
## horizons(loafercreek)$foohz2 <- 1:nrow(site(loafercreek))
## site(loafercreek)$foo <- 2
## horizons(loafercreek)$foo <- 2
## site(your.spc) <- your.site.data
## horizons(your.spc) <- merge(horizons(your.spc), your.horizon.data, by = "hzid", all.x=TRUE)
data("loafercreek")
data("gopheridge")
## loafergopher <- union(list(loafercreek, gopheridge))
## plot(loafercreek[1:10, ], label = "pedon_id")
groupedProfilePlot(loafercreek[1:10, ], groups = "taxonname", 
                  label = "pedon_id", id.style="side",
                  axis.line.offset = -1.5, y.offset=25)
## coordinates(loafercreek) <- ~ x_std + y_std
## proj4string(loafercreek) <- '+proj=longlat +datum=WGS84'
## loafercreek.spdf <- as(loafercreek, 'SpatialPointsDataFrame')
## plot(loafercreek.spdf)
## rgdal::writeOGR(obj = loafercreek.spdf, layer = 'output.shapefile.name', driver = 'ESRI Shapefile')
