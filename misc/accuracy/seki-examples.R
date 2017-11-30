library(raster)
library(reshape2)
library(latticeExtra)
library(grid)

## TODO: this stuff should be integrate into one of the AQP suite
source('local-functions.R')

# # get some data from the latest SEKI model
# s <- stack('S:/NRCS/Archive_Dylan_Beaudette/SEKI/supervised_classification/2017-ideas/output/predictions/class-wise-probability-latest.tif')
# 
# # sub-sample
# x <- sampleRegular(s, 10000)
# # NA in the results are areas without predictions
# x <- na.omit(x)
# 
# # cache for later
# save(x, file='SEKI-example-data.Rda', compress = 'xz')

