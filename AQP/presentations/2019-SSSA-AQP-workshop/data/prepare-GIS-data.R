

rs <- raster::stack('E:/gis_data/prism/effective_precipitation_800m.tif', 'e:/gis_data/prism/rain_fraction_mean_800m.tif')

e <- data.frame(id=profile_id(g), raster::extract(rs, as(g, 'SpatialPoints')), stringsAsFactors = FALSE)

write.csv(e, 'transect-GIS-data.csv', row.names = FALSE)

