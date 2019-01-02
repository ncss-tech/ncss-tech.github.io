

rs <- raster::stack(
  'E:/gis_data/prism/effective_precipitation_800m.tif', 
  'E:/gis_data/prism/rain_fraction_mean_800m.tif',
  'E:/gis_data/prism/final_MAAT_800m.tif',
  'E:/gis_data/prism/final_MAP_mm_800m.tif'
  )

names(rs) <- c('effective.ppt_800', 'rain.fraction_800', 'MAAT_800', 'MAP_800')

e <- data.frame(id=profile_id(g), raster::extract(rs, as(g, 'SpatialPoints')), stringsAsFactors = FALSE)

write.csv(e, 'transect-GIS-data.csv', row.names = FALSE)

