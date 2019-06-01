
Make a new folder somewhere on a local disk and save the files posted [here](https://github.com/ncss-tech/ncss-tech.github.io/tree/master/AQP/presentations/2019-SSSA-AQP-workshop/data). These will be a helpful backup in case of degraded network availability during the workshop. If you have not used GitHub before, the following links may be more convenient.

   * [granite transect data](https://raw.githubusercontent.com/ncss-tech/ncss-tech.github.io/master/AQP/presentations/2019-SSSA-AQP-workshop/data/dahlgren-granitics.csv)
   * [andesite transect site data](https://raw.githubusercontent.com/ncss-tech/ncss-tech.github.io/master/AQP/presentations/2019-SSSA-AQP-workshop/data/rasmussen-andisitic-lahar-site.csv)
   * [andesite transect pedon data](https://raw.githubusercontent.com/ncss-tech/ncss-tech.github.io/master/AQP/presentations/2019-SSSA-AQP-workshop/data/rasmussen-andisitic-lahar.csv)
   * [sampled raster data](https://raw.githubusercontent.com/ncss-tech/ncss-tech.github.io/master/AQP/presentations/2019-SSSA-AQP-workshop/data/transect-GIS-data.csv)
   * [data processing R code](https://raw.githubusercontent.com/ncss-tech/ncss-tech.github.io/master/AQP/presentations/2019-SSSA-AQP-workshop/data/processing.R)
   * [raster sampling R code](https://raw.githubusercontent.com/ncss-tech/ncss-tech.github.io/master/AQP/presentations/2019-SSSA-AQP-workshop/data/prepare-GIS-data.R)
   * [MPS vs. slab R code](https://github.com/ncss-tech/ncss-tech.github.io/blob/master/AQP/presentations/2019-SSSA-AQP-workshop/data/mps-vs-slab.R)



New Features

http://soilmap2-1.lawr.ucdavis.edu/dylan/soilweb/sde/test.php?series=keswick


Inspired by Kyle's presentation last meeting (NASIS reports to support OSD work), I made some progress on a long-standing TODO item for my interest (obsession?) with a soil series database.

In short, it is now possible to retrieve summaries of several key data sources that can assist with the refinement of and comparison between soil series concepts. These summaries were developed from the current FY SSURGO snapshot, using "normalized" component names. This means component names were cleaned of any cruft (e.g. taxadjunct) and used as grouping criteria. For better / worse (let me know what you think) family level components were also included.

Summaries include:
.	hillslope position
.	geomorphic component
.	mountain position
.	parent material kind
.	parent material origin
.	MLRA "membership"
.	acreage / number of components
.	annual and monthly climate variables (from PRISM)
.	"siblings" ~ geographically associated soils (based on series the co-occur in the same map unit)
.	"cousins" (siblings of siblings)
.	competing series (via SC database)

You can get at all of these summaries for one or more series in R using new functionality found in the fetchOSD() and siblings() functions, from the soilDB package. Since these are just wrappers around a web-service, you can theoretically get the same summaries via URL with results returned as JSON-much like SDA.

Try it: 
https://casoilresource.lawr.ucdavis.edu/api/soil-series.php?q=siblings&s=bordengulch

Static versions of some of the summaries are posted here, in case you need all of them in one place:
https://github.com/ncss-tech/SoilTaxonomy/tree/master/databases


I have started to draft instructions on how to use the data, methods, key assumptions, possible uses, etc. here:

.	siblings: https://ncss-tech.github.io/AQP/soilDB/siblings.html
.	evaluating competing series: http://ncss-tech.github.io/AQP/soilDB/competing-series.html
.	climate summaries: http://ncss-tech.github.io/AQP/soilDB/series-climate-summary-eval.html
.	overview with examples: http://ncss-tech.github.io/AQP/soilDB/soil-series-query-functions.html




The project page is hosted here:
http://ncss-tech.github.io/AQP/

Many of the original ideas were outlined in this paper:
http://dx.doi.org/10.1016/j.cageo.2012.10.020


GitHub Pages:
https://github.com/ncss-tech/aqp

https://github.com/ncss-tech/soilDB

https://github.com/ncss-tech/sharpshootR


Function Reference:
http://ncss-tech.github.io/aqp/docs/reference/index.html

http://ncss-tech.github.io/soilDB/docs/reference/index.html

http://ncss-tech.github.io/sharpshootR/docs/reference/index.html




I wanted to share an excellent tutorial on aqp / soilDB concepts by one of my colleagues, Andrew Brown.

http://ncss-tech.github.io/AQP/demos/profileApply/loafercreek.html

This tutorial is based on one of the sample data sets included with the soilDB package, representing soil morphology data. These are all estimates made in the field, complete with some mistakes. The tutorial does a fine job of explaining some of the data-distilling strategies used by our soil scientists.

A fine prequel to this tutorial is the SoilProfileCollection tutorial. Just about everything you need to created / modify / crunch collections of soil data.

http://ncss-tech.github.io/AQP/aqp/aqp-intro.html


Alain asked about one of the demos I gave during the recent workshop in San Diego, related to soil series concepts. Here is a link to that (draft) document. This is my working copy of documentation pertaining to soil series concepts, supporting data, and comparisons among series.

http://ncss-tech.github.io/AQP/soilDB/soil-series-query-functions.html


For those interested in the SCAN/SNOTEL network, be sure to have a look at the fetchSCAN documentation here. This is a rough demonstration of ideas pertaining to getting / aggregating / displaying SCAN/SNOTEL data.

http://ncss-tech.github.io/AQP/soilDB/fetchSCAN-demo.html


Finally, if you have ideas, suggestions, or bugs to report please use the GitHub issues pages for each package:

https://github.com/ncss-tech/soilDB/issues

https://github.com/ncss-tech/aqp/issues

https://github.com/ncss-tech/sharpshootR/issues


