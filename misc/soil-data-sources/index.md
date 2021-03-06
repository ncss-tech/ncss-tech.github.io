---
title: "USDA-NRCS / NCSS Soil Data Sources"
author: "D.E. Beaudette and Soil Survey Staff"
date: "April 2019"
output:
  html_vignette:
    mathjax: null
    jquery: null
    smart: no
---


## Digital Soil Survey Data

[Main link](http://www.nrcs.usda.gov/wps/portal/nrcs/main/soils/survey/geo/)

   * [SSURGO](ssurgo.html) (1:24,000 scale) detailed soil survey
   * STATSGO (1:250,000 scale) generalized soil survey
   

## Official Soil Series Descriptions (OSDs)


## Laboratory Characterization Data

These pages provide search and download facilities for soil characterization data.

  * [Basic Query / Home Page](http://ncsslabdatamart.sc.egov.usda.gov/)
  * [Advanced Query](http://ncsslabdatamart.sc.egov.usda.gov/advquery.aspx)


## Soil Data Access (SDA)


## Official Series Descriptions (OSD)

These pages provide search tools for the OSD pages and associated Soil Classification (SC) database.

  * [OSD Search](https://soilseries.sc.egov.usda.gov/osdname.aspx)
  * [Multiple OSD Search](https://soilseries.sc.egov.usda.gov/osdlist.aspx)
  * [OSD Search via Soil Classification (SC) Database](https://soilseries.sc.egov.usda.gov/osdquery.aspx)

Search the OSD documents using pattern matching within the entire OSD, or limited to sections of the OSD. See each page below for examples on how to create search patterns.

  * [search by section](https://casoilresource.lawr.ucdavis.edu/osd-search/)


## SCAN / SNOTEL

Above and below ground sensor data. The [interactive map](http://www.wcc.nrcs.usda.gov/webmap/index.html#elements=&networks=SCAN&states=!&counties=!&hucs=&minElevation=&maxElevation=&elementSelectType=all&activeOnly=true&hucLabels=false&stationLabels=&overlays=&hucOverlays=&mode=stations&openSections=dataElement,parameter,date,elements,location,networks&controlsOpen=true&popup=&base=esriNgwm&lat=45.06&lon=-101.95&zoom=4&dataElement=PREC&parameter=PCTAVG&frequency=DAILY&duration=null&customDuration=&dayPart=E&year=2016&month=6&day=22&monthPart=E) is the simplest way to search for data.

   * [SCAN](http://www.wcc.nrcs.usda.gov/scan/)
   * [SNOTEL](http://www.wcc.nrcs.usda.gov/snow/)
   

## Block Diagrams

This is a collection of block diagrams from historic soil survey documents that have been scanned and indexed.

[Historic Block Diagram Archive](http://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/geo/?cid=nrcs142p2_054317)
   


## SoilWeb: UC Davis / NRCS Collaboration

SoilWeb is an interface to the SSURGO, KSSL, SC, Block Diagram, and OSD databases, updated with each (fiscal year) SSURGO release.

[SoilWeb Home](https://casoilresource.lawr.ucdavis.edu/soilweb-apps)


## AQP

[Project Homepage](http://ncss-tech.github.io/AQP/)


### R-Based Interfaces

 * SDA
    + [`SDA_query`](http://ncss-tech.github.io/soilDB/docs/reference/SDA_query.html)
    + ['SDA_query_features`](http://ncss-tech.github.io/soilDB/docs/reference/SDA_query_features.html)

  * SSURGO/KSSL via SoilWeb
    + [`fetchKSSL`](http://ncss-tech.github.io/soilDB/docs/reference/fetchKSSL.html)
    + [`fetchOSD`](http://ncss-tech.github.io/soilDB/docs/reference/fetchOSD.html)
    + [`siblings`](http://ncss-tech.github.io/soilDB/docs/reference/siblings.html)
    + [`OSDquery`](http://ncss-tech.github.io/soilDB/docs/reference/OSDquery.html) 
    + [`seriesExtent`](http://ncss-tech.github.io/soilDB/docs/reference/seriesExtent.html)
    
  * SCAN/SNOTEL
    + [`fetchSCAN`](http://ncss-tech.github.io/soilDB/docs/reference/fetchSCAN.html)
    + [`SCAN_SNOTEL_metadata`](http://ncss-tech.github.io/soilDB/docs/reference/SCAN_SNOTEL_metadata.html)

