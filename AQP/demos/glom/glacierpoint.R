#' ---
#' title: "Glacierpoint - glom() Demo"
#' subtitle: 'v0.2.0 - last updated 2020/01/25'
#' author: "Andrew Brown; andrew.g.brown@usda.gov"
#' version: "0.2.0"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' In this demo, we use the _aqp_ functions `glom()` and `profileApply()` to analyze a depth-subset of horizons from each profile in a _SoilProfileCollection_.
#' 
#' `glom()` returns the horizon(s) that intersect the specified depth [_z1_] or depth interval [_z1_, _z2_].
#' 
#' ### Checking for cambic horizon criteria
#' 
#' First, the script queries components by component name, then gets records populated in the _Component Diagnostic Features_ (`codiagfeatures`) table. 
#' 
#' `glom()` is used to get all horizons that are _within the cambic depth range_. A custom function `has_sandy_cambic()` is defined to use `glom()` on individual profiles (`p`). 
#' 
#' A regular expression pattern is used to check the `glom()`-ed horizon textures field (`p$texture`) for stratified and/or sandy texture classes that are _not_ allowed in cambic horizons.
#' 
#' If either stratified or sandy textures are found, the function returns `TRUE` for that profile. Otherwise, the function returns `FALSE`.
#' 
#' ***
#' 
#' This demonstration script loads some components from Yosemite National Park (_CA790_). In this simple example, we look at some components correlated to _Glacierpoint_.
#' 
#' The _Glacierpoint_ series (_Sandy-skeletal, isotic Xeric Dystrocryepts_) is defined as having __either__ an Umbric epipedon __OR__ a cambic horizon. It appears that cambic horizons were populated even for some of the components whose RV textures do not meet cambic criteria.
#' 
## ----message=F, warning=F------------------------------------------------
library(aqp)
library(soilDB)

# get components from Yosemite National Park for a sandy-skeletal, cryic, inceptisol
my.spc <- fetchSDA(WHERE="compname = 'Glacierpoint'")

#' 
#' Now, access the SPC `@diagnostics` slot to figure out which components have cambics, and then make a subset SPC `my.spc.sub` with _just those components_.
#' 
## ----message=F, warning=F------------------------------------------------
# copy contents of diagnostics slot
diagz <- diagnostic_hz(my.spc)

# create a subset SPC with just profiles that have a cambic in their diagnostics
my.spc.sub <- my.spc[profile_id(my.spc) %in% 
                       diagz[grepl(diagz$featkind, pattern="[Cc]ambic"), idname(my.spc)],]

#' 
#' A total of `r length(my.spc.sub)` components out of `r length(my.spc)` have cambic horizons. 
#' 
#' Let's inspect those `r length(my.spc.sub)` components visually. 
#' 
#' Use `plotSPC()` to show sand fraction percentage by horizon -- with cambic top and bottom depths shown with brackets using `aqp::addDiagnosticBracket()`.
#' 
## ----message=F, warning=F------------------------------------------------
# set margins
par(mar=c(0,1,3,1))

# plot spc using total sand RV % as color
plotSPC(my.spc.sub, 
        color='sandtotal_r', 
        width=0.15, 
        cex.names = 1)

# add red brackets showing cambic horizon top and bottom depth
addDiagnosticBracket(my.spc.sub, 
                     kind='Cambic horizon', 
                     top='featdept_r',
                     bottom='featdepb_r',
                     col='red', 
                     offset=0,
                     lwd=10)

#' 
#' ***
#' 
#' Now that we have the data and have done some preliminary inspection, let's define a regular expression pattern and a function to check our cambic criteria. 
#' 
#' Define a custom function `has_sandy_cambic()` to test horizons _within the cambic depth range_ for each profile. The default argument `sandy.texture.pattern` matches stratified or sandy textures. 
#' 
## ----message=F, warning=F------------------------------------------------
# has_sandy_cambic() has one required argument p, which is a single profile SPC 
has_sandy_cambic <- function(p, sandy.texture.pattern = "SR|-S$|^S$|COS$|L[^V]FS$|[^L]VFS$|LS$|LFS$") {
  d <- diagnostic_hz(p)
  
  # match the pattern cambic or Cambic
  d.cambic.idx <- grep(d$featkind, pattern="[Cc]ambic")
  cambic_top <- d[d.cambic.idx, 'featdept_r']
  cambic_bot <- d[d.cambic.idx, 'featdepb_r']
  
  # if top or bottom depth are NULL there is no cambic
  if(!is.null(cambic_top) | !is.null(cambic_bot)) {
    
    # glom horizons in p from cambic_top to cambic_bot 
    cambic_horizons <- glom(p, cambic_top, cambic_bot, as.data.frame = TRUE)
    
    print(paste0("COKEY: ",unique(p$cokey)," has cambic (",
                 cambic_top," - ",cambic_bot," cm) with ", paste0(cambic_horizons$hzname, collapse=","), 
                 " horizons of texture class: ", paste0(cambic_horizons$texture, collapse=", ")))
  
    # use a regular expression to match the sandy textures
    if(any(grepl(cambic_horizons$texture, 
                 pattern = sandy.texture.pattern, 
                 ignore.case = TRUE))) {
      return(TRUE) 
    }  
  }
  return(FALSE)
}

#' 
#' Apply the function to each profile.
#' 
## ----message=F, warning=F------------------------------------------------
# apply function has_sandy_cambic() to each profile
my.spc.sub$sandy_cambic <- profileApply(my.spc.sub, has_sandy_cambic)

#' 
#' Inspect results. Print out the offending unique profile IDs (`cokey`).
#' 
## ----message=F, warning=F------------------------------------------------
# do any components have a sandy cambic?
if(any(my.spc.sub$sandy_cambic)) {
  print(paste0("The following components have sandy cambics: ",
        paste0(profile_id(my.spc.sub)[my.spc.sub$sandy_cambic], collapse=",")))
} else {
  print("All components with cambics have appropriate textures.")
}

#' 
#' ----------------------------
#' This document is based on `aqp` version `r utils::packageDescription("aqp", field="Version")`, `soilDB` version `r utils::packageDescription("soilDB", field="Version")`, and `sharpshootR` version `r utils::packageDescription("sharpshootR", field="Version")`.
