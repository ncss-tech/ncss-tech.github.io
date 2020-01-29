#' ---
#' title: "Particle Size Control Section Validation with aqp::estimatePSCS()"
#' author: "Andrew G. Brown; andrew.g.brown@usda.gov"
#' date: "Last updated: 2020/01/25"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ---- echo=FALSE, message=FALSE------------------------------------------
# load aqp and soilDB libraries
library(aqp)
library(soilDB)

# load example dataset mineralKing from CA792
data("mineralKing", package = 'soilDB')
spc <- mineralKing
hzidname(spc) <- "phiid"

## alternately, get your selected set!
#spc <- fetchNASIS()

## or, get KSSL pedons (you will need to adjust variable names)
#spc <- fetchKSSL(series="cecil")

# optional: plot all profiles, not just the ones with errors
plot_all=FALSE

#' 
## ---- echo=FALSE, message=FALSE------------------------------------------
# Use `aqp` to calculate PSCS

# We make an estimate of what the particle size control section should be for each profile in the SoilProfileCollection spc. We then make a table of pedon record IDs and estimated PSCS top and bottom depth.
# estimate PSCS for each profile, return data.frame result that
# can be left-joined via site<- setter
site(spc) <- profileApply(spc, frameify=TRUE, function(p) { 
    pscs <- estimatePSCS(p)
    return(data.frame(peiid=profile_id(p),
                      calcpscstop=pscs[1],
                      calcpscsbot=pscs[2]))
  })


#' 
#' 
#' 
## ---- echo=FALSE, message=FALSE, fig.align='center'----------------------
# Compare Stored v.s. Calculated

# `psctopdepth` and `pscbottomdepth` contain the top and bottom depths picked by `soilDB:::.pickBestTaxHistory` for `fetchNASIS`.
# `calcpscstop` and `calcpscsbot` are the values estimated by `aqp` based on applying control section key criteria.

res <- site(spc)[,c('pedon_id','psctopdepth','calcpscstop',
                             'pscbotdepth','calcpscsbot')]

# determine whether top, bottom, both match
res$topmatch <- res[,'calcpscstop'] == res[,'psctopdepth']
res$botmatch <- res[,'calcpscsbot'] == res[,'pscbotdepth']
res$match <- res$topmatch & res$botmatch

# make table of non-matchingS results
idx <- which(is.na(res$match) | !res$match)

# always make a table of non-matching or non-populated
if(length(idx))
  knitr::kable(res[idx,], row.names = FALSE)

# yay!
if(!length(idx))
  print("All stored PSCS top and bottom depths match calculated.")

# if all profiles are to be plotted, plot them
if(plot_all)
  idx <- 1:length(spc)

#' 
## ---- echo=FALSE, message=FALSE, results="asis", fig.width=12, fig.asp = .62----
if(length(idx)) {
  
  cat("Blue bracket depicts calculated PSCS, and red line in center of profile is stored PSCS.")
  
  chunks <- makeChunks(idx, size = 5)
  names(chunks) <- idx
  
  chunks
  
  for(i in 1:max(chunks)) {
    j <- as.numeric(names(chunks)[which(chunks == i)])
    
    #subset
    spc.sub <- spc[j,]
    
    # make a profile plot
    plotSPC(spc.sub, 
            label="pedon_id",
            axis.line.offset = -2,
            cex.names=1.1)
    
    addBracket(data.frame(peiid=profile_id(spc.sub),
                          top=spc.sub$calcpscstop,
                          bottom=spc.sub$calcpscsbot),
                          col="blue", lwd=3,
                          tick.length = 0.1, offset = -0.1)
    
    addBracket(data.frame(peiid=profile_id(spc.sub),
                          top=spc.sub$psctopdepth,
                          bottom=spc.sub$pscbotdepth),
                          col="red", lwd=3,
                          tick.length = 0, offset = 0)
  }
}
  

