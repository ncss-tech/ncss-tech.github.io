knitr::opts_chunk$set(echo = TRUE)

library(aqp)
library(soilDB)

# get components from Yosemite National Park for a sandy-skeletal, cryic, inceptisol
my.spc <- fetchSDA(WHERE="compname = 'Glacierpoint'")

# copy contents of diagnostics slot
diagz <- diagnostic_hz(my.spc)

# create a subset SPC with just profiles that have a cambic in their diagnostics
my.spc.sub <- my.spc[profile_id(my.spc) %in% 
                       diagz[grepl(diagz$featkind, pattern="[Cc]ambic"), idname(my.spc)],]

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

# apply function has_sandy_cambic() to each profile
my.spc.sub$sandy_cambic <- profileApply(my.spc.sub, has_sandy_cambic)

# do any components have a sandy cambic?
if(any(my.spc.sub$sandy_cambic)) {
  print(paste0("The following components have sandy cambics: ",
        paste0(profile_id(my.spc.sub)[my.spc.sub$sandy_cambic], collapse=",")))
} else {
  print("All components with cambics have appropriate textures.")
}
