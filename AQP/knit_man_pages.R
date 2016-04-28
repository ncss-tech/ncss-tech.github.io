library(knitr)

opts_chunk$set(message=FALSE, warning=FALSE, background='#F7F7F7', dpi=100, fig.align='center', fig.height=6, fig.width=6, dev='png', tidy=FALSE)
options(width=100, stringsAsFactors=FALSE)

setwd('E:/working_copies/aqp/www/aqp-html-manual/')
knit_rd('aqp')
unlink('figure', recursive=TRUE)


setwd('E:/working_copies/aqp/www/soilDB-html-manual/')
knit_rd('soilDB')
unlink('figure', recursive=TRUE)


setwd('E:/working_copies/aqp/www/sharpshootR-html-manual/')
knit_rd('sharpshootR')
unlink('figure', recursive=TRUE)
