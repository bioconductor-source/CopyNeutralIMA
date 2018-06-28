## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment='#>')

## ------------------------------------------------------------------------
library(minfi)
library(conumee)

RGSetTCGA <- read.450k.url()
RGSetTCGA

## ------------------------------------------------------------------------
MSetTCGA <- preprocessIllumina(RGSetTCGA)
data(tcgaBRCA.sentrix2name)

## ------------------------------------------------------------------------
sampleNames(MSetTCGA) <- tcgaBRCA.sentrix2name[sampleNames(MSetTCGA)]
MSetTCGA

