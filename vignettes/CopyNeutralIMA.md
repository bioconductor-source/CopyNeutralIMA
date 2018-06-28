---
title: "The CopyNeutralIMA vignette."
author: "Moritz Przybilla, Xavier Pastor Hostench (Computational Oncology, German Cancer Research Center, Heidelberg, Germany)"
date: "2018-06-28"
output: rmarkdown::html_vignette
bibliography: REFERENCES.bib
vignette: >
  %\VignetteIndexEntry{CopyNeutralIMA}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
<style>
body {
text-align: justify}
</style>



# Contents 
1. Introduction
1. Description
1. Data
1. Example 

## 1. Overview
This document presents an overview of the *CopyNeutralIMA* package. The *CopyNeutralIMA* package provides reference samples for perfomring copy-number variation (CNV) analysis using Ilumina Infinium 450k or EPIC DNA methylation arrays. Several R bioconductor packages do genomic copy number profiling, including *conumee*, *ChAMP* or *CopyNumber450k* (deprecated). In order to achieve information about the copy number alterations, a set of copy neutral samples is required as a reference. The package *CopyNumber450kData*, usually used to provide the reference, is now deprecated. Additionally, there has never been an effort to provide reference samples for the EPIC arrays. To fill this gap of lacking reference samples, we here introduce the *CopyNeutralIMA* package. 

## 2. Description

In this package we provide a set of 51 InfiniumHumanMethylation450k and 13 InfiniumHumanMethylationEPIC samples. The provided samples constist of material from healthy individuals with nominally no copy number abberations. Users of *conumee* or other copy number profiling packages may use this data package as reference genomes.

## 3. Data

We selected the data from different studies accessible in the Gene Expression Omnibus. In particular, for 450k arrays samples from **GSE49618**, **GSE61441** and **GSE106089** were chosen. For EPIC arrays, normal or control samples from series **GSE86831/GSE86833**, **GSE98990** and **GSE100825** were chosen. The data can be found at [GEO](https://www.ncbi.nlm.nih.gov/geo/) using the respective identifier. 

## 4. Example using *CopyNeutralIMA* in combination with *conumee*

First, we load the data we want to analyse. To make it easier we will use the same example as in the **conumee** package.


```r
library(minfi)
#> Loading required package: BiocGenerics
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, cbind, colMeans,
#>     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
#>     Find, get, grep, grepl, intersect, is.unsorted, lapply,
#>     lengths, Map, mapply, match, mget, order, paste, pmax,
#>     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
#>     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
#>     tapply, union, unique, unsplit, which, which.max, which.min
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> Loading required package: SummarizedExperiment
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: DelayedArray
#> Loading required package: matrixStats
#> 
#> Attaching package: 'matrixStats'
#> The following objects are masked from 'package:Biobase':
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: 'DelayedArray'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
#> The following object is masked from 'package:base':
#> 
#>     apply
#> Loading required package: Biostrings
#> Loading required package: XVector
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:DelayedArray':
#> 
#>     type
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> Loading required package: bumphunter
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: locfit
#> locfit 1.5-9.1 	 2013-03-22
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
library(conumee)
#> Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
#> Loading required package: IlluminaHumanMethylation450kmanifest
#> Loading required package: IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#> 
#> Attaching package: 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19'
#> The following objects are masked from 'package:IlluminaHumanMethylation450kanno.ilmn12.hg19':
#> 
#>     Islands.UCSC, Locations, Manifest, Other,
#>     SNPs.132CommonSingle, SNPs.135CommonSingle,
#>     SNPs.137CommonSingle, SNPs.138CommonSingle,
#>     SNPs.141CommonSingle, SNPs.142CommonSingle,
#>     SNPs.144CommonSingle, SNPs.146CommonSingle,
#>     SNPs.147CommonSingle, SNPs.Illumina
#> Loading required package: IlluminaHumanMethylationEPICmanifest

RGSetTCGA <- read.450k.url()
#> downloading 6042324037_R05C02_Grn.idat
#>  - done.
#> downloading 6042324037_R05C02_Red.idat - done.
#> downloading 6042324037_R06C01_Grn.idat - done.
#> downloading 6042324037_R06C01_Red.idat - done.
RGSetTCGA
#> class: RGChannelSet 
#> dim: 622399 2 
#> metadata(0):
#> assays(2): Green Red
#> rownames(622399): 10600313 10600322 ... 74810490 74810492
#> rowData names(0):
#> colnames(2): 6042324037_R05C02 6042324037_R06C01
#> colData names(0):
#> Annotation
#>   array: IlluminaHumanMethylation450k
#>   annotation: ilmn12.hg19
```

After loading the data we normalize it:

```r
MSetTCGA <- preprocessIllumina(RGSetTCGA)
data(tcgaBRCA.sentrix2name)
```

We rename the samples with the known identifiers:

```r
sampleNames(MSetTCGA) <- tcgaBRCA.sentrix2name[sampleNames(MSetTCGA)]
MSetTCGA
#> class: MethylSet 
#> dim: 485512 2 
#> metadata(0):
#> assays(2): Meth Unmeth
#> rownames(485512): cg00050873 cg00212031 ... ch.22.47579720R
#>   ch.22.48274842R
#> rowData names(0):
#> colnames(2): TCGA-AR-A1AU-01A-11D-A12R-05
#>   TCGA-AR-A1AY-01A-21D-A12R-05
#> colData names(0):
#> Annotation
#>   array: IlluminaHumanMethylation450k
#>   annotation: ilmn12.hg19
#> Preprocessing
#>   Method: Illumina, bg.correct = TRUE, normalize = controls, reference = 1
#>   minfi version: 1.24.0
#>   Manifest version: 0.4.0
```

Now we load our control samples, from the same array type as our test samples and normalize them:




