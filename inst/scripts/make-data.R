## The following code is used to download normal tissue samples from GEO. 
## Afterwards, the data are processed using minfi and 
## subsequently stored in an RGChannelSet which is exported to the data folder. 

library(minfi)
library(GEOquery)
library(devtools)
library(GEOmetadb)

# Create paths to the most used directories
pathToRawData <- file.path(getwd(), "data-raw")

# Generate data.frame with metadata for selected samples #
## connect with GEOmetadb
sqlite_file <- file.path(pathToRawData, 'GEOmetadb.sqlite')
if (!file.exists(sqlite_file)) getSQLiteFile(destdir=pathToRawData)
con <- dbConnect(SQLite(), sqlite_file)

select_string <- "select id, gsm, series_id, title, description, source_name_ch1, characteristics_ch1 from gsm where series_id in"

## store series_id, title, gpl and gsm from the previously investigated datasets in a data.frame
meta <- dbGetQuery(con, paste(select_string, "('GSE49618', 'GSE61441', 'GSE32148', 'GSE106089')"))
rownames(meta) <- meta$gsm

## produce two vectors with GSM identifiers for EPIC and 450k arrays
meta <- meta[grep('[Nn]ormal', meta$description),]
gsm <- meta$gsm[meta$series_id=='GSE49618']
gsm_other <- meta$gsm[!meta$gsm %in% gsm]
gsm_450k <- c(gsm, gsm_other[seq(1, length(gsm_other), length.out=30)])
meta450k <- meta[gsm_450k,]

# Download EPIC identifiers and metadata
meta <- dbGetQuery(con, paste(select_string, "('GSE86831,GSE86833', 'GSE100825', 'GSE98990')"))

# produce two vectors with GSM identifiers for EPIC and 450k arrays
multipleregex <- c("(Normal)", "(CTR)", "(NAF)", "(Guthrie)", "(Genomic)")
meta <- meta[grep(paste(multipleregex, collapse = "|"), meta$title),]
gsm_epic <- meta$gsm
## disconnect from GEOmetadb
dbDisconnect(con)

## Produce RGChannelSetExdended Object from GEO #
make_GEO_RGSet <- function(gsm, pdata, tmpdir) {
	## Load sample files from GEO
	dir.create(tmpdir)
	sapply(gsm, getGEOSuppFiles, baseDir = tmpdir)

	## save all files with "idat.gz" in a list of idatfiles for EPIC and 450k data ##
	idat <- list.files(tmpdir, pattern = "idat.gz$", full = TRUE, recursive = TRUE)

	## apply the gunzip function on all the files in the lists and overwrite them
	sapply(idat, gunzip, overwrite = TRUE)

	## RGChannelSet for data from Infinium Methlyation450k
	rgset <- read.metharray.exp(tmpdir, targets=NULL, recursive=T, extended=T)
	smplNames <- gsub('_.*', '', sampleNames(rgset))
	sampleNames(rgset) <- smplNames
	pData(rgset) <- DataFrame(pdata)
	return(rgset)
}

IlluminaHumanMethylation450k <- make_GEO_RGSet(gsm_450k, meta450k, file.path(pathToRawData, '450k'))
IlluminaHumanMethylationEPIC <- make_GEO_RGSet(gsm_epic, meta, file.path(pathToRawData, 'epic'))

devtools::use_data(IlluminaHumanMethylation450k)
devtools::use_data(IlluminaHumanMethylationEPIC)
