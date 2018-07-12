meta <- data.frame(
  Title = c("Copy neutral samples from IlluminaHumanMethylation450k arrays in GEO",
            "Copy neutral samples from IlluminaHumanMethylationEPIC arrays in GEO"),
  Description = c("51 IlluminaHumanMethylation450k samples. The provided controls consist of material from healthy individuals with nominally no copy number abberations.",
                  "13 IlluminaHumanMethylationEPIC samples. The provided controls consist of material from healthy individuals with nominally no copy number abberations."),
  BiocVersion = rep("3.8",2),
  Genome = rep("hg19",2),
  SourceType = rep("IDAT",2),
  SourceUrl = rep("http://www.ncbi.nlm.nih.gov/geo/",2),
  SourceVersion = format(Sys.time(), "%b %d %Y"),
  Species = rep("Homo sapiens",2),
  TaxonomyId = rep(9606,2),
  Coordinate_1_based = rep(TRUE,2),
  DataProvider = rep("GEO",2),
  Maintainer = rep("copyneutralima@compbio-dev.com",2),
  RDataClass = rep("RGChannelSetExtended",2),
  DispatchClass = rep("Rda",2),
  RDataPath = file.path("CopyNeutralIMA/data", c("IlluminaHumanMethylation450k.rda","IlluminaHumanMethylationEPIC.rda"))
)

# write metadata.csv file
write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
