#' Accessor to the data objects.
#'
#' getCopyNeutralRGSet simplifies the access to the data of the package in ExperimentHub. The allowed values matched of those in the array definition of the RGChannelSet objects from package 'minfi'.
#' If 'ima' is set to 'IlluminaHumanMethylation450k' it will return the object with index 'EH1453' in ExperimentHub; if set to 'IlluminaHumanMethylationEPIC' it will return the object with index 'EH1454'.
#'
#' @param ima a character string specifying for which array type to retrieve data. Valid values are 'IlluminaHumanMethylation450k' and 'IlluminaHumanMethylationEPIC'.
#'
#' @return A \code{\link[minfi]{RGChannelSet-class}} object
#'
#' @export
#' 
#' @import ExperimentHub
#'
#' @examples
#' rgset_450k <- getCopyNeutralRGSet('IlluminaHumanMethylation450k')
#' rgset_450k
getCopyNeutralRGSet <- function(ima=c('IlluminaHumanMethylation450k', 'IlluminaHumanMethylationEPIC'))
{
    ima <- match.arg(ima)
	eh <- ExperimentHub()
	# load the reference dataset
	rgset <- loadResources(eh, package='CopyNeutralIMA', filterBy=ima)
	return(rgset[[1]])
}
