#' Accessor to the data objects.
#'
#' @param ima a character string identifying for which array type retrieve data. Valid values are 'IlluminaHumanMethylation450k' and 'IlluminaHumanMethylationEPIC'.
#'
#' @return A RGChannelSetExtended object
#'
#' @export
#' 
#' @import ExperimentHub
#'
#' @examples
#' rgset_450k <- getCopyNeutralRGSet('IlluminaHumanMethylation450k')
#' rgset_450k
getCopyNeutralRGSet <- function(ima)
{
	valid <- c(paste0('IlluminaHumanMethylation', c('450k', 'EPIC')))
	if (! ima %in% valid)
		stop(paste('Valid array types are:', paste(valid, collapse=', ')))
	library(ExperimentHub)
	eh <- ExperimentHub()
	# load the reference dataset provided by the CopyNeutral450k package
	rgset <- loadResources(eh, package='CopyNeutralIMA', filterBy=ima)
	return(rgset)
}
