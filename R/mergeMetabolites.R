#' mergeMetabolites
#'
#' This function merge metabolites identified more than once or merge intensities of differents derivatives.
#' @export
#' @param metabolites List of compound names as written in the pre_anno file.
#' @return None.
#' @importFrom grDevices dev.off pdf png tiff
#' @importFrom graphics boxplot grid legend par text
#' @importFrom methods as new
#' @importFrom stats cor sd t.test var
#' @importFrom utils read.csv write.csv write.table
#' @importFrom pheatmap pheatmap
#' @examples
#' \donttest{
#' \dontrun{
#' load(system.file("extdata", "n.RData", package = "PipMet"))
#' load(system.file("extdata", "metadata.RData", package = "PipMet"))
#' load(system.file("extdata", "colors.RData", package = "PipMet"))
#' heatmap(n, metadata, myDir = "~/", colors)
#' }
#' }
mergeMetabolites <- function(metabolites) {}
