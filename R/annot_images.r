#' Annotation images
#'
#' This function creates network, hierarchy and heatmap plots of the spectra (annotated or not).
#' @export
#' @param pslist List of spectra.
#' @param myDir Path to the directory of work.
#' @param pictures Logical. If pictures should be plotted or not.
#' @return A list with (1) 'pseudoespectrum' object annotated resulted from 'addAnnotations' (CluMSID package) and (2) a table with annotation.
#' @importFrom grDevices dev.off pdf png tiff
#' @importFrom graphics par
#' @importFrom utils read.csv
#' @importFrom CluMSID addAnnotations distanceMatrix networkplot HCplot
#' @examples
#' \dontrun{
#' annot_images(pslist, myDir)
#' }
#'
annot_images <- function(pslist, myDir, pictures = c(TRUE, FALSE)) {
  r <- read.csv("pre_anno.csv", sep = ",", na.string = c("NA, "))
  if (length(colnames(r)) <= 1) {
    r <- read.csv("pre_anno.csv", sep = ";")
  }

  ### add annotations from 'pre_anno.csv' file (if there is annotation)
  if (sum(is.na(r$annotation)) == nrow(r)) {
    apslist <- pslist
  } else {
    apslist <- addAnnotations(featlist = pslist, annolist = r)
  }
  pseudodistmat <- distanceMatrix(apslist, mz_tolerance = 0.02) ### calculates distance matrix; takes a while

  if (pictures == TRUE) {
  # creates a folder for images
  dir.create("Statistics")
  setwd("Statistics")

  # network plot - not working
  # tiff
  # tiff("network_plot_0.7.tiff", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
  # networkplot(pseudodistmat, highlight_annotated = TRUE, show_labels = TRUE, exclude_singletons = TRUE, min_similarity = 0.7)
  # dev.off()
  # png
  # png("network_plot_0.7.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
  # networkplot(pseudodistmat, highlight_annotated = TRUE, show_labels = TRUE, exclude_singletons = TRUE, min_similarity = 0.7)
  # dev.off()

  # hierarchy plot
  # tiff
  tiff("hierarchy_plot.tiff", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
  par(mar = c(4.5, 4.2, 1, 0.5))
  HCplot(pseudodistmat, h = 0.7, cex = 0.5)
  dev.off()
  # png
  png("hierarchy_plot.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
  par(mar = c(4.5, 4.2, 1, 0.5))
  HCplot(pseudodistmat, h = 0.7, cex = 0.5)
  dev.off()

  # heatmap plot
  # tiff
  tiff("heatmap.tiff", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
  par(mar = c(4.5, 4.2, 1, 0.5))
  HCplot(pseudodistmat, type = "heatmap", cexRow = 0.4, cexCol = 0.4, margins = c(7, 7))
  dev.off()
  # png
  png("heatmap.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
  par(mar = c(4.5, 4.2, 1, 0.5))
  HCplot(pseudodistmat, type = "heatmap", cexRow = 0.4, cexCol = 0.4, margins = c(7, 7))
  dev.off()
  }
  
  # set to main folder
  setwd(myDir)

  # return results
  return(list(apslist, r))
}
