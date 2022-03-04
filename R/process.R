#' process
#'
#' This function computes the second step of workflow (peak picking, RT correction and peak grouping) and build images for visualization of data using functions from xcms package.
#' @export
#' @param myDir Path to the directory of work.
#' @param raw_data A 'MSnExp' from MSnbase package.
#' @param metadata A matrix or data.frame with metadata information about samples. Include, at least 'sample' and 'file' columns with name of sample and its path, respectively. More information can be added in new columns, such as 'group', 'class', 'biorep' and 'tecrep'.
#' @param colors A list with colors generated from "read_data()".
#' @param EIC Numeric. 1 = there are ions to monitor through the processing. 2 = there are none. Default to 2.
#' @param ions List with sublist mz = mz (numeric) of the monitored ion and rt = retention time of monitored ion (numeric). To the 'rt' will be added and subtracted 5 seconds. Default to null.
#' @return A 'xcmsSet' or a 'XCMSnExp' object with detected, grouped and filled peaks with retention time corrected.
#' @importFrom grDevices dev.off pdf png tiff
#' @importFrom graphics boxplot grid legend par text
#' @importFrom methods as new
#' @importFrom svDialogs dlgInput
#' @importFrom xcms MatchedFilterParam findChromPeaks refineChromPeaks FilterIntensityParam adjustRtime ObiwarpParam groupChromPeaks PeakDensityParam fillChromPeaks ChromPeakAreaParam plotChromPeakImage chromPeaks plotAdjustedRtime chromatogram fillPeaks
#' @importFrom utils choose.dir memory.limit menu read.csv select.list write.csv write.table
#' @examples
#' \dontrun{
#' read_data()
#' xdata4 <- process(raw_data, metadata, myDir, colors)
#' }

process <- function(raw_data, metadata, myDir, colors, EIC = 2, ions = NULL) {
  # create and set folder for images
  dir.create("peakProcessing_results")
  setwd("peakProcessing_results")

  # peak picking
  mfp <- MatchedFilterParam(fwhm = 5, binSize = 0.5, steps = 2, mzdiff = 0.5, snthresh = 2, max = 500)
  xdata <- findChromPeaks(raw_data, param = mfp)
  save(xdata, file = "xdata.RData")

  # apply intensity filter? how much?
  filt <- menu(c("Yes", "No"), graphics = TRUE, title = "Apply intensity filter?")
  if (filt == 1) {
    filter <- dlgInput("Intensity threshold ", "0")$res
    xdata <- refineChromPeaks(xdata, param = FilterIntensityParam(threshold = as.integer(filter), nValues = 1, value = "maxo"))
  }

  # retention time correction
  xdata2 <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
  save(xdata2, file = "xdata2.RData")

  # grouping peaks
  xdata3 <- groupChromPeaks(xdata2, param = PeakDensityParam(sampleGroups = xdata2$group, bw = 0.5, minSamples = 1, maxFeatures = 500, minFraction = 0.4))
  save(xdata3, file = "xdata3.RData")

  # fill missing peaks
  xdata4 <- fillChromPeaks(xdata3, param = ChromPeakAreaParam())

  # imagens dos dados em processamento e pós processamento (padrões, cromatogramas de íon extraído)


  if (exists("xdata4")) {
    # heatmap of identified peaks per region of chromatogram
    # tiff
    tiff("plotChromPeakImage.tiff", units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
    par(mar = c(5, 9, 4, 1) + .1)
    plotChromPeakImage(xdata4)
    dev.off()
    # png
    png("plotChromPeakImage.png", units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
    par(mar = c(5, 9, 4, 1) + .1)
    plotChromPeakImage(xdata4)
    dev.off()

    # boxplot of log2 intensities per sample
    ints <- split(log2(chromPeaks(xdata4)[, "into"]),
      f = chromPeaks(xdata4)[, "sample"]
    )
    names(ints) <- metadata$all2
    for (i in 1:length(colors)) {
      # tiff
      tiff(paste0(names(colors)[i], "_boxplotLog2Postprocessed.tiff"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
      par(mar = c(7, 5, 3, 1) + .1, cex.axis = 1)
      boxplot(ints,
        varwidth = TRUE, col = colors[[i]][[2]][colors[[i]][[1]]],
        ylab = expression(log[2] ~ intensity), main = "Peak intensities", las = 3, xaxt = "n"
      )
      grid(nx = NA, ny = NULL)
      text(seq_along(metadata$all2), par("usr")[3],
        labels = metadata$all2, srt = 45, adj = c(1.1, 1.1), xpd = TRUE, cex = 0.7
      )
      dev.off()
      # png
      png(paste0(names(colors)[i], "_boxplotLog2Postprocessed.png"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
      par(mar = c(7, 5, 3, 1) + .1, cex.axis = 1)
      boxplot(ints,
        varwidth = TRUE, col = colors[[i]][[2]][colors[[i]][[1]]],
        ylab = expression(log[2] ~ intensity), main = "Peak intensities", las = 3, xaxt = "n"
      )
      grid(nx = NA, ny = NULL)
      text(seq_along(metadata$all2), par("usr")[3],
        labels = metadata$all2, srt = 45, adj = c(1.1, 1.1), xpd = TRUE, cex = 0.7
      )
      dev.off()
    }

    # chromatogram postprocessed
    bpc_after <- chromatogram(xdata4, aggregationFun = "max", include = "none")
    for (i in 1:length(colors)) {
      # tiff
      tiff(paste0(names(colors)[i], "_postprocessedChromatogram.tiff"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
      par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
      plot(bpc_after, col = colors[[i]][[2]][colors[[i]][[1]]])
      legend("topright", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 1, bg = "transparent")
      plotAdjustedRtime(xdata4, col = colors[[i]][[2]][colors[[i]][[1]]])
      legend("bottomright", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 1, bg = "transparent")
      dev.off()
      # png
      png(paste0(names(colors)[i], "_postprocessedChromatogram.png"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
      par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
      plot(bpc_after, col = colors[[i]][[2]][colors[[i]][[1]]])
      legend("topright", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 1, bg = "transparent")
      plotAdjustedRtime(xdata4, col = colors[[i]][[2]][colors[[i]][[1]]])
      legend("bottomright", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 1, bg = "transparent")
      dev.off()
    }

    # cromatograma de íons extraído
    if (EIC == 1) {
      dir.create("Monitoring ions")
      setwd("Monitoring ions")

      for (ii in 1:length(ions)) {
        crom <- chromatogram(xdata4, rt = c(ions[[ii]][["rt"]] - 5, ions[[ii]][["rt"]] + 5), mz = ions[[ii]][["mz"]], include = "none")
        for (i in 1:length(colors)) {
          # tiff
          tiff(paste0(names(colors)[i], "_", ii, "_postPross_EIC.tiff"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
          plot(crom, col = colors[[i]][[2]][colors[[i]][[1]]])
          legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
          dev.off()
          # png
          png(paste0(names(colors)[i], "_", ii, "_postPross_EIC.png"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
          plot(crom, col = colors[[i]][[2]][colors[[i]][[1]]])
          legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
          dev.off()
        }
      }
    }
  } else {
    xdata4 <- as(xdata3, "xcmsSet")
    xdata4 <- fillPeaks(xdata4)
  }



  # set to main folder
  setwd(myDir)

  # return results
  return(xdata4)
}
