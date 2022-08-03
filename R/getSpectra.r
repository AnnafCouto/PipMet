#' getSpectra
#'
#' This function computes the first subsetp of third step of the workflow. It defines the spectra to proper annotate in NIST MS Search Software.
#' @export
#' @param xdata4 A 'xcmsSet' or 'XCMSnExp' object.
#' @param example Logical. If is example, pop-ups won't appear. Default to FALSE.
#' @param raw_data A 'XCMSnExp' object.
#' @param colors A list with colors generated from "read_data()".
#' @return A list with 'xsAnnotate' object with peaks grouped by retention time and correlation peaks information, a 'pseudospectrum' object, a spectra list in .msp format and the ion mode of data acquisition ('negative' or 'positive').
#' @importFrom grDevices dev.off pdf png tiff
#' @importFrom methods as
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr
#' @importFrom CluMSID writeFeaturelist specplot
#' @importFrom metaMS construct.msp write.msp
#' @importFrom utils memory.limit
#' @importFrom svDialogs dlg_message dlg_input dlg_list
#' @examples
#' \donttest{
#' \dontrun{
#' load(system.file("extdata", "xdata4.RData", package = "PipMet"))
#' spectra <- getSpectra(xdata4, example = TRUE)
#' }
#' }

getSpectra <- function(xdata4, example = FALSE, raw_data, colors) {
  if (class(xdata4)[1] == "XCMSnExp") {
    xset <- as(xdata4, "xcmsSet")
  } else {
    xset <- xdata4
  }
  if (example == TRUE) {
    polarity <- "positive"
  } else {
    polarity <- dlg_list(c("positive", "negative"), multiple = FALSE, title = "Polarity:")$res # ask user for the polarity
  }
  an <- xsAnnotate(xset, polarity = polarity)
  anF <- groupFWHM(an, perfwhm = 1)
  rm(an, xset)
  memory.limit(100000)
  anIC <- groupCorr(anF, calcIso = FALSE)
  rm(anF)

  # extract spectra with minimum of 5 peaks
  pslist <- extractSpectra(anIC)

  # create 'pre_anno.csv' where the user annotates the spectra
  # contains id and retention time of each spectra
  writeFeaturelist(pslist)

  # get information about data acquisition
  column <- dlg_list(c("polar", "non-polar"), multiple = FALSE)$res
  prog <- dlg_list(c("isothermal", "ramp", "custom"), multiple = FALSE)$res
  Instrument_type <- dlg_input("Type of instrument of acquisition", "GC-EI-Q")$res

  # creates a .msp file for spectra
  spectra <- list()
  for (i in 1:length(pslist)) {
    x <- data.frame(cbind(pslist[[i]]@spectrum[, 1], (pslist[[i]]@spectrum[, 2] / max(pslist[[i]]@spectrum[, 2])))) # standardazing intensities by dividing the intensity of each peak from a spectrum by the maximum intensity of that spectra.
    colnames(x) <- c("mz", "into")
    spectra[[i]] <- x
  }
  result <- construct.msp(spectra, extra.info = NULL)
  for (i in 1:length(result)) {
    result[[i]]$id <- pslist[[i]]@id
    result[[i]]$rt <- pslist[[i]]@rt
    result[[i]]$Name <- paste0("Unknown ", pslist[[i]]@id)
    result[[i]]$Date <- as.character(Sys.Date())
    result[[i]]$Instrument_type <- Instrument_type
    result[[i]]$Comments <- paste0("Column class: ", paste0("Standard ", column), "; ", "ProgramType: ", prog)
    result[[i]]$Ion_mode <- polarity
  }
  write.msp(result, "spectra.msp", newFile = TRUE)

  # generate a pdf file with the spectrum, its chromatogram and the XIC of the 6 most intense m/z individually,
  # for every one of the generated spectrum
  if (!example == TRUE) {
    y <- dlg_message("Plot spectra, with EIC of each of the 6 most intense m/z?", "yesno")$res
    if (y == "yes") {
      z <- as.numeric(menu(colnames(metadata), graphics = TRUE, title = "Choose conditions (from metadata table) to group samples: "))
      z <- colnames(metadata)[z]
      rt <- list()
      for (i in 1:length(pslist)) {
        rt[[i]] <- as.numeric(pslist[[i]]@rt)
      }
      pdf("EIC_XIC.pdf")
      for (i in 1:length(rt)) {
        par(mfrow = c(2, 1))
        specplot(pslist[[i]])
        x <- paste0("Unknown ", pslist[[i]]@id)
        crom <- chromatogram(raw_data, rt = c(rt[[i]] - 8, rt[[i]] + 8))
        plot(crom, col = colors[[z]][[2]][colors[[z]][[1]]], main = paste0("Pre-processing - ", x))
        legend("right", legend = names(colors[[z]][[2]]), col = colors[[z]][[2]], fill = colors[[z]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
        # plot EIC from the 6 most intense ion-fragm
        par(mfrow = c(2, 3))
        f <- pslist[[i]]@spectrum[order(pslist[[i]]@spectrum[, 2], decreasing = TRUE)]
        sort(pslist[[i]]@spectrum[, 1], decreasing = TRUE)
        for (ii in 1:6) {
          crom <- chromatogram(raw_data, rt = c(rt[[i]] - 8, rt[[i]] + 8), mz = c(as.numeric(f[[ii]]) - 0.6, as.numeric(f[[ii]]) + 0.6))
          plot(crom, col = colors[[z]][[2]][colors[[z]][[1]]], main = paste0("EIC - mz", f[[ii]]))
          legend("right", legend = names(colors[[z]][[2]]), col = colors[[z]][[2]], fill = colors[[z]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
        }
      }
      dev.off()
    }
  }

  # done (with this step)
  dlg_message("Annotatation step: The files 'pre_anno.csv' and 'spectra.msp' were created in you directory. Upload the file 'spectra.msp' in NIST MS Search and annotate the spectra in the file 'pre_anno', in the column 'Annotation', according to the spectra 'id'. After, press 'ok'.")$res
  return(list(anIC = anIC, pslist = pslist, result = result, polarity = polarity))
}
