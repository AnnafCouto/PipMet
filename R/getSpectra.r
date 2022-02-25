#' Spectra definition
#'
#' This function computes the first subsetp of third step of the workflow. It defines the spectra to proper annotate in NIST MS Search Software.
#' @export
#' @param xdata4 A 'xcmsSet' or 'XCMSnExp' object.
#' @importFrom grDevices dev.off pdf png tiff
#' @importFrom methods as
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr
#' @importFrom CluMSID writeFeaturelist
#' @importFrom metaMS construct.msp write.msp
#' @importFrom utils memory.limit
#' @importFrom tcltk tkmessageBox
#' @examples
#' read_data()
#' process()
#' spectra_def()
getSpectra <- function(xdata4) {
  if (class(xdata4)[1] == "XCMSnExp") {
    xset <- as(xdata4, "xcmsSet")
  } else {
    xset <- xdata4
  }
  polarity <- menu(c("Positive", "Negative"), graphics = TRUE, title = "Polarity:")
  if (polarity == 1) {
    an <- xsAnnotate(xset, polarity = "positive")
  } else {
    an <- xsAnnotate(xset, polarity = "negative")
  }
  anF <- groupFWHM(an, perfwhm = 1)
  rm(an, xset)
  memory.limit(100000)
  anIC <- groupCorr(anF, calcIso = FALSE)
  rm(anF)

  # extract spectra with minimum of 5 peaks
  pslist <- exctractSpectra(anIC, min_peaks = 5)

  # create 'pre_anno.csv' where the user annotates the spectra
  # contains id and retention time of each spectra
  writeFeaturelist(pslist)

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
    result[[i]]$Formula <- "Unknown"
    result[[i]]$monoMW <- "Unknown"
    result[[i]]$CAS <- "Unknown"
    result[[i]]$ChemSpiderID <- "Unknown"
    result[[i]]$Class <- "Unknown"
    result[[i]]$Date <- Sys.Date()
  }
  write.msp(result, "spectra.msp", newFile = TRUE)
  tkmessageBox(title = "Annotatation step", message = "The files 'pre_anno.csv' and 'spectra.msp' were created in you directory. Upload the file 'spectra.msp' in NIST MS Search and annotate the spectra in the file 'pre_anno', in the column 'Annotation', according to the spectra 'id'. After, press 'ok'.", icon = "info", type = "ok")
  return(list(anIC, pslist, result))
}
