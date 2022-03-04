#' GC_dataProcess
#'
#' Main function for GC data processing.
#' @keywords metadata
#' @export
#' @param myDir Path to working directory. Default to none. Only for example used.
#' @param sample_dir Path to sample directory. Default to none. Only for example used.
#' @param metadata Path to .csv file or an R data.frame object containing metadata. At least 'sample' and 'file' columns must be included. Default to none. Only for example used.
#' @param extensao Extension of mass spectrometry files to read. Only accepted '.mzML' and '.mzXML'. Default to none. Only for example used.
#' @return A list containing (1) the path of working folder, (2) the metadata table, (3) the annotated pseudospectra list, (4) a OnDiskMSnExp object, (5) a XCMSnExp or xcmsSet object, (6) a xsAnnotate object, (7) a list of colors used and (8) the normalized instensities matrix
#' @importFrom methods as new
#' @importFrom svDialogs dlgInput
#' @importFrom ddpcr quiet
#' @importFrom utils choose.dir menu read.csv write.csv write.table
#' @importFrom metaMS addRI write.msp
#' @importFrom BiocParallel register SerialParam SnowParam MulticoreParam
#' @importFrom tcltk tkmessageBox
#' @examples
#' \dontrun{
#' result <- GC_dataProcess()
#' }

GC_dataProcess <- function(myDir = NULL, sample_dir = NULL, metadata = NULL, extensao = c('.mzML', '.mzXML')) {

  # ask for monitoring ions infos - CHECAR SE FUNCIONA
  EIC <- menu(c("Yes", "No"), graphics = TRUE, title = "Would you like to monitor EICs?")
  ions <- list()
  if (EIC == 1) {
    okay <- 1
    ei <- 1
    while (okay == 1) {
      ions[[ei]] <- vector(mode = "list", length = 2)
      names(ions[[ei]]) <- c("mz", "rt")
      # names(ions)[ei] <- as.character(dlgInput(paste0('Name monitoring ion ', ei, ' :'), 'First')$res)
      ions[[ei]][["mz"]] <- as.integer(dlgInput(paste0("Mz of EIC ", ei, " :"), "0")$res)
      ions[[ei]][["rt"]] <- as.integer(dlgInput(paste0("Rt of EIC ", ei, " (automatically will be add +/- 5s to Rt):"), "0")$res)
      okay <- menu(c("Yes", "No"), graphics = TRUE, title = "Would you like to monitor another one?")
      ei <- ei + 1
    }
  }

  # ask for parallelization mode
  parallel <- menu(c("Serial Param (disable)", "Snow Param", "MultiCore Param"), graphics = TRUE, title = "Choose parallelization mode:")
  if (parallel == 1) {
    register(SerialParam(), default = TRUE)
  }
  if (parallel == 2) {
    register(SnowParam(), default = TRUE)
  }
  if (parallel == 3) {
    register(MulticoreParam(), default = TRUE)
  }

  # ask informations and read files
  quiet(read <- read_data(EIC, ions))
  colors <- read[[1]]
  metadata <- read[[2]]
  raw_data <- read[[3]]
  myDir <- read[[4]]
  rm(read)

  # process samples
  quiet(xdata4 <- process(raw_data, metadata, myDir, colors))

  # define spectra and create .msp files
  quiet(spectra <- getSpectra(xdata4))
  anIC <- spectra[[1]]
  pslist <- spectra[[2]]
  result <- spectra[[2]]

  # add retention index info if needed
  ri <- menu(c("Yes", "No"), graphics = TRUE, title = 'Would you like to add retention index data to your spectra? For that, you must provide a .csv file with column "rt" and "RI" in you directory.')
  if (ri == 1) {
    RI <- read.csv("RI.csv", sep = ";", dec = ",")
    result <- addRI(result, RI)
    write.msp(result, "spectra.msp", newFile = TRUE)
    tkmessageBox(title = "Retention index", message = "The retention index for the spectra was calculated and added to the .msp file.", icon = "info", type = "ok")
    rm(RI)
  }
  rm(spectra, result, ri)

  # update annotated spectra and plot images
  quiet(apslist <- annot_images(pslist, myDir))

  # normalize, choose peaks and plot images
  quiet(n <- normalize(anIC, pslist, metadata, myDir))

  # plot volcanos
  okay <- 1
  while (okay == 1) {
    quiet(mat <- vol_lvl1(n, metadata, myDir))
    volDir <- mat[[2]]
    mat <- mat[[1]]
    okay <- menu(c("Repeat", "Next"), graphics = TRUE, title = "Plot volcano 1 level-comparison again?")
  }

  okay <- 1
  while (okay == 1) {
    quiet(vol_lvl2(mat, n, metadata, myDir, volDir))
    okay <- menu(c("Repeat", "Next"), graphics = TRUE, title = "Plot volcano 2 level-comparison again?")
  }

  # plot PCA
  quiet(PCA_(mat, metadata, myDir, colors))

  # plot heatmaps
  quiet(heatmap(mat, n, metadata, myDir, colors))

  dlg_message("Processing done!")$res

  return (list (myDir, metadata, apslist, raw_data, xdata4, anIC, colors, n))
}
