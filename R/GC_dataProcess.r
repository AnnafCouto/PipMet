#' GC_dataProcess
#'
#' Main function for GC data processing.
#' @keywords metadata
#' @export
#' @param myDir Path to working directory. Default to none. Only for example used.
#' @param sample_dir Path to sample directory. Default to none. Only for example used.
#' @param metadata Path to .csv file or an R data.frame object containing metadata. At least 'sample' and 'file' columns must be included. Default to none. Only for example used.
#' @param extensao Extension of mass spectrometry files to read. Only accepted '.mzML' and '.mzXML'. Default to none. Only for example used.
#' @param pictures Logical. If pictures should be plotted or not.
#' @param example Logical. If is example, pop-ups won't appear.
#' @return A list containing (1) the path of working folder, (2) the metadata table, (3) the annotated pseudospectra list, (4) a OnDiskMSnExp object, (5) a XCMSnExp or xcmsSet object, (6) a xsAnnotate object, (7) a list of colors used and (8) the normalized instensities quantification table
#' @importFrom methods as new
#' @importFrom svDialogs dlgInput dlg_message
#' @importFrom ddpcr quiet
#' @importFrom utils choose.dir menu read.csv write.csv write.table
#' @importFrom metaMS addRI write.msp
#' @importFrom BiocParallel register SerialParam SnowParam MulticoreParam
#' @importFrom tcltk tkmessageBox
#' @examples
#' \dontrun{
#' result <- GC_dataProcess(
#'   sample_dir = system.file("extdata", package = "PipMet"),
#'   metadata = system.file("extdata", "metadata.csv", package = "PipMet"),
#'   extensao = ".mzML",
#'   myDir = "~/",
#'   example = TRUE,
#'   pictures = TRUE
#' )
#' }
#'
GC_dataProcess <- function(myDir = NULL, sample_dir = NULL, metadata = NULL, extensao = c(".mzML", ".mzXML"), pictures = TRUE, example = FALSE) {

  # ask for monitoring ions infos - CHECAR SE FUNCIONA
  if (!example == TRUE & pictures == TRUE) {
    EIC <- menu(c("Yes", "No"), graphics = TRUE, title = "Would you like to monitor EICs?")
    ions <- list()
    if (EIC == 1) {
      okay <- 1
      ei <- 1
      while (okay == 1) {
        ions[[ei]] <- vector(mode = "list", length = 2)
        names(ions[[ei]]) <- c("mz", "rt")
        # names(ions)[ei] <- as.character(dlgInput(paste0('Name monitoring ion ', ei, ' :'), 'First')$res)
        ions[[ei]][["mz"]] <- as.integer(dlgInput(paste0("Mz of EIC ", ei, ":"), "0")$res)
        ions[[ei]][["rt"]] <- as.integer(dlgInput(paste0("Rt of EIC ", ei, " (automatically will be added +/- 5s to Rt):"), "0")$res)
        okay <- menu(c("Yes", "No"), graphics = TRUE, title = "Would you like to monitor another one?")
        ei <- ei + 1
      }
    }
  } else {
    EIC <- 2
  }

  # ask for parallelization mode
  if (!example == TRUE) {
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
  } else {
    register(SerialParam(), default = TRUE)
  }

  # ask informations and read files
  quiet(read <- read_data(EIC, ions, sample_dir = sample_dir, metadata = metadata, extensao = extensao, myDir = myDir, pictures = pictures, example = example))
  colors <- read[[1]]
  metadata <- read[[2]]
  raw_data <- read[[3]]
  myDir <- read[[4]]
  rm(read)

  # process samples
  quiet(xdata4 <- process(raw_data, metadata, myDir, colors, EIC, ions, pictures))

  # define spectra and create .msp files
  quiet(spectra <- getSpectra(xdata4, example, raw_data, colors))
  anIC <- spectra[[1]]
  result <- spectra[[3]]
  pslist <- spectra[[2]]

  # add retention index info if needed
  if (!example == TRUE) {
    ri <- menu(c("Yes", "No"), graphics = TRUE, title = 'Would you like to add retention index data to your spectra? For that, you must provide a .csv file with column "rt" and "RI" in you directory.')
    if (ri == 1) {
      RI <- read.csv("RI.csv")
      result <- addRI(result, RI)
      write.msp(result, "spectra.msp", newFile = TRUE)
      tkmessageBox(title = "Retention index", message = "The retention index for the spectra was calculated and added to the .msp file.", icon = "info", type = "ok")
      rm(RI)
    }
    rm(spectra, result, ri)
  }

  # update annotated spectra and plot images
  quiet(annot <- annot_images(pslist, myDir, pictures))
  apslist <- annot$apslist
  pre_anno <- annot$r

  # normalize, choose peaks and plot images
  quiet(n <- normalize_data(anIC, pslist, metadata, myDir, pre_anno, example))

  if (pictures == TRUE) {
    if (!example == TRUE) {
      # plot volcanos
      okay <- 1
      while (okay == 1) {
        quiet(x <- try(vol_lvl1(n, metadata, myDir)))
        okay <- menu(c("Repeat", "Next"), graphics = TRUE, title = "Plot volcano 1 level-comparison again?")
      }

      okay <- 1
      while (okay == 1) {
        quiet(x <- try(vol_lvl2(n, metadata, myDir, volDir)))
        okay <- menu(c("Repeat", "Next"), graphics = TRUE, title = "Plot volcano 2 level-comparison again?")
      }

      # plot PCA
      quiet(PCA_(n, metadata, myDir, colors))
    }

    # plot heatmaps
    quiet(heatmap(n, metadata, myDir, colors))
  }

  # save session
  save.image(paste0(Sys.Date(), ".R"))

  dlg_message("Processing done!")$res

  return(list(myDir = myDir, metadata = metadata, apslist = apslist, raw_data = raw_data, xdata4 = xdata4, anIC = anIC, colors = colors, quantification_table = n))
}
