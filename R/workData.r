#' workData
#'
#' Main function for GC data processing.
#' @keywords metadata
#' @export
#' @param myDir Path to working directory. Default to none.
#' @param sample_dir Path to sample directory. Default to none.
#' @param metadata Path to .csv file or an R data.frame object containing metadata. At least 'sample' and 'file' columns must be included. Default to none.
#' @param extension Extension of mass spectrometry files to read. Only accepted '.mzML' and '.mzXML'. Default to none.
#' @param pictures Logical. If pictures should be plotted or not.
#' @param example Logical. If is example, pop-ups won't appear.
#' @param filter Numeric. Intensity threshold for the peak detection. Default to NULL. When NULL, the user will be asked for a number. Set filter = 0 for no intensity filtering.
#' @param peakMonitor Logical. Are there peak to monitor throuhout the workflow? Default to FALSE.
#' @param parallel Character. Sort of parallelization for code to perfom. Supported are "Serial Param", "Snow Param", "MultiCore Param". For more information, check the BiocParallel R package. If parallel = NULL, the user will be asked.
#' @param pic_extension Character. Pictures format to generate. Supported = '.tiff', '.png'. Default to c('.tiff', '.png').
#' @param group Character. Name from 'metadata' column names to group the samples. Default to 'group'.
#' @param derivatization Character. Kind of derivatization the samples were prepared with. Supported are 'Trimethylsilyl' and 'None'. If NULL, the user will be asked. Default to 'NULL'.
#' @param cores Numeric. Number of cores to be used in Snow Param. Default to NULL. If NULL, the user will be asked. Set cores = 0 to Serial Param.
#' @return A list containing (1) the path of working folder, (2) the metadata table, (3) the annotated pseudospectra list, (4) a OnDiskMSnExp object, (5) a XCMSnExp or xcmsSet object, (6) a xsAnnotate object, (7) a list of colors used and (8) the normalized instensities quantification table
#' @importFrom methods as new
#' @importFrom svDialogs dlgInput dlg_message
#' @importFrom ddpcr quiet
#' @importFrom utils choose.dir menu read.csv write.csv write.table
#' @importFrom metaMS addRI write.msp
#' @importFrom parallel detectCores
#' @importFrom BiocParallel register SerialParam SnowParam MulticoreParam
#' @importFrom tcltk tkmessageBox
#' @examples
#' \donttest{
#' \dontrun{
#' result <- workData(
#'   sample_dir = system.file("extdata", package = "PipMet"),
#'   metadata = system.file("extdata", "metadata.csv", package = "PipMet"),
#'   extension = ".mzXML",
#'   myDir = "~/",
#'   example = TRUE,
#'   pictures = TRUE
#' )
#' }
#' }
workData <- function(myDir = NULL, sample_dir = NULL, metadata = NULL, extension = NULL, pictures = TRUE, example = FALSE, filter = NULL, peakMonitor = NULL, pic_extension = c('.tiff', '.png'), parallel = NULL, group = 'group', derivatization = NULL, cores = 1) {

  # ask for monitoring ions infos - CHECAR SE FUNCIONA
  if (!example == TRUE & pictures == TRUE) {
    if (is.null (peakMonitor)) {
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
    }
  } #else {
    #EIC <- 2
  #}

  # ask for parallelization mode
  if (!example == TRUE) {
    if (is.null (parallel)) {
      parallel <- dlg_list(c("Serial Param", "Snow Param", "MultiCore Param"), multiple = FALSE, title = "Choose parallelization mode:")$res
    }
    #parallel <- menu(c("Serial Param (disable)", "Snow Param", "MultiCore Param"), graphics = TRUE, title = "Choose parallelization mode:")
    if (parallel == 'Serial Param') {
      register(SerialParam(), default = TRUE)
    }
    if (parallel == 'Snow Param') {
      if (is.null(cores)) {cores <- dlgInput(paste0("Number of cores used (you have ", detectCores()," cores available):"),  detectCores()-2)$res}
      register(SnowParam(workers = as.numeric(cores)), default = TRUE)
    }
    if (parallel == 'MultiCore Param') {
      register(MulticoreParam(), default = TRUE)
    }
  } else {register(SerialParam(), default = TRUE)}

  # ask informations and read files
  quiet(read <- read_data(peakMonitor = peakMonitor, ions = ions, sample_dir = sample_dir, metadata = metadata, extension = extension, myDir = myDir, pictures = pictures, example = example, pic_extension = pic_extension))
  colors <- read[[1]]
  metadata <- read[[2]]
  raw_data <- read[[3]]
  myDir <- read[[4]]
  rm(read)

  # process samples
  quiet(xdata4 <- process(raw_data, metadata, myDir, colors, peakMonitor, ions, pictures, filter = filter, pic_extension = pic_extension, group = group))

  # define spectra and create .msp files
  message (paste0('Grouping peaks into spectra...'))
  quiet(spectra <- getSpectra(xdata4, example, raw_data, colors))
  anIC <- spectra[[1]]
  result <- spectra[[3]]
  pslist <- spectra[[2]]
  polarity <- spectra[[4]]

  # add retention index info if needed
  if (!example == TRUE) {
    ri <- menu(c("Yes", "No"), graphics = TRUE, title = 'Would you like to add retention index data to your spectra? For that, you must provide a .csv file with column "rt" and "RI" in you directory.')
    if (ri == 1) {
      RI <- read.csv(choose.files())
      message (paste0('Calculating retention index...'))
      result <- addRI(result, RI)
      write.msp(result, "spectra.msp", newFile = TRUE)
      dlg_message("The retention index for the spectra was calculated and added to the .msp file.")$res
    }
    rm(spectra, result, ri)
  }

  # update annotated spectra and plot images
  message (paste0('Annotating the spectra ...'))
  quiet(annot <- annot_images(pslist, myDir, pictures, pic_extension = pic_extension))
  apslist <- annot$apslist
  pre_anno <- annot$r

  # normalize, choose peaks and plot images
  message (paste0('Normalizing data...'))
  quiet(n <- normalize_data(anIC, pslist, metadata, myDir, pre_anno, example, pic_extension = pic_extension, derivatization))

  message (paste0('Statistics pictures...'))
  if (pictures == TRUE & nrow(metadata) > 1) {
    if (!example == TRUE) {
      # plot volcanos
      okay <- 1
      while (okay == 1) {
        quiet(volDir <- try(vol_lvl1(n, metadata, myDir, pic_extension = pic_extension)))
        okay <- menu(c("Repeat", "Next"), graphics = TRUE, title = "Plot volcano 1 level-comparison again?")
      }

      okay <- 1
      while (okay == 1) {
        quiet(x <- try(vol_lvl2(n, metadata, myDir, volDir, pic_extension = pic_extension)))
        okay <- menu(c("Repeat", "Next"), graphics = TRUE, title = "Plot volcano 2 level-comparison again?")
      }
    }
    # plot PCA
    quiet(PCA_(n, metadata, myDir, colors, pic_extension = pic_extension, example))

    # plot heatmaps
    quiet(heatmap(n, metadata, myDir, colors, pic_extension = pic_extension))
  }

  # save session
  save.image(paste0(Sys.Date(), ".RData"))

  dlg_message("Processing done!")$res

  if (dlg_message("Would you like to create a in-house database with the identified spectra?", type = "yesno")$res == "yes") {
    create_database(apslist, polarity)
  }

  return(list(myDir = myDir, metadata = metadata, apslist = apslist, raw_data = raw_data, xdata4 = xdata4, anIC = anIC, colors = colors, quantification_table = n))
}
