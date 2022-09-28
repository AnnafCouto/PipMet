#' workData
#'
#' Main function for GC data processing. It read files in '.mzXML' and '.mzML' formats, apply algorithms for for peacking, grouping, retention time correction and filling missing peaks. The next step proposes spectra based on retention time window and correlation information. After annotation, the function performs quantification by normalizing data and apply t test and principal components analysis, plotting pictures for the evaluation of the steps performed. At the end, the user can choose to generate a internal library with the identified compounds to be uploaded into NIST MS Search software.
#' @keywords metadata in-house database preprocessing spectra spectrum peakpicking grouping RI rt
#' @export
#' @param myDir Path to working directory. Default to none.
#' @param sample_dir Path to sample directory. Default to none.
#' @param metadata Path to .csv file or an R data.frame object containing metadata. At least 'sample' and 'file' columns must be included. Default to none.
#' @param extension Extension of mass spectrometry files to read. Only accepted '.mzML' and '.mzXML'. Default to none.
#' @param pictures Logical. If pictures should be plotted or not.
#' @param example Logical. If example = TRUE, the metadata and other needed files will be loaded from package files.
#' @param filter Numeric. Intensity threshold for the peak detection. Default to NULL. When NULL, the user will be asked for a number. Set filter = 0 for no intensity filtering.
#' @param peakMonitor Logical. Are there peak to monitor throuhout the workflow? If NULL, the user will be asked. Default to NULL.
#' @param parallel Character. Sort of parallelization for code to perfom. Supported are "Serial Param", "Snow Param", "MultiCore Param". For more information, check the BiocParallel R package. Default to NULL. If parallel = NULL, the user will be asked.
#' @param pic_extension Character. Pictures format to generate. Supported = '.tiff', '.png'. Default to c('.tiff', '.png').
#' @param group Character. Name from 'metadata' column names to group the samples. Default to 'group'.
#' @param derivatization Character. Kind of derivatization the samples were prepared with. Supported are 'Trimethylsilyl' and 'None'. If NULL, the user will be asked. Default to 'NULL'.
#' @param cores Numeric. Number of cores to be used in Snow Param. Default to NULL. If NULL, the user will be asked. Set cores = 0 to Serial Param.
#' @param column_set Character. Polarity of column used for the chromatography: 'polar', 'non-polar'. If NULL, the user will be asked. Default to NULL.
#' @param prog Character. Configuration of temperature in data acquisition: "isothermal", "ramp", "custom". If NULL the user will be asked. Default to NULL.
#' @param RI Logical or path to retention index .csv file. Addition of retention index information to the spectra. If RI = TRUE, the user will be asked to provide a .csv file with 'rt' and 'RI' columns. If RI = path to the .csv files, the retention index will be calculated. Default to FALSE.
#' @param ion_mode Character. Ion mode acquisition 'positive' or 'negative'. If NULL, the user will be asked. Default to NULL.
#' @param plot_eic Logical. Plot the EIC of each of the 6 most intense m/z in the spectra. Default to NULL.
#' @param lib_build Logical. For lib_build == TRUE, the identified compounds will be gathered into a new internal library. Default to NULL.
#' @param replicate Character or Logical. If FALSE, there is no replicate informations in the metadata table. Otherwise, inform the name in the metadata column containing replicate information. If NULL, the user will be asked. Default to NULL.
#' @param removeCompounds Logical. If TRUE, the user may choose from a pop-up identified compounds to remove from the quantification table. If NULL, the user will be asked. Default to 'NULL'.
#' @param mergeCompounds Logical. If TRUE, the user may choose from a pop-up what compounds identified should be representated as one with intensities summed. Exemple: derivatizations derivatives. If NULL, the user will be asked. Default to 'NULL'.
#' @param info Logical. If TRUE, the user cann add more information to the creatin of the database, such as CAS number, PubMed and ChemSpider ID and InChiKey, among others. If NULL, the user will be asked. Default to NULL. Only required if lib_build = TRUE or NULL.
#' @param pre_anno Path to pre_anno.csv annotation file or pre_anno table. If NULL, the pre_anno.csv in the folder will be read. Default to NULL.
#' @param Ri Character. Retention index calculation method. Supported are 'lee', 'linear', 'kovats' and 'alcane'. If NULL, the user will be asked. Default ot NULL. Only if lib_build = TRUE or NULL.
#' @param min_peaks Numeric. Minimal number of peaks a spectrum must have to be considered a viable spectrum. Default to 5.
#' @param sample_names Character. Name of metadata column to name the samples. Default to NULL. If NULL, the user will be asked.
#' @return A list containing (1) the path of working folder, (2) the metadata table, (3) the annotated pseudospectra list, (4) a OnDiskMSnExp object, (5) a XCMSnExp or xcmsSet object, (6) a xsAnnotate object, (7) a list of colors used and (8) the normalized instensities quantification table
#' @importFrom methods as new
#' @importFrom svDialogs dlgInput dlg_message dlg_list
#' @importFrom ddpcr quiet
#' @importFrom utils choose.dir menu read.csv write.csv write.table
#' @importFrom metaMS addRI write.msp
#' @importFrom parallel detectCores
#' @importFrom BiocParallel register SerialParam SnowParam MulticoreParam
#' @importFrom tcltk tkmessageBox
#' @importFrom fritools is_path
#' @examples
#' \donttest{
#' \dontrun{
#' result <- workData(
#'   extension = ".mzXML",
#'   myDir = "~/",
#'   example = TRUE,
#'   pictures = TRUE 
#' )
#' }
#' }
workData <- function(myDir = NULL, sample_dir = NULL, metadata = NULL, extension = NULL, pictures = TRUE, example = FALSE, filter = NULL, peakMonitor = NULL, pic_extension = c('.tiff', '.png'), parallel = NULL, group = NULL, derivatization = NULL, cores = 1, column_set = NULL, prog = NULL, ion_mode = NULL, plot_eic = NULL, lib_build = NULL, RI = NULL, replicate = NULL, mergeCompounds = NULL, removeCompounds = NULL, info = NULL, Ri = NULL, pre_anno = NULL, min_peaks = 5, sample_names = NULL) {

  # ask for monitoring ions infos - CHECAR SE FUNCIONA
  if (pictures == TRUE) {
    if (is.null (peakMonitor)) {
      peakMonitor <- dlg_list(c("Yes", "No"), multiple = FALSE, title = "Monitor EICs?")$res
    } 
    ions <- list()
    if (peakMonitor == 'Yes' | peakMonitor == TRUE) {
      okay <- 1
      ei <- 1
      while (okay == 1) {
        ions[[ei]] <- vector(mode = "list", length = 2)
        names(ions[[ei]]) <- c("mz", "rt")
        ions[[ei]][["mz"]] <- as.integer(dlgInput(paste0("Mz of EIC ", ei, ":"), "0")$res)
        ions[[ei]][["rt"]] <- as.integer(dlgInput(paste0("Rt of EIC ", ei, " (automatically will be added +/- 5s to Rt):"), "0")$res)
        okay <- menu(c("Yes", "No"), graphics = TRUE, title = "Monitor another one?")
        ei <- ei + 1
      }
    }
  }

  # ask for parallelization mode
    if (is.null (parallel)) {
      parallel <- dlg_list(c("Serial Param", "Snow Param", "MultiCore Param"), multiple = FALSE, title = "Parallelization mode:")$res
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

  # ask informations and read files
  message ('Reading files...')
  quiet(read <- read_data(peakMonitor = peakMonitor, ions = ions, sample_dir = sample_dir, metadata = metadata, extension = extension, myDir = myDir, pictures = pictures, example = example, pic_extension = pic_extension, sample_names = sample_names))
  colors <- read[[1]]
  metadata <- read[[2]]
  raw_data <- read[[3]]
  myDir <- read[[4]]

  # process samples
  message ('Pre-processing files...')
  quiet(processed <- process(raw_data = raw_data, metadata = metadata, myDir = myDir, colors = colors, peakMonitor = peakMonitor, ions = ions, pictures = pictures, filter = filter, pic_extension = pic_extension, group = group))
  xdata4 <- processed$xdata4
  group <- processed$group

  # define spectra and create .msp files
  message (paste0('Grouping peaks into spectra...'))
  quiet(spectra <- getSpectra(xdata4, raw_data, min_peaks, colors, column_set = column_set, prog = prog, ion_mode = ion_mode, plot_eic = plot_eic))
  anIC <- spectra[[1]]
  result <- spectra[[3]]
  pslist <- spectra[[2]]
  ion_mode <- spectra[[4]]

  dlg_message("Annotatation step: The files 'pre_anno.csv' and 'spectra.msp' were created in you directory. Upload the file 'spectra.msp' in NIST MS Search and annotate the spectra in the file 'pre_anno', in the column 'Annotation', according to the spectra 'id'. After, press 'ok'.")$res

  # add retention index info if needed
    if (is.null(RI)) {RI <- dlg_message("Add retention index information?", "yesno")$res}
    if (RI == TRUE | RI == 'yes') {
      RI <- read.csv(choose.files())
    }
    if (is_path(RI)) {
      RI <- read.csv(RI)
      message ('Calculating retention index...')
      result <- addRI(result, RI)
      write.msp(result, "spectra_RI.msp", newFile = TRUE)
      dlg_message("The retention index for the spectra was calculated and added to the .msp file.")
      #rm(spectra, result)
    }

  # update annotated spectra and plot images
  message (paste0('Annotating the spectra ...'))
  quiet(annot <- annot_images(pslist, myDir, pictures, pic_extension = pic_extension, pre_anno = pre_anno))
  apslist <- annot$apslist
  pre_anno <- annot$r

  # normalize, choose peaks and plot images, only if more than one sample

  message (paste0('Normalizing data...'))
  quiet(n <- normalize_data(anIC, pslist, metadata, myDir, pre_anno, pic_extension = pic_extension, derivatization, mergeCompounds, removeCompounds, group))
  
  if (nrow(metadata) > 1 & nrow(pre_anno) > 1) {
    if (is.null(replicate)) {
      if (!'tec_rep' %in% colnames(metadata)) {
        replicate <- dlg_list(c(colnames(metadata), 'No information'), multiple = FALSE, title = "Replicate column:")$res}
        } else {replicate <- 'tec_rep'}
    
    message (paste0('Statistics pictures...'))

    if (pictures == TRUE & nrow(metadata) > 1) {
      pic <- dlg_list(c('Volcano - level 1', 'Volcano - level 2', 'PCA - All spectra', 'PCA - Identified spectra only', 'Heatmaps', 'All', 'None'), multiple = TRUE, title = "Statistics pictures:")$res
      # plot volcanos
      okay <- 1
      if ('Volcano - level 1' %in% pic | 'All' %in% pic) { # Volcano level 2 can only be made if volcano - level 1 is too
        while (okay == 1) {
          quiet(volDir <- try(vol_lvl1(n, metadata, myDir, pic_extension = pic_extension)))
          okay <- menu(c("Repeat", "Next", ""), graphics = TRUE, title = "Repeat?")
        }
        okay <- 1
        if ('Volcano - level 2' %in% pic | 'All' %in% pic) {
          while (okay == 1) {
            quiet(x <- try(vol_lvl2(n, metadata, myDir, volDir, pic_extension = pic_extension)))
            okay <- menu(c("Repeat", "No"), graphics = TRUE, title = "Repeat Volcanos?")
          }
        }
      }

      # plot PCA
      if ('PCA - All spectra' %in% pic | 'All' %in% pic) {quiet(PCA_general(n, metadata, myDir, colors, pic_extension = pic_extension, example = example))}
      if ('PCA - Identified spectra only' %in% pic | 'All' %in% pic) {quiet(PCA_identified(n, metadata, myDir, colors, pic_extension = pic_extension, example = example))}

      # plot heatmaps
      if ('Heatmaps' %in% pic | 'All' %in% pic) {quiet(heatmap(n, metadata, myDir, colors, pic_extension = pic_extension, replicate = replicate, sample_names = NULL))}
    }
  }

  dlg_message("Processing done!")$res

  # check create_database first
  if (is.null(lib_build)) {lib_build <- dlg_message("Would you like to create a in-house database with the identified spectra?", type = "yesno")$res}
  if (lib_build == TRUE | lib_build == 'yes') {
    create_database(apslist, column_set = column_set, prog = prog, ion_mode = ion_mode, info = info, Ri = Ri)
  }

  return(list(myDir = myDir, metadata = metadata, apslist = apslist, raw_data = raw_data, xdata4 = xdata4, anIC = anIC, colors = colors, quantification_table = n))
}
