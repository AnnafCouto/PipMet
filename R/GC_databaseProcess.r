#' GC_databaseProcess
#'
#' Main function for GC-MS standards processing for internal library development.
#' @keywords standards, library, internal
#' @export
#' @param myDir Path to working directory. Default to none.
#' @param sample_dir Path to sample directory. Default to none.
#' @param metadata Path to .csv file containing metadata. "compound", "formula", "exact.mass", "rt", "file", "CAS", "ChemSpider", "class", "RI", "InChIKey" columns must be included. Default to none.
#' @param extension Extension of mass spectrometry files to read. Only accepted '.mzML' and '.mzXML'. Default to none.
#' @param example Logical. If is example, pop-ups won't appear.
#' @importFrom methods as
#' @importFrom ddpcr quiet
#' @importFrom utils choose.dir menu read.csv write.csv
#' @importFrom metaMS addRI write.msp construct.msp
#' @importFrom MSnbase readMSData
#' @importFrom ProtGenerics filterRt
#' @importFrom tools file_path_as_absolute
#' @importFrom xcms MatchedFilterParam findChromPeaks
#' @importFrom BiocParallel register SerialParam
#' @importFrom webchem nist_ri as.cas
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr
#' @importFrom tcltk tkmessageBox
#' @importFrom ddpcr quiet
#' @importFrom svDialogs dlgInput dlg_message
#' @examples
#' \dontrun{
#' GC_databaseProcess(
#'   sample_dir = system.file("extdata", package = "PipMet"),
#'   lib_metadata = system.file("extdata", "lib_metadata.csv", package = "PipMet"),
#'   extension = ".mzXML",
#'   myDir = "~/",
#'   example = TRUE
#' )
#' }
#'
GC_databaseProcess <- function(myDir = NULL, sample_dir = NULL, lib_metadata = NULL, extension = NULL, example = FALSE) {

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

  # get directory for samples and project and set library metadata
  if (!example == TRUE) {
    # samples directory
    if (is.null(sample_dir) | missing(sample_dir)) {
      sample_dir <- choose.dir(default = getwd(), caption = "Please, select the Samples directory, should be C:/Users/_/Samples")
    }
    setwd("~/")
    # project directory
    if (is.null(myDir) | missing(myDir)) {
      myDir <- dlgInput("Name your library", Sys.info()["user"])$res
      dir.create(myDir, showWarnings = FALSE)
    }
    setwd(myDir)
    myDir <- getwd()
    # files extension
    if (is.null(extension)) {
      extension <- dlg_list(c(".mzML", ".mzXML"), multiple = TRUE, title = "Files extension:")$res
    }
    # library metadata set up
    if (is.null(lib_metadata)) {
      metd <- dlg_message("Metadata table already exists?", "yesno")$res
      if (metd == "yes") {
        lib_metadata <- read.csv(choose.files(), na.string = c("NA", ""), colClasses = "character", sep = ",")
        if (!"file" %in% colnames(lib_metadata)) {
          lib_metadata$file <- paste0(sample_dir, "/", lib_metadata$sample, extension)
        }
      } else {
        files <- list.files(sample_dir, full.names = TRUE, pattern = extension, recursive = TRUE)
        lib_metadata <- matrix(nrow = length(files), ncol = 10)
        colnames(lib_metadata) <- c("compound", "formula", "exact.mass", "rt", "file", "CAS", "ChemSpider", "class", "RI", "InChIKey")
        lib_metadata[, "file"] <- files
        lib_metadata[, "Compound"] <- sub(basename(files), pattern = extension, replacement = "", fixed = TRUE)
        write.csv(lib_metadata, "lib_metadata.csv", row.names = FALSE)
        dlg_message("A file 'lib_metadata.csv' was created in your directory. Fill the sheet before continuing. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.", "yesno")
        while (file.exists("lib_metadata.csv") == FALSE) {
          write.csv(lib_metadata, "lib_metadata.csv", row.names = FALSE)
          dlg_message("A file 'lib_metadata.csv' was created in you directory. Fill the sheet before continuing. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.", type = "ok")
          lib_metadata <- read.csv("lib_metadata.csv", na.string = c("NA", ""), colClasses = "character", sep = ",", dec = ".")
        }
        lib_metadata <- read.csv("lib_metadata.csv", na.string = c("NA", ""), colClasses = "character", sep = ",", dec = ".")
      }
    }
    # if example = TRUE
  } else {
    sample_dir <- system.file("extdata", package = "PipMet")
    lib_metadata <- read.csv(system.file("extdata", "lib_metadata.csv", package = "PipMet"), na.string = c("NA", ""), colClasses = "character", sep = ",", dec = ".")
    extension <- ".mzXML"
    setwd(sample_dir)
    for (i in 1:nrow(lib_metadata)) {lib_metadata$file[i] <- file_path_as_absolute(lib_metadata$file[i])}
    myDir <- "~/LibraryBuilding_example"
    setwd('~/')
    dir.create(myDir)
    setwd(myDir)
  }

  # remove empty lines/columns from lib_metadata
  lib_metadata <- Filter(function(x) !all(is.na(x)), lib_metadata) # remove empty columns
  lib_metadata <- lib_metadata[!apply(is.na(lib_metadata), 1, all), ] # remove empty rows

  # begin processing
  # get polarity of data acquisiton
  if (example == TRUE) {polarity <- 'positive'} else {polarity <- dlg_list(c("positive", "negative"), multiple = FALSE, title = "Polarity:")$res}
  raw_data <- list()
  for (i in 1:nrow(lib_metadata)) {
    # read - OK
    quiet(raw_data[[i]]<- readMSData(lib_metadata$file[i], mode = 'onDisk'))
    # filter RT
    quiet(raw_data[[i]] <- filterRt(raw_data[[i]], c(as.numeric(lib_metadata[i, "rt"]) - 10, as.numeric(lib_metadata[i, "rt"]) + 10)))
  }
  # peak picking
  quiet(xdata <- lapply(raw_data, FUN = findChromPeaks, param = MatchedFilterParam(fwhm = 5, binSize = 0.5, steps = 2, mzdiff = 0.5, snthresh = 2, max = 500)))
  # spectra defition
  quiet(xset <- lapply(xdata, FUN = as, Class = "xcmsSet"))
  # create a xsAnnotate object
  an <- lapply(xset, FUN = xsAnnotate, polarity = polarity)
  # group by rt
  anF <- lapply(an, FUN = groupFWHM, perfwhm = 1)
  # group by correlation information
  anIC <- lapply(anF, FUN = groupCorr, calcIso = FALSE)
  # extract to spectra list
  pslist <- lapply(anIC, FUN = extractSpectra, min_peaks = 5)

  names(pslist) <- names(anIC) <- names(an) <- names(anF) <- names(xset) <- names(xdata) <- names(raw_data) <- lib_metadata[, "compound"]

  # get information on retention index, column and temperature program
  Ri <- dlg_list(c("kovats", "linear", "alkane", "lee"), multiple = FALSE, title = "Retention time index")$res
  column <- dlg_list(c("polar", "non-polar"), multiple = FALSE)$res
  prog <- dlg_list(c("isothermal", "ramp", "custom"), multiple = FALSE)$res
  Instrument_type <- dlg_input("Type of instrument of acquisition", "GC-EI-Q")$res

  # retrieve RI from CAS numbers (from NIST database), if RI is not present
  lib_metadata[, "RI"] <- NA
  for (i in 1:nrow(lib_metadata)) {
    lib_metadata[i, "CAS"] <- gsub(pattern = "-", replacement = "", lib_metadata[i, "CAS"]) # remove hifens from CAS
    if (is.na(lib_metadata[i, "RI"])) {
      lib_metadata[i, "RI_"] <- mean(nist_ri(as.cas(lib_metadata[i, "CAS"]), from = "cas", type = Ri, polarity = column, temp_prog = prog)$RI)
    }
  }

  # create .msp structure and files
  dir.create("spectra")
  setwd("spectra")
  for (i in 1:length(pslist)) {
    spectra <- list()
    for (ii in 1:length(pslist[[i]])) {
      x <- cbind(pslist[[i]][[ii]]@spectrum[, 1], (pslist[[i]][[ii]]@spectrum[, 2] / max(pslist[[i]][[ii]]@spectrum[, 2]))) # padronizo dividindo todas as intensidades de um mesmo espectro pela maior intensidade no mesmo (fica tipo 1 e 0,X ou seja, porcentagens)
      x <- data.frame(x)
      colnames(x) <- c("mz", "into")
      spectra[[ii]] <- x
    }
    result <- construct.msp(spectra, extra.info = NULL)
    for (ii in 1:length(result)) {
      result[[ii]]$Name <- paste0(pslist[[i]][[ii]]@id, " - Candidate to ", names(pslist)[i])
      result[[ii]]$id <- pslist[[i]][[ii]]@id
      result[[ii]]$rt <- pslist[[i]][[ii]]@rt
    }
    metaMS::write.msp(result, paste0(names(pslist)[i], "_", Sys.Date(), ".msp"), newFile = TRUE)
  }

  # alerts the user that a .msp file for each component of the lib_metadata list was created, so the user uploads once at a time to NIST MSSearch, validate and answer the folowwing pop-up asking which index of mass spectra is the component
  # tcltk::tkmessageBox(title = "Library development", message = "For each file processed, a '.msp' file was created in the files directory. Upload to NIST MS Search and select the mass spectrum correspondent to the standard.", icon = "info", type = "ok")
  dlg_message("For each file processed, a '.msp' file was created in the files directory. Upload to NIST MS Search and select the mass spectrum correspondent to the standard.")$res
  
  # ppslist is the pslist corrected
  ppslist <- pslist

  # validating spectra
  done <- 2
  while (done == 2) {
    for (i in 1:length(pslist)) {
      id <- dlg_list(c(1:length(pslist[[i]]), "None"), multiple = FALSE, title = names(pslist)[[i]])$res
      if (id == "None") {
        ppslist[[i]] <- "NA"
      } else {
        ppslist[[i]] <- pslist[[i]][[as.numeric(id)]]
      }
    }
    done <- menu(c("Continue", "Do it again"), graphics = TRUE, title = "Validation done!")
  }

  # generate final .msp file
  spectra <- list()
  for (i in 1:length(ppslist)) {
    x <- cbind(ppslist[[i]]@spectrum[, 1], (ppslist[[i]]@spectrum[, 2] / max(ppslist[[i]]@spectrum[, 2]))) # padronizo dividindo todas as intensidades de um mesmo espectro pela maior intensidade no mesmo (fica tipo 1 e 0,X ou seja, porcentagens)
    x <- data.frame(x)
    colnames(x) <- c("mz", "into")
    spectra[[i]] <- x
    rm(x)
  }
  result <- metaMS::construct.msp(spectra, extra.info = NULL)
  for (i in 1:length(result)) {
    result[[i]]$id <- ppslist[[i]]@id
    result[[i]]$rt <- ppslist[[i]]@rt
    result[[i]]$Name <- names(ppslist)[i]
    result[[i]]$Formula <- lib_metadata[grep(names(ppslist)[i], lib_metadata[, 1], value = FALSE), "formula"]
    result[[i]]$MW <- lib_metadata[grep(names(ppslist)[i], lib_metadata[, 1], value = FALSE), "exact.mass"]
    result[[i]]$CAS <- lib_metadata[grep(names(ppslist)[i], lib_metadata[, 1], value = FALSE), "CAS"]
    result[[i]]$ChemSpiderID <- lib_metadata[grep(names(ppslist)[i], lib_metadata[, 1], value = FALSE), "ChemSpiderID"]
    result[[i]]$Class <- lib_metadata[grep(names(ppslist)[i], lib_metadata[, 1], value = FALSE), "class"]
    result[[i]]$Date <- as.character(Sys.Date())
    result[[i]]$Instrument_type <- Instrument_type
    result[[i]]$Comments <- paste0("Column class: ", paste0("Standard ", column), "; ", "ProgramType: ", prog)
    result[[i]]$Ion_mode <- polarity
    if ("RI" %in% colnames(lib_metadata)) {
      result[[i]]$RI <- lib_metadata[grep(names(ppslist)[i], lib_metadata[, 1], value = FALSE), "RI"]
    }
    if ("InChIKey" %in% colnames(lib_metadata)) {
      result[[i]]$InChIKey <- lib_metadata[grep(names(ppslist)[i], lib_metadata[, 1], value = FALSE), "InChIKey"]
    }
  }
  names(result) <- names(ppslist)
  metaMS::write.msp(result, paste0("Library_", Sys.Date(), ".msp"), newFile = TRUE)
  save.image(paste0(Sys.Date(), ".RData"))
  dlg_message("Done! Internal library development finalized. A '.msp' file was created in your directory for conversion to NIST MS Search Library.")
}
