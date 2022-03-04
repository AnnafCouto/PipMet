#' lib_build
#'
#' Main function for GC internal library building. This function process files from internal standards or FAMEs and build a .msp library to be uploaded in NIST MSSearch software. A table of metadata is required containing information about each of the compounds to be processed, such as CAS and SpiderChemID numbers and retention index. If the user does not provide retention index information, the function will search in NIST database for a retention index through CAS number. Therefore, the CAS number must be filled.
#' @keywords library, internal standards, FAMEs
#' @export
#' @return None.
#' @importFrom methods as new
#' @importFrom svDialogs dlgInput
#' @importFrom ddpcr quiet
#' @importFrom utils choose.dir menu read.csv select.list write.csv write.table
#' @importFrom metaMS addRI write.msp
#' @importFrom BiocParallel register SerialParam SnowParam MulticoreParam
#' @importFrom tcltk tkmessageBox
#' @importFrom svDialogs dlgInput dlg_message
#' @importFrom webchem nist_ri as.cas
#' @importFrom MSnbase readMSData filterRt
#' @importFrom xcms MatchedFilterParam findChromPeaks refineChromPeaks FilterIntensityParam
#' @importFrom methods as
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr
#' @importFrom metaMS construct.msp write.msp addRI
#' @importFrom tools file_path_as_absolute
#' @examples
#' \dontrun{
#' lib_build()
#' }
#'
lib_build <- function() {
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
  register(parallel, default = TRUE)

  # get files information
  myDir <- choose.dir(default = getwd(), caption = "Please, select the Samples directory, should be C:/Users/_/Samples")
  setwd(myDir)
  myDir <- getwd()
  extensao <- menu(c(".mzML", ".mzXML"), graphics = TRUE, title = "Files extension:")
  if (extensao == 1) {
    extensao <- ".mzML"
  } else {
    extensao <- ".mzXML"
  }

  # create, fill and read (metadata table)
  files <- list.files(myDir, full.names = TRUE, pattern = extensao, recursive = TRUE)
  metadata <- matrix(nrow = length(files), ncol = 6)
  colnames(metadata) <- c("Name", "Formula", "monoMW", "RTman", "StdFile", "Class", "CAS", "ChemSpiderID", "RI", "Instrument", )
  metadata[, "StdFile"] <- files
  metadata[, "Name"] <- sub(basename(files), pattern = extensao, replacement = "", fixed = TRUE)
  write.csv(metadata, "Lib_Metadata.csv", row.names = FALSE)
  tkmessageBox(title = "Metadata", message = "A file 'Lib_Metadata.csv' was created in your directory. Fill the sheet before continuing. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.", icon = "info", type = "ok")
  while (file.exists("Lib_Metadata.csv") == FALSE) {
    dlg_message("A file 'Lib_Metadata.csv' was created in you directory. Fill the sheet before continuing. Make sure to fill the CAS number of the compounds in order to look for missing RI in literature. After filling the sheet, press 'ok'.", type = "ok")
  }
  metadata <- read.csv("Lib_Metadata.csv", na.string = c("NA", ""), colClasses = "character", sep = ";")
  if (!sum((is.na(metadata))) == 0) {
    metadata <- metadata[, -which(is.na(metadata), arr.ind = TRUE)[, 2]]
  } # remove empty columns

  for (i in nrow(metadata)) {
    metadata[i, 5] <- file_path_as_absolute(metadata[i, 5])
  }

  # more info
  col_pol <- select.list(c("polar", "non-polar"), preselect = NULL, multiple = FALSE, title = "Polarity of columns", graphics = TRUE)
  temp <- select.list(c("isothermal", "ramp", "custom"), preselect = NULL, multiple = FALSE, title = "Temperatura program", graphics = TRUE)
  rt <- select.list(c("kovats", "linear", "alkane", "lee"), preselect = NULL, multiple = FALSE, title = "Rentention Index type", graphics = TRUE)
  colEnergy <- dlgInput("Colision Energy (eV): ", 70)$res

  # check if RI is present
  for (i in 1:length(metadata)) {
    metadata[i, "RI"] <- mean(nist_ri(as.cas(metadata[i, "CAS"]),
      from = "cas",
      type = rt,
      polarity = col_pol,
      temp_prog = temp,
      verbose = getOption("verbose")
    )$RI)
  }

  # read files and process
  dados_brutos <- list()
  for (i in 1:nrow(metadata)) {
    dados_brutos[[i]] <- readMSData(metadata[i, 5], mode = "onDisk")
  }
  xdata <- lapply(dados_brutos, FUN = findChromPeaks, param = MatchedFilterParam(fwhm = 5, binSize = 0.5, steps = 2, mzdiff = 0.5, snthresh = 2, max = 500))
  # filter RT time based on de time provided for each of the standards +/- 10 seconds
  for (i in 1:nrow(metadata)) {
    xdata[[i]] <- filterRt(xdata[[i]], c(metadata[i, 4] - 10, metadata[i, 4] + 10))
  }
  for (i in 1:length(xdata)) {
    xset[[i]] <- as(xdata[[i]], "xcmsSet")
  }
  # group peak in spectra using retention time and peak correlation information
  polarity <- menu(c("Positive", "Negative"), graphics = TRUE, title = "Polarity:")
  if (polarity == 1) {
    an <- lapply(xset, FUN = xsAnnotate, polarity = "positive")
  } else {
    an <- lapply(xset, FUN = xsAnnotate, polarity = "negative")
  }
  an <- lapply(xset, FUN = xsAnnotate, polarity = "positive")
  anF <- lapply(an, FUN = groupFWHM, perfwhm = 1)
  anIC <- lapply(anF, FUN = groupCorr, calcIso = FALSE)
  pslist <- apply(anIC, FUN = exctractSpectra)
  names(pslist) <- names(anIC) <- names(an) <- names(anF) <- names(xset) <- names(xdata) <- names(dados_brutos) <- metadata[, 1]

  # standardize the intensities in the spectrum
  dir.create("Validation")
  setwd("Validation")
  for (ii in 1:length(pslist)) {
    spectra <- list()
    for (i in 1:length(pslist[ii])) {
      x <- cbind(pslist[[i]]@spectrum[, 1], (pslist[[i]]@spectrum[, 2] / max(pslist[[i]]@spectrum[, 2]))) # pslistronizo dividindo todas as intensidades de um mesmo espectro pela maior intensidade no mesmo (fica tipo 1 e 0,X ou seja, porcentagens)
      x <- data.frame(x)
      colnames(x) <- c("mz", "into")
      spectra[[i]] <- x
      rm(x)
    }
    result <- construct.msp(spectra, extra.info = NULL)
    write.msp(result, paste0(" - ", i, "Database_[date].msp"), newFile = TRUE)
  }

  dlg_message("In 'Validation' folder, you will find .msp files to upload to NIST MSSearch and define wich spectra (by index) represents the compound.", type = "ok")

  # validate the spectra
  done <- 2
  while (done == 2) {
    for (i in 1:length(pslist)) {
      pslist[[i]] <- pslist[[i]][[as.integer(menu(c(1:length(pslist[[i]])), graphics = TRUE, title = paste0("Wich is ", names(pslist)[i])))]]
    }
    done <- menu(c("Continue", "Do it again"), graphics = TRUE, title = "Validation done!")
  }

  # write .msp structure for each spectrum
  spectra <- list()
  for (i in 1:length(pslist)) {
    x <- cbind(pslist[pslist]@spectrum[, 1], (pad[[i]]@spectrum[, 2] / max(pad[[i]]@spectrum[, 2]))) # padronizo dividindo todas as intensidades de um mesmo espectro pela maior intensidade no mesmo (fica tipo 1 e 0,X ou seja, porcentagens)
    x <- data.frame(x)
    colnames(x) <- c("mz", "into")
    spectra[[i]] <- x
    rm(x)
  }
  result <- construct.msp(spectra, extra.info = NULL)
  rm(i)
  for (i in 1:length(result)) {
    result[[i]]$id <- pslist[[i]]@id
    result[[i]]$rt <- pslist[[i]]@rt
    result[[i]]$Name <- names(pslist)[i]
    result[[i]]$Formula <- metadata[grep(names(pslist)[i], metadata[, 1], value = FALSE), "Formula"]
    result[[i]]$monoMW <- metadata[grep(names(pslist)[i], metadata[, 1], value = FALSE), "monoMW"]
    result[[i]]$CAS <- metadata[grep(names(pslist)[i], metadata[, 1], value = FALSE), "CAS"]
    result[[i]]$ChemSpiderID <- pslimetadata[grep(names(pad)[i], metadata[, 1], value = FALSE), "ChemSpiderID"]
    if ("Class" %in% colnames(metadata)) {
      result[[i]]$Class <- metadata[i, "Class"]
    } else {
      result[[i]]$Class <- "Standard"
    }
    result[[i]]$Date <- as.character(Sys.Date())
    result[[i]]$RI <- metadata[grep(names(pslist)[i], metadata[, "RI"], value = FALSE), "CAS"]
    result[[i]]$CollisionEnergy <- paste0(colEnergy, 'eV')
  }
  names(result) <- names(pslist)
  write.msp(result, "Database_[date].msp", newFile = TRUE)

  # done!
  dlg_message("Library built!", type = "ok")
}
