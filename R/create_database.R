#' In-house database from processed files
#'
#' Function for creating a in-house database from processed GC-MS experiment files. The user may choose to add more information, such as Formula, Monoisotopic Mass, CAS, ChemSpider, InChIKey, PubChem ID, Class and retention index. The function is capable of search and retrieving most of the information from the online databases such as PubChem and NIST uatomatically through 'webchem' package or add manually the information.
#' @keywords database spectra
#' @export
#' @param pslist List of spectra with annotation. Only annotated spectra will be use for the database construction.
#' @param polarity Character. Ion mode of data acquisition ('negative' or 'positive').
#' @return None.
#' @importFrom svDialogs dlg_message dlg_list
#' @importFrom utils read.csv write.csv
#' @importFrom metaMS write.msp
#' @importFrom webchem get_cid pc_prop cts_convert nist_ri
#' @examples
#' \dontrun{
#' load(system.file("extdata", "pslist.RData", package = "PipMet"))
#' create_database(pslist)
#' }
#'
create_database <- function(pslist, polarity) {

  # get only annotated spectra
  pslist <- list()
  count <- 1
  for (i in 1:length(x$apslist)) {
    if (!isEmpty(x$apslist[[i]]@annotation)) {
      pslist[count] <- x$apslist[[i]]
      count <- count + 1
    }
  }

  # add information from internet?
  if (dlg_message("Would you like to add more informations about the compounds?", type = "yesno")$res == "yes") {
    dataInfo <- matrix(nrow = length(pslist), ncol = 10)
    colnames(dataInfo) <- c("Name", "formula", "exact.mass", "rt", "CAS", "ChemSpider", "InChIKey", "PubChem ID", "Class", "RI")
    dataInfo <- as.data.frame(dataInfo)
    for (i in 1:length(pslist)) {
      dataInfo[i, 1] <- pslist[[i]]@annotation
      dataInfo[i, "rt"] <- pslist[[i]]@rt
      dataInfo[i, "RI"] <- pslist[[i]]@RI
      dataInfo[i, "PubChem ID"] <- get_cid(dataInfo[i, 1])$cid
      f <- pc_prop(dataInfo[i, "PubChem ID"])
      if (!is.na(f["CID"])) {
        dataInfo[i, "formula"] <- f[1, "MolecularFormula"]
        dataInfo[i, "InChIKey"] <- f[1, "InChIKey"]
        dataInfo[i, "exact.mass"] <- f[1, "MonoisotopicMass"]
      }
    }
    dataInfo[, "CAS"] <- unlist(cts_convert(dataInfo[, "PubChem ID"], "PubChem CID", "cas", match = "first"))
    # ask information for retention index search
    if (dlg_message("Look for retention index in NIST?", type = "yesno")$res == "yes") {
      Ri <- dlg_list(c("kovats", "linear", "alkane", "lee"), multiple = FALSE, title = "Retention time index")$res
      column <- dlg_list(c("polar", "non-polar"), multiple = FALSE)$res
      prog <- dlg_list(c("isothermal", "ramp", "custom"), multiple = FALSE)$res
      Instrument_type <- dlg_input("Type of instrument of acquisition", "GC-EI-Q")$res
      for (i in 1:length(pslist)) {
        dataInfo[i, "RI"] <- nist_ri(dataInfo[i, "Name"], from = "name", type = Ri, polarity = column, temp_prog = prog)
      }
    }
    # write informations for checking
    write.csv(dataInfo, file = "DatabaseInfo.csv", row.names = FALSE)
    dlg_message("Check and fill the 'DatabaseInfo.csv' file in myDir and press 'OK'.")$res
    dataInfo <- read.csv("DatabaseInfo.csv", na.string = c("NA", ""), colClasses = "character", sep = ",")
  }

  # create spectra in .msp file format
  spectra <- list()
  for (i in 1:length(pslist)) {
    x <- cbind(pslist[[i]]@spectrum[, 1], (pslist[[i]]@spectrum[, 2] / max(pslist[[i]]@spectrum[, 2]))) # pslistronizo dividindo todas as intensidades de um mesmo espectro pela maior intensidade no mesmo (fica tipo 1 e 0,X ou seja, porcentagens)
    x <- data.frame(x)
    colnames(x) <- c("mz", "into")
    spectra[[i]] <- x
    rm(x)
  }
  result <- metaMS::construct.msp(spectra, extra.info = NULL)
  for (i in 1:length(result)) {
    result[[i]]$rt <- dataInfo[i, "rt"]
    result[[i]]$Name <- dataInfo[i, "Name"]
    result[[i]]$Formula <- dataInfo[i, "formula"]
    result[[i]]$MW <- dataInfo[i, "exact.mass"]
    result[[i]]$CAS <- dataInfo[i, "CAS"]
    result[[i]]$ChemSpiderID <- dataInfo[i, "ChemSpiderID"]
    result[[i]]$InChIKey <- dataInfo[i, "InChIKey"]
    result[[i]]$PubChemID <- dataInfo[i, "PubChem ID"]
    result[[i]]$Class <- dataInfo[i, "Class"]
    result[[i]]$Date <- as.character(Sys.Date())
    result[[i]]$RI <- dataInfo[i, "RI"]
    result[[i]]$Instrument_type <- Instrument_type
    result[[i]]$Comments <- paste0("Column class: ", paste0("Standard ", column), "; ", "ProgramType: ", prog)
    result[[i]]$Ion_mode <- polarity
  }
  names(result) <- dataInfo[, "Name"]
  metaMS::write.msp(result, file = paste0("Database_", Sys.Date(), ".msp"), newFile = TRUE)
  dlg_message("Database creation done!")$res
}
