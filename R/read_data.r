#' Read data
#'
#' This function creates a sheet for the user to fill with the experiment design in order to create categories for the files to be processed. Also creates images of preprocessed data.
#' @keywords metadata
#' @export
#' @param EIC Numeric. 1 = there are ions to monitor through the processing. 2 = there are none. Default to 2.
#' @param ions List with sublist mz = mz (numeric) of the monitored ion and rt = retention time of monitored ion (numeric). To the 'rt' will be added and subtracted 5 seconds. Default to null.
#' @importFrom grDevices dev.off png tiff
#' @importFrom graphics boxplot legend par
#' @importFrom methods as new
#' @importFrom utils choose.dir menu read.csv select.list write.csv write.table
#' @importFrom svDialogs dlgInput dlg_message
#' @importFrom xcms chromatogram
#' @importFrom MSnbase readMSData fromFile intensity bin
#' @importFrom pheatmap pheatmap
#' @importFrom tcltk tkmessageBox
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' read <- read_data()
read_data <- function(EIC = 2, ions = NULL) {

  # ask user about samples and folders path and create a new folder named after a "Project"
  sample_dir <- choose.dir(default = getwd(), caption = "Please, select the Samples directory, should be C:/Users/_/Samples")
  setwd(sample_dir)
  sample_dir <- getwd()
  myDir <- dlgInput("Name your project", Sys.info()["user"])$res
  dir.create(myDir)
  setwd(myDir)
  myDir <- getwd()
  extensao <- menu(c(".mzML", ".mzXML"), graphics = TRUE, title = "Files extension:")
  if (extensao == 1) {
    extensao <- ".mzML"
  } else {
    extensao <- ".mzXML"
  }

  # set metadata table up
  files <- list.files(sample_dir, full.names = TRUE, pattern = extensao, recursive = TRUE)
  metadata <- matrix(nrow = length(files), ncol = 6)
  colnames(metadata) <- c("sample", "group", "class", "tec_rep", "bio_rep", "file")
  metadata[, "file"] <- files
  metadata[, "sample"] <- sub(basename(files), pattern = extensao, replacement = "", fixed = TRUE)
  write.csv(metadata, "metadata.csv", row.names = FALSE)
  tkmessageBox(title = "Phenodata", message = "A file 'metadata.csv' was created in your directory. Fill the sheet before continuing. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.", icon = "info", type = "ok")
  while (file.exists("metadata.csv") == FALSE) {
    dlg_message("A file 'metadata.csv' was created in you directory. Fill the sheet before continuing. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.", type = "ok")
  }
  metadata <- read.csv("metadata.csv", na.string = c("NA", ""), colClasses = "character", sep = ";")
  if (!sum((is.na(metadata))) == 0) {
    metadata <- metadata[, -which(is.na(metadata), arr.ind = TRUE)[, 2]]
  } # remove empty columns

  # create 'metadata$all' 1 and 2 for identification
  x <- colnames(metadata)
  for (i in c("sample", "tec_rep", "bio_rep", "file", "all", "all2")) {
    if (i %in% x) {
      x <- x[-which(x == i)]
    }
  }
  for (i in 1:nrow(metadata)) {
    f <- metadata[i, x[[1]]]
    for (ii in 2:length(x)) {
      f <- paste0(f, "_", metadata[i, x[[ii]]])
    }
    metadata[i, "all"] <- f
    metadata$all2[i] <- paste0(i, " - ", metadata$all[i])
  }

  # colors for each column in metadata except 'sample' and 'tec_rep'
  colors <- vector(mode = "list", length = length(x))
  names(colors) <- x
  for (i in 1:length(x)) {
    colors[[i]] <- list(metadata[, x[i]], paste0(RColorBrewer::brewer.pal(length(unique(metadata[, x[[i]]])), "Set1")[1:length(unique(metadata[, x[[i]]]))], "60"))
    names(colors[[i]][[2]]) <- c(unique(metadata[, x[i]]))
    names(colors[[i]]) <- c(x[[i]], paste0(x[[i]], "_colors"))
  }

  # read data into R
  dados_brutos <- readMSData(metadata$file, pdata = new("NAnnotatedDataFrame", metadata), mode = "onDisk")

  # images of pre-processing
  dir.create("Visualization_results")
  setwd("Visualization_results")
  bpc <- chromatogram(dados_brutos, aggregationFun = "max")
  tic <- chromatogram(dados_brutos, aggregationFun = "sum")

  # chromatograms
  for (i in 1:length(colors)) {
    # tiff
    tiff(paste0(names(colors)[i], "_chromatograms.tiff"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
    par(mfrow = c(2, 1))
    plot(bpc, col = colors[[i]][[2]][colors[[i]][[1]]], main = "Base Peak Chromatogram")
    legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
    plot(tic, col = colors[[i]][[2]][colors[[i]][[1]]], main = "Total Ion Current Chromatogram")
    legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
    dev.off()
    # png
    png(paste0(names(colors)[i], "_chromatograms.png"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
    par(mfrow = c(2, 1))
    plot(bpc, col = colors[[i]][[2]][colors[[i]][[1]]], main = "Base Peak Chromatogram")
    legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
    plot(tic, col = colors[[i]][[2]][colors[[i]][[1]]], main = "Total Ion Current Chromatogram")
    legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
    dev.off()
  }

  # boxplot of total ion current
  tic_por_arquivo <- split(tic(dados_brutos), f = fromFile(dados_brutos))
  for (i in 1:length(colors)) {
    tic_por_arquivo <- split(tic(dados_brutos), f = fromFile(dados_brutos))
    # tiff
    tiff(paste0(names(colors)[i], "_ticBoxplot.tiff"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
    boxplot(tic_por_arquivo, col = colors[[i]][[2]][colors[[i]][[1]]], ylab = "intensity", xlab = "sample", main = "Total ion current")
    dev.off()
    # png
    png(paste0(names(colors)[i], "_ticBoxplot.png"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
    boxplot(tic_por_arquivo, col = colors[[i]][[2]][colors[[i]][[1]]], ylab = "intensity", xlab = "sample", main = "Total ion current")
    dev.off()
  }
  rm(tic_por_arquivo)

  # cluster
  tic_bin <- bin(tic, binSize = 1)
  cl <- do.call(cbind, lapply(tic_bin, intensity))
  cl[cl == 0] <- NA
  cormat <- cor(log2(cl), use = "pairwise.complete.obs")
  colnames(cormat) <- rownames(cormat) <- metadata$all2
  # for each set of colors (conditions of experiment)
  for (i in 1:length(colors)) {
    ann <- data.frame(colors[[i]][[1]])
    colnames(ann) <- names(colors)[i]
    rownames(ann) <- metadata$all2
    ant <- list(colors[[i]][[2]])
    names(ant) <- names(colors)[i]
    # tiff
    tiff(paste0(names(colors)[i], "_cluster.tiff"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
    pheatmap(cormat, annotation = ann, annotation_color = ant, border_color = "NA", cluster_rows = FALSE, )
    dev.off()
    # png
    png(paste0(names(colors)[i], "_cluster.png"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
    pheatmap(cormat, annotation = ann, annotation_color = ant, border_color = "NA", cluster_rows = FALSE, )
    dev.off()
  }
  rm(tic_bin)

  # extracted ion chromatogram based on mz and rt asked previously by user
  dir.create("Monitoring ions")
  setwd("Monitoring ions")

  if (EIC == 1) {
    for (ii in 1:length(ions)) {
      crom <- chromatogram(dados_brutos, rt = c(as.numeric(ions[[ii]][["rt"]] - 5), as.numeric(ions[[ii]][["rt"]] + 5)), mz = as.numeric(ions[[ii]][["mz"]]))
      for (i in 1:length(colors)) {
        # tiff
        tiff(paste0(names(colors)[i], "_", ii, "_prePross_EIC.tiff"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
        plot(crom, col = colors[[i]][[2]][colors[[i]][[1]]])
        legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
        dev.off()

        # png
        png(paste0(names(colors)[i], "_", ii, "_prePross_EIC.png"), units = "cm", width = 16, height = 16, res = 1500, bg = "NA")
        plot(crom, col = colors[[i]][[2]][colors[[i]][[1]]])
        legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]], fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
        dev.off()
      }
    }
  }

  # return to main folder
  setwd(myDir)

  # return results
  return(list(colors, metadata, dados_brutos, myDir))
}
