#' Normalize
#'
#' This function normalize data.
#' @export
#' @param myDir Path to the directory of work.
#' @param pslist List of spectra.
#' @param metadata A matrix or data.frame with metadata information about samples. Include, at least 'sample' and 'file' columns with name of sample and its path, respectively. More information can be added in new columns, such as 'group', 'class', 'biorep' and 'tecrep'.
#' @param anIC A 'xsAnnotate' CAMERA object with grouped spectra.
#' @param pre_anno A table with annotations for spectra. or path to .csv file.
#' @param pic_extension Character. Pictures format to generate. Supported = '.tiff', '.png'. Default to c('.tiff', '.png').
#' @param example Logical. If is example, pop-ups won't appear. Default to FALSE.
#' @param derivatization Character. Kind of derivatization the samples were prepared with. Supported are 'Trimethylsilyl' and 'None'. If NULL, the user will be asked. Default to 'NULL'.
#' @return A matrix of all spectra, with their annotation (if available), most intense peak m/z and its intensities in every sample.
#' @importFrom grDevices dev.off pdf png tiff boxplot.stats
#' @importFrom graphics boxplot grid legend par text
#' @importFrom svDialogs dlg_message
#' @importFrom methods as new
#' @importFrom stats cor sd t.test var
#' @importFrom utils choose.dir menu read.csv select.list write.csv write.table
#' @importFrom tools file_path_as_absolute
#' @import NormalyzerDE
#' @import ggplot2
#' @examples
#' \donttest{
#' \dontrun{
#' pre_anno <- read.csv(system.file("extdata", "pre_anno.csv", package = "PipMet"))
#' load(system.file("extdata", "pslist.RData", package = "PipMet"))
#' load(system.file("extdata", "anIC.RData", package = "PipMet"))
#' load(system.file("extdata", "metadata.RData", package = "PipMet"))
#' normalized <- normalize_data(
#'   anIC,
#'   pslist,
#'   metadata,
#'   myDir = "~/",
#'   pre_anno = pre_anno,
#'   example = TRUE
#' )
#' }
#' }
normalize_data <- function(anIC, pslist, metadata, myDir, pre_anno, example, pic_extension = c('.tiff', '.png'), derivatization = NULL) {

  # check if representative ions are ok
  okay <- 2
  while (okay == 2) {
    setwd(myDir)
    quant <- matrix(nrow = nrow(pre_anno), ncol = 5 + (2 * nrow(metadata)))
    colnames(quant) <- c("id", "Fragment Ion (m/z Quant)", "Compound Name", "Chemical Formula", "Metabolic Class", metadata$sample, metadata$sample)
    quant[, "id"] <- as.numeric(pre_anno[, "id"])
    quant[, "Compound Name"] <- pre_anno[, "annotation"]

    # set to folder
    if (!dir.exists("Statistics") == TRUE) {
      dir.create("Statistics")
    }
    setwd("Statistics")

    # ask about derivatizations
    if (example == TRUE) {
      derivatization <- 'Trimethylsilyl'
    } else {
      if (is.null(derivatization)) {
        #derivatization <- menu(c('Trimethylsilyl','None'), graphics = TRUE, title = "Sort of derivatization: ")
        derivatization <- dlg_list(c('Trimethylsilyl','None'), multiple = FALSE, title = "Sort of derivatization: ")$res
      }
    }

    for (i in 1:nrow(quant)) {
      temp <- anIC@pspectra[[as.integer(quant[i, 1])]]
      if (derivatization == 'Trimethylsilyl' && 73 %in% round(pslist[[i]]@spectrum[, 1])) { # if derivatization='Trimethylsilyl', ion m/z 73 is present
        x <- sort(pslist[[i]]@spectrum[-which(round(pslist[[i]]@spectrum[, 1]) == 73), 2], decreasing = TRUE) # remove m/z 73 from possibilities of representative ion
      } else {
        x <- sort(pslist[[i]]@spectrum[, 2], decreasing = TRUE)
      }
      y <- rbind(anIC@groupInfo[temp, ])
      z <- rbind(y[which(y == x[1], arr.ind = TRUE)[, 1], ], y[which(y == x[2], arr.ind = TRUE)[, 1], ])
      quant[i, 2] <- paste0(z[1, 1], ", ", z[2, 1])
      quant[i, 6:(5 + nrow(metadata))] <- z[1, which(colnames(z) == "X1"):which(colnames(z) == paste0("X", nrow(metadata)))]
      quant[i, (6 + nrow(metadata)):(5 + 2 * nrow(metadata))] <- z[2, which(colnames(z) == "X1"):which(colnames(z) == paste0("X", nrow(metadata)))]
    }
    mat <- quant[, 6:ncol(quant)]
    mat <- apply(mat, c(1, 2), FUN = as.numeric)

    # write non normalized data into .csv file
    write.csv(quant, "NotNormalized_quantification.csv", row.names = FALSE, na = "")

    # Normalyzer: evaluation and picking of normalization method for the spectra intensities (only most intense peak)
    x <- as.data.frame(mat)
    x[x == 0] <- NA
    x[x > 0 & x < 1] <- 0
    write.table(x[, 1:nrow(metadata)], "data.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(metadata, "design.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

    # set up for NormalyzerDE
    designFp <- file_path_as_absolute("design.tsv")
    dataFp <- file_path_as_absolute("data.tsv")

    # ask condition to compare from the metadata table
    if (example == TRUE) {
      x <- 2
    } else {
      x <- menu(colnames(metadata), graphics = TRUE, title = "Choose conditions (from metadata table) to group for normalization: ")
    }
    normalyzer(jobName = "Normalyzer_results", designPath = designFp, dataPath = dataFp, outputDir = myDir, sampleColName = "sample", groupColName = colnames(metadata)[x], requireReplicates = FALSE)

    # Pick method and apply
    fill <- list.files(paste0(myDir, "/Normalyzer_results"), full.names = TRUE, pattern = "-normalized.txt", recursive = TRUE)
    norms <- sub("-normalized.txt", replacement = "", fixed = TRUE, x = basename(fill))
    # ask method for normalization. Standard method for the example is CycLoess.
    if (example == TRUE) {
      dlg_message("Normalization done! Check the report from NormalyzerDE. As example, the chosen one is CycLoess method.")$res
      bestNormMat <- 1
      okay <- 1
    } else {
      bestNormMat <- menu(norms, graphics = TRUE, title = "Choose the best normalization")
    }
    mat[mat == 0] <- NA
    mat[mat > 0 & mat < 1] <- 0
    if (norms[bestNormMat] == "CycLoess") {
      mat <- performCyclicLoessNormalization(mat)
    }
    if (norms[bestNormMat] == "GI") {
      mat <- globalIntensityNormalization(mat)
    }
    if (norms[bestNormMat] == "log2") {
      mat <- log2(mat)
    }
    if (norms[bestNormMat] == "mean") {
      mat <- meanNormalization(mat)
    }
    if (norms[bestNormMat] == "median") {
      mat <- medianNormalization(mat)
    }
    if (norms[bestNormMat] == "Quantile") {
      mat <- performQuantileNormalization(mat)
    }
    if (norms[bestNormMat] == "RLR") {
      mat <- performGlobalRLRNormalization(mat)
    }
    if (norms[bestNormMat] == "VSN") {
      mat <- performVSNNormalization(mat)
    }
    n <- as.data.frame(cbind(quant[, 1:5], mat))

    # Confirm if the chosen ions are indeed from the same spectrum
    # if everything is ok, use only first ion for the rest of statistics
    y <- data.frame()
    for (ii in 6:(5 + nrow(metadata))) {
      for (i in 1:nrow(n)) {
        y[i, (ii - 5)] <- as.double(n[i, ii]) / as.double(n[i, (ii + nrow(metadata))])
      }
    }
    y[, nrow(metadata) + 1] <- apply(y[, 1:nrow(metadata)], MARGIN = 1, FUN = var, na.rm = TRUE)
    y[, nrow(metadata) + 2] <- apply(y[, 1:nrow(metadata)], MARGIN = 1, FUN = sd, na.rm = TRUE)
    rownames(y) <- n$id
    colnames(y) <- c(metadata$sample, "Variance", "sd")
    a <- ggplot(y, aes(x = y[, nrow(metadata) + 1])) +
      geom_histogram() +
      ggtitle("Variance Distribution") +
      xlab("Variance") +
      ylab("Frequency") +
      theme(rect = element_rect(fill = "transparent")) # histograma da variância. O ideal é que tenha pouca variância alta
    b <- ggplot(y, aes(x = y[, nrow(metadata) + 2])) +
      geom_histogram() +
      ggtitle("Standard Deviation") +
      xlab("Standard Deviation") +
      ylab("Frequency") +
      theme(rect = element_rect(fill = "transparent")) # histograma do desvio padrão.

    # plot histogram of variance and standard deviation and boxplots
    # png
    if ('.png' %in% pic_extension) {
    png("variance_dp.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
    gridExtra::grid.arrange(a, b, ncol = 1)
    dev.off()
    png("boxplot.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
    par(mfrow = c(2,1))
    boxplot(y$Variance, main = "Variance", horizontal = TRUE)
    boxplot(y$sd, main = "Standard Deviation", horizontal = TRUE)
    dev.off()}

    # tiff
    if ('.tiff' %in% pic_extension) {
    tiff("variance_dp.tiff", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
    gridExtra::grid.arrange(a, b, ncol = 1)
    dev.off()
    tiff("boxplot.tiff", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
    par(mfrow = c(2,1))
    boxplot(y$Variance, main = "Variance", horizontal = TRUE)
    boxplot(y$sd, main = "Standard Deviation", horizontal = TRUE)
    dev.off()}

    n <- cbind(n[, 1:(5 + nrow(metadata))], y[, c("Variance", "sd")]) # as everything is ok, use only first most intense ions for representative
    write.csv(n, "Normalized_quantification.csv", row.names = FALSE, na = "")

    if (example == FALSE) {
      res <- menu(c("Conclude normalization", "Re-do normalization", "Remove all outliers", "Remove spectra with variance bigger than defined"), graphics = TRUE, title = "Choose the best normalization")
      if (res == 1) {
        okay <- 1
      }
      if (res == 2) {
        okay <- 2
      }
      if (res == 3) {
        y_1 <- y[-which(y$Variance >= boxplot.stats(y$Variance)$stats[5]), ]
        n_1 <- n[-which(y$Variance >= boxplot.stats(y$Variance)$stats[5]), ]

        # png
        if ('.png' %in% pic_extension) {
        png("boxplot_withoutOutliers.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
        par(mfrow = c(2,1))
        boxplot(y_1$Variance, main = "Variance", horizontal = TRUE)
        boxplot(y_1$sd, main = "Standard Deviation", horizontal = TRUE)
        dev.off()}
        # tiff
        if ('.tiff' %in% pic_extension) {
        tiff("boxplot_withoutOutliers.tiff", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
        par(mfrow = c(2,1))
        boxplot(y_1$Variance, main = "Variance", horizontal = TRUE)
        boxplot(y_1$sd, main = "Standard Deviation", horizontal = TRUE)
        dev.off()}

        n_1 <- cbind(n_1[, 1:(5 + nrow(metadata))], y_1[, c("Variance", "sd")]) # as everything is ok, use only first most intense ions for representative
        write.csv(n_1, "Normalized_quantification_withoutOutliers.csv", row.names = FALSE, na = "")
        a <- ggplot(y_1, aes(x = y[, nrow(metadata) + 1])) +
          geom_histogram() +
          ggtitle("Variance Distribution") +
          xlab("Variance") +
          ylab("Frequency") +
          theme(rect = element_rect(fill = "transparent")) # histograma da variância. O ideal é que tenha pouca variância alta
        b <- ggplot(y_1, aes(x = y[, nrow(metadata) + 2])) +
          geom_histogram() +
          ggtitle("Standard Deviation") +
          xlab("Standard Deviation") +
          ylab("Frequency") +
          theme(rect = element_rect(fill = "transparent")) # histograma do desvio padrão.
        # png
        if ('.png' %in% pic_extension) {
        png("variance_dp_withoutOutliers.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
        gridExtra::grid.arrange(a, b, ncol = 1)
        dev.off()}
        # tiff
        if ('.tiff' %in% pic_extension) {
        tiff("variance_dp_withoutOutliers.tiff", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
        gridExtra::grid.arrange(a, b, ncol = 1)
        dev.off()}        
      }

      if (res == 4) {
        th <- as.double(dlgInput("Variance threshold ", "0")$res)
        y_2 <- y[-which(y$Variance >= th), ]
        n_2 <- n[-which(y$Variance >= th), ]
        # png
        if ('.png' %in% pic_extension) {
        png(paste0("boxplot_varianceUpTo", th, ".png"), units = "cm", width = 16, height = 16, res = 900, bg = "NA")
        par(mfrow = c(2,1))
        boxplot(y_2$Variance, main = "Variance", horizontal = TRUE)
        boxplot(y_2$sd, main = "Standard Deviation", horizontal = TRUE)
        dev.off()}
        # tiff
        if ('.tiff' %in% pic_extension) {
        tiff(paste0("boxplot_varianceUpTo", th, ".tiff"), units = "cm", width = 16, height = 16, res = 900, bg = "NA")
        par(mfrow = c(2,1))
        boxplot(y_2$Variance, main = "Variance", horizontal = TRUE)
        boxplot(y_2$sd, main = "Standard Deviation", horizontal = TRUE)
        dev.off()}
        n_2 <- cbind(n_2[, 1:(5 + nrow(metadata))], y_2[, c("Variance", "sd")]) # as everything is ok, use only first most intense ions for representative
        write.csv(n_2, paste0("Normalized_quantification_varianceUpTo", th, ".csv"), row.names = FALSE, na = "")
      }
      okay <- menu(c("OK, keep going", "No, re-do normalization"), graphics = TRUE, title = "Are the results ok?")
    }
  }

  # set to main folder
  setwd(myDir)

  # return results
  return(n)
}
