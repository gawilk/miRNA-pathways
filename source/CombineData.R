# A bunch of functions for combining TCGA data into RData objects for analysis
 
#====================================================================================
# create ExpressionSet from R object data file
#====================================================================================

CreateEset <- function(path, clinical.file, type = "RNA", cancer = "BRCA") {
  # Creates bioconductor ExpressionSet for gene expression, 
  # miRNA Seq, or protein expression data.
  # 
  # Args: 
  #   path: path where data set is located
  #   clinical.file: path for clinical data file (like "...clinical_patient_brca.txt") 
  #   type: data type, either "gene", "RNA", or "miRNA"
  #   cancer: cancer type
  #
  # Returns: 
  #   ExpressionSet with assayData being genes/miRNA/RNA vs samples
  #   Includes clinical pData: samples vs covariates
  require(Biobase)
  obj <- load(path)
  obj.df <- get(obj)
  data.mat <- as.matrix(obj.df)
  class(data.mat) <- "numeric"
  clinical <- read.delim(clinical.file, header = FALSE, skip = 3, sep = "\t")
  header <- scan(clinical.file, nlines = 1, skip = 1, what = character())
  names(clinical) <- header   
  clinical.header <- c("tcga_barcode", header)
  covariates <- matrix(nrow = ncol(data.mat), ncol = length(clinical.header))
  rownames(covariates) <- colnames(data.mat)
  colnames(covariates) <- clinical.header
  for (i in 1:length(rownames(covariates))) {
    fullbarcode <- rownames(covariates)[i]
    patid <- substr(fullbarcode, 1, 12)
    if (sum(clinical$bcr_patient_barcode == patid) > 0) {
      clinical.row <- clinical[clinical$bcr_patient_barcode == patid, ]
      covariates[i, 1] <- fullbarcode
      covariates[i, 2:length(clinical.header)] <- as.matrix(clinical.row)
    } else {
      covariates[i, 1] <- fullbarcode
      covariates[i, 2:length(clinical.header)] <- rep("[Not Available]", 
                                                      length(clinical.header) - 1)
    }
  }
  sanity.check <- all(colnames(data.mat)==rownames(covariates))
  if (sanity.check) {
    eset <- new("ExpressionSet", exprs = data.mat)
    covariates.adf <- as(as.data.frame(covariates), "AnnotatedDataFrame")
    phenoData(eset) <- covariates.adf
    experimentData(eset) <- new("MIAME",
                                name = paste("TCGA", cancer, type, "level 3 data"),
                                lab = "TCGA",
                                url = "http://cancergenome.nih.gov/")  
    return(eset) 
  }
}

#====================================================================================
#RNA Seq (for gene expression)
#====================================================================================

CombineRNASeqData <- function(directory, output, scaled = TRUE, filemap) {
  # Combines all TCGA level 3 RNA Seq values for a given cancer type
  # Note: only considers files containing "*.genes.normalized_results" 
  #       ignores isoforms and other measures
  #
  # Args: 
  #   directory: path of directory where all RNA Seq .txt files live (1 file per sample)
  #   output: string indicating which txt file in RSEM output to use: 
	#			      "*.genes.normalized_results" files for normalized_count
	#			      "*.genes.results" for scaled_esimatet or raw_count
  #            must include "*" in front of string
  #   scaled: use column scaled_estimate (tau) if TRUE, raw_count if FALSE (default TRUE)
	#		filemap: path of sample map file, maps filenames to TCGA sample barcodes
  #
  # Returns: 
  #   The data frame of combined data: RNASeq genes assayed (rows) by samples (columns)
  curwd <- getwd()
  setwd(directory)
  filenames <- list.files()
  files <- filenames[grep(output, filenames)]
  df.files <- lapply(files, read.delim, sep = "\t", header = TRUE)
  setwd(curwd)
  samplemap <- read.delim(filemap, header = TRUE, sep = "\t")
  names(df.files) <- files
  if (output == "*.genes.results") {
    if (scaled) {
      df.files <- lapply(df.files, subset, select = c(gene_id, scaled_estimate))
	  } else {
      df.files <- lapply(df.files, subset, select = c(gene_id, raw_count))
    }
  }
  for (i in 1:length(df.files)) {
    file <- names(df.files)[i]
    barcode <- samplemap$barcode.s.[samplemap$filename == file]
    colnames(df.files[[i]]) <- c("gene_id", as.character(barcode))
  }
  RNASeq <- Reduce(function(x, y) merge(x, y, by = "gene_id", sort = FALSE), df.files)
  RNASeq$gene_id <- gsub("\\Q|\\E[0-9]*", "", RNASeq$gene_id)
  RNASeq <- RNASeq[!grepl("\\Q?\\E", RNASeq$gene_id), ] 
  RNASeq <- RNASeq[!duplicated(RNASeq$gene_id), ]
  rownames(RNASeq) <- RNASeq$gene_id
  RNASeq <- subset(RNASeq, select = -c(gene_id))
  return(RNASeq)
}

#====================================================================================
# miRNA Seq
#====================================================================================

CombinemiRNASeqData <- function(file, ref = "FALSE") {
  # Combines all TCGA level 3 miRNA Seq values for a given cancer.
  # Note: Expression values are reads_per_million_miRNA_mapped 
  #       which is the normalized read count reported for each miRNA_ID
  #
  # Args: 
  #   file: path where all miRNA Seq *.txt files are located (1 file per sample)
  #   ref: if TRUE, specify reference genome (hg19)
  #
  # Returns: 
  #   The data frame of combined data: miRNA IDs assayed by samples
  setwd(file)
  filenames <- list.files()
  n <- length(filenames)
  if (ref) {
    files.miRNA <- filenames[grep("hg19.mirna.quantification", filenames)]
    barcodes <- gsub(".hg19.mirna.quantification.txt", "", files.miRNA)
  } else {
    files.miRNA <- filenames[grep("mirna.quantification", filenames)]
    barcodes <- gsub(".mirna.quantification.txt", "", files.miRNA)
  }
  df.list <- lapply(files.miRNA, read.delim, sep = "\t", header = TRUE)
  keeps <- c("miRNA_ID", "reads_per_million_miRNA_mapped")
  trimmed <- lapply(df.list, function(x) x[keeps])
  n <- length(trimmed)
  for (i in 1:n) {
    colnames(trimmed[[i]]) <- c("miRNA_ID", barcodes[i])
  }
  miRNASeq <- Reduce(function(x, y) merge(x, y, by = "miRNA_ID", sort = FALSE), trimmed)
  miRNASeq <- miRNASeq[!duplicated(miRNASeq$miRNA_ID), ]
  rownames(miRNASeq) <- miRNASeq$miRNA_ID
  miRNASeq <- subset(miRNASeq, select = -c(miRNA_ID))
  return(miRNASeq)
}
