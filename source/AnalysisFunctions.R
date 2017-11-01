# Analysis functions to source

#==============================================================================
# get genes in pathway
#==============================================================================

getPathwayGenes <- function(pathID = NULL) {
  # Args:
  #  pathID: pathway ID(s), if NULL then get all genes in all pathways
  # Returns:
  #   a list of pathway(s), each entry is vector of genes in pathway
  require(Biobase)
  require(org.Hs.eg.db)
  require(annotate)
  require(KEGG.db)
  keggpath <- org.Hs.egPATH2EG
  keggpath.mapped <- mappedkeys(keggpath)
  paths <- as.list(keggpath[keggpath.mapped])
  if (is.null(pathID)) {
    pathways <- sapply(paths, function(X) {
      eg2symb <- getSYMBOL(X, "org.Hs.eg")
      symbols <- as.character(eg2symb)                       
    })
    return(pathways)
  } else {
    pathways <- sapply(as.character(pathID), function(X) {
      eg2symb <- getSYMBOL(paths[[X]], "org.Hs.eg")
      symbols <- as.character(eg2symb)
    }, simplify = FALSE)
    return(pathways)
  }
}

#==============================================================================
# subset gene expression data by KEGG paths
#==============================================================================

subsetGEbyKEGG <- function(eset, pathIDs = NULL) {
  # Args:
  #   eset: gene expression eset
  #   pathIDs: pathways to subset, if NULL then all pathways
  # Returns:
  #   A list of KEGG pathways and their gene expression matrices
  #   (note: matrices are samples by genes, NAs are omitted)
  require(Biobase)
  require(org.Hs.eg.db)
  require(annotate)
  require(KEGG.db)
  pathways <- getPathwayGenes(pathID = pathIDs)
  features <- featureNames(eset) 
  indices <- sapply(pathways, match, features)
  pathway.matrices <- lapply(indices, function(x) {
    eset.noNA <- eset[na.omit(x), ]
    mat.noNA <- t(na.omit(exprs(eset.noNA)))
  })
  return(pathway.matrices)
}

#====================================================================================
# find "optimal" k through Isomap
#====================================================================================

findKisomap <- function(X, K, scale = TRUE, plot = TRUE) {
  # Args:
  #   X: numeric matrix
  #   K: sequence of k values to consider
  #   scale: scale columns of X? TRUE or FALSE
  #   plot: plot Isomap spectra vs. PCA spectra? TRUE or FALSE
  # Returns:
  #   A list of the optimal "k" and the spectral gap ratio
  require(vegan)
  if (scale) {
    X <- scale(X, center = TRUE, scale = TRUE)
  }
  mat <- as.matrix(dist(X))
  PCA <- prcomp(X, scale = FALSE, center = TRUE)
  PCAspect <- PCA$sdev^2 / sum(PCA$sdev^2)
  PCAratio <- diff(PCAspect)[1] / diff(PCAspect)[2]
  ratios <- c()
  spectrum <- list()
  for (i in 1:length(K)) {
    runISOMAP <- tryCatch({
      iso <- isomap(mat, ndim = 3, k = K[i])
      spect <- iso$eig / sum(iso$eig)
      spectrum[[i]] <- spect
      gaps <- diff(spect)
      x <- gaps[1] / gaps[2]
    }, warning = function(war) {
      print(paste("MY WARNING: ", war))
      iso <- isomap(mat, ndim = 3, k = K[i])
      spect <- iso$eig / sum(iso$eig)
      spectrum[[i]] <- spect
      gaps <- diff(spect)
      x <- gaps[1] / gaps[2]
    }, error = function(err) {
      print(paste("MY ERROR: ", err))
      spectrum[[i]] <- NA
      x <- NA
    }, finally = {}) 
    ratios <- c(ratios, runISOMAP)    
  }
  if (plot) {
    colvect <- rainbow(length(K))
    plot(PCAspect[1:6], ylab = "value", xlab = "spectra", ylim = c(0, 1), pch = 19)
    for (j in 1:length(spectrum)) {
      lines(spectrum[[j]], xlab = "", ylab = "", col = colvect[j], type = "b")
    }
    legend("topright", paste("k =", K), col = colvect, pch = 1, cex = 0.6)
  }
  spectdiffs <- ratios / PCAratio
  spectgaps <- cbind(K, spectdiffs)
  index <- which(spectdiffs == max(spectdiffs, na.rm = TRUE), arr.ind = TRUE)
  return(list(table = spectgaps, k = K[index]))
}

#==============================================================================
# separation of phenotypes along PAS (p-vals) 
#==============================================================================

phenoSep <- function(embeddings, sampClasses, classA = "01", classB = "11") {
  # Args:
  #   embeddings: pathway embeddings (samples by embedded coordinates)
  #   sampClasses: vector of sample classes
  #   classA: category of samples to compare, default tumor "01"
  #   classB: another sample category, default tissue-normal "11"
  # Returns:
  #   A vector, each element is p-value of coordinate embedding for phenotype separation
  if (nrow(embeddings) != length(sampClasses)) {
    stop("dimensions don't match")
  }
  out <- apply(embeddings, 2, function(coordinate) {
    pheno <- t.test(coordinate[sampClasses == classA], coordinate[sampClasses == classB])
    pheno
  })
  return(out)
}

#====================================================================================
# find common samples in miRNASeq and pathway data
#====================================================================================

findCommonSamples <- function(miRNASeq, pathMatrix, barcode.len = 15, na.rm = FALSE) {
  # Args:
  #   miRNASeq: miRNASeq matrix, samples (rows) by miRNAs (columns) already labeled
  #   pathMatrix: matrix of pathway values, samples (rows) by pathways (columns) already labeled
  #   barcode.len: number of characters in TCGA barcode to match, default 15
  #   na.rm: if TRUE, remove columns in pathway matrix that have NAs, default FALSE
  # Returns:
  #   A list containing filtered miRNASeq + pathMatrix matrices and samples in common
  mirna <- miRNASeq[!duplicated(substr(rownames(miRNASeq), 1, barcode.len)), ]
  rownames(mirna) <- substr(rownames(mirna), 1, barcode.len)
  paths <- pathMatrix[!duplicated(substr(rownames(pathMatrix), 1, barcode.len)), ]
  rownames(paths) <- substr(rownames(paths), 1, barcode.len)
  commonSamples <- intersect(rownames(mirna), rownames(paths))
  mirna <- mirna[commonSamples, ]
  paths <- paths[commonSamples, ]
  if (na.rm) {
    na.cols <- as.logical(apply(paths, 2, function(x) sum(is.na(x))))
    paths <- paths[, !na.cols]
  }
  return(list(mirna = mirna, paths = paths, commonSamples = commonSamples))
}

#====================================================================================
# compare correlations across samples classes between two matrices
#====================================================================================

compareCorls <- function(X, Y, sampClasses, classA = "01", classB = "11") {
  # Args:
  #   X & Y: matrices to compare correlations, samples must be rows
  #   sampClasses: vector signifying which samples belong to which class 
  #   classA: class of group A, default is "01" for tumor samples
  #   classB: class of group B, default is "11" for adjacent normals
  # Returns:
  #   Matrix of miRNAs and pathways, each element is correlation differences between two classes
  #   (note: check that X & Y are aligned & have same length as sampClasses)
  stopifnot(all(rownames(X) == rownames(Y)), length(sampClasses) == nrow(X))
  # "pairwise" so we don't need to remove NA columns;
  # "spearman" to use rank corls
  if (sum(is.na(X)) == 0 && sum(is.na(Y)) == 0) {
    corl.classA <- cor(X[sampClasses == classA, ], Y[sampClasses == classA, ], use = "everything", method = "spearman")
    corl.classB <- cor(X[sampClasses == classB, ], Y[sampClasses == classB, ], use = "everything", method = "spearman")
  } else {
    corl.classA <- cor(X[sampClasses == classA, ], Y[sampClasses == classA, ], use = "pairwise", method = "spearman")
    corl.classB <- cor(X[sampClasses == classB, ], Y[sampClasses == classB, ], use = "pairwise", method = "spearman")
  }
  return(corl.classA - corl.classB)
}

#====================================================================================
# compute p-values for each miRNA-pathway pair
#====================================================================================

computePvals <- function(X, Y, trueCorlDiffs, sampClasses, classA = "01", classB = "11",
                         resamplings = 1000) {
  # Args:
  #   X & Y: matrices to compute class corl diffs 
  #          samples in common must be rows (same as compareCorls above)
  #   trueCorlDiffs: matrix of class corl diffs for each miRNA-pathway pair
  #   sampClasses: classes that each sample belongs to
  #   resamplings: number of resamplings to construct null distribution
  # Returns:
  #   matrix of pvalues for each miRNA-pathway pair
  largerRandDiffCount <- 0 * trueCorlDiffs
  for (i in 1:resamplings) {
    randClasses <- sample(sampClasses) # scramble sample classes
    randCorlDiffs <- compareCorls(X, Y, randClasses, classA = classA, classB = classB)
    largerRandDiffs <- abs(randCorlDiffs) >= abs(trueCorlDiffs)
    largerRandDiffs[is.na(largerRandDiffs)] <- 0  #change NA's to 0
    largerRandDiffCount <- largerRandDiffCount + largerRandDiffs
  }
  return((largerRandDiffCount + 1)/(resamplings + 1))
}

#====================================================================================
# wrapper function to compute miRNA-PAS class-corl diffs
#====================================================================================

computemiRNAPSSCorls <- function(miRNA, PSS, resamplings) {
  #wrapper function to compute miRNA-pathway class corl diffs
  # Args:
  #   miRNA: miRNA gene expression matrix (samples by miRNAs)
  #   PSS: matrix of PSS (samples by pathway 1st coordinate PSS)
  #   resamplings: how many permutations to perform
  # Returns:
  #   A list of miRNA-pathway corl diffs, including:
  #     trueCorlDiffs: the actual class-corl diffs
  #     pvals: pval of the class-corl diffs
  #     commonSamples: vector of common samples 
  #get common samples between matrices
  commonData <- findCommonSamples(miRNA, PSS)
  classes <- substr(commonData$commonSamples, 14, 15)
  trueCorlDiffs <- with(commonData, compareCorls(mirna, paths, classes))
  pvals <- with(commonData, computePvals(mirna, paths, trueCorlDiffs, 
                                         classes, resamplings = resamplings))
  return(list(trueCorlDiffs = trueCorlDiffs, pvals = pvals, 
              commonSamples = commonData$commonSamples))
}

#====================================================================================
# finds most significant miRNA-pathway pairs by corl diff
#====================================================================================

sortPairs <- function(corllist, significance = NULL) {
  #finds most significant miRNA-pathway pairs
  # Args:
  #   corllist: output of computemiRNAPSSCorls (list including 
  #               trueCorlDiffs, pvals, commonSamples)
  #   significance: max pval to include
  # Returns:
  #   df of miRNA-pathway pairs sorted by corl diff 
  #   (only signifiance miRNA-pathway pairs reported)
  if (is.null(significance)) {
    corls <- with(corllist, trueCorlDiffs[pvals <= min(pvals)])
  } else {
    corls <- with(corllist, trueCorlDiffs[pvals <= significance])
  }
  corlsorts <- sort(abs(corls), decreasing = TRUE)
  indices <- sapply(corlsorts, function(X) {
    vect <- with(corllist, which(abs(trueCorlDiffs) == X, arr.ind = TRUE))
    vect[1, ]
  })
  df <- with(corllist, 
             data.frame(miRNA = rownames(pvals)[indices["row", ]], 
                        pathway = colnames(pvals)[indices["col", ]], 
                        corldiff = corlsorts))
  return(df)
}

#====================================================================================
# pathway intersecting genes
#====================================================================================

PathwayIntersectGenes <- function(assayedGenes = NULL) {
  #Args:
  #   assayedGenes: default NULL, if given a vector of genes then computes 
  #                 overlap using only genes that have been assayed!
  #Returns: 
  #   matrix: ijth entry is number of genes pathways i and j have in common,
  #           normalized by the size of the smallest pathway in that pair 
  require(Biobase)
  require(org.Hs.eg.db)
  require(annotate)
  require(KEGG.db)
  allPaths <- getPathwayGenes()
  mat <- matrix(0, length(allPaths), length(allPaths), 
                dimnames = list(names(allPaths), names(allPaths)))
  if (!is.null(assayedGenes)) {
    allPaths <- sapply(allPaths, function(X) intersect(X, assayedGenes))
  }
  path.len <- sapply(allPaths, length)
  path.combs <- combn(names(allPaths), 2)
  stat <- apply(path.combs, 2, function(X) {
    genes <- intersect(allPaths[[X[1]]], allPaths[[X[2]]])
    minpath <- min(path.len[X[1]], path.len[X[2]])
    return(length(genes) / minpath)
  })
  for (i in 1:ncol(path.combs)) {
    mat[path.combs[1, i], path.combs[2, i]] <- stat[i]
  }
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  diag(mat) <- 1
  return(mat)
}

#==============================================================================
# pathway enrichment of mirna targets
#==============================================================================

stringMatch <- function(mirnaID, humanMirna) {
  # Args:
  #  mirnaID: mirna ID string to match
  #  humanMirna: vector of all humanMirna IDs from TargetScan
  # Returns:
  #   string that matches within humanMirna set
  #   matches as close as possible to input mirnaID
  mirnaID <- gsub("mir", "miR", mirnaID)
  if (mirnaID %in% humanMirna) return(mirnaID)
  index <- grep(mirnaID, humanMirna)
  if (length(index) == 1) {
    return(humanMirna[index])
  } else if (length(index) > 1) {
    mirnas <- humanMirna[index]
    return(mirnas[grep(paste0(mirnaID, "-"), mirnas)])
  } else if (length(index) == 0) {
    mirnaID <- gsub("(-)[1-2]$", "", mirnaID)
    if (mirnaID %in% humanMirna) return(mirnaID)
    mirnaID <- gsub("[a-z]$", "", mirnaID)
    if (mirnaID %in% humanMirna) return(mirnaID)
  } else {
    return(NA)
  }
}

computePathwayEnrichment <- function(mirnaID, pathwayID, KEGGgenes) {
  # runs hypergeometric test to see if pathway is enriched with miRNA gene targets
  # population is all KEGGgenes that have been assayed
  # successes are targets from TargetScan that have been assayed
  # Args:
  #  mirnaID: mirna ID string
  #  pathwayID: pathwayID string
  #  KEGGgenes: all unique KEGG genes that HAVE BEEN ASSAYED!
  # Returns:
  #   pvalue of pathway enrichment 
  #   (probability of pathway enrichment at x level or higher)
  require(org.Hs.eg.db)
  require(annotate)
  require(targetscan.Hs.eg.db)
  humanMirnaIdx <- grep("hsa", mappedkeys(targetscan.Hs.egMIRNA))
  humanMirna <- mappedkeys(targetscan.Hs.egMIRNA)[humanMirnaIdx]
  mirnas <- stringMatch(mirnaID = mirnaID, humanMirna = humanMirna)
  if (all(is.na(mirnas))) {
    warning(paste(mirnas, "not found in TargetScan"))
    return(NA)
  }
  mirTargets <- getmirTargets(mirnas)
  pathGenes <- intersect(unlist(getPathwayGenes(as.character(pathwayID))), 
                         KEGGgenes)
  x <- length(intersect(mirTargets, pathGenes))
  m <- length(intersect(mirTargets, KEGGgenes))
  n <- length(setdiff(KEGGgenes, mirTargets))
  k <- length(pathGenes)
  return(phyper(q = x - 1, m = m, n = n, k = k, lower.tail = FALSE))
}

#==============================================================================
# get Isomap spectra
#==============================================================================

getSpectra <- function(Iso, plot.spectra = TRUE) {
  # Args:
  #   Iso: Isomap output from vegan package
  #   plot.spectra: plot the spectra? TRUE or FALSE
  # Returns:
  #   Isomap spectra by total proportion
  spect <- Iso$eig/sum(Iso$eig)
  if (plot.spectra) {
    plot(spect, xlab = "Index", ylab = "spectra", main = "Isomap spectra")
  }
  return(spect)
}