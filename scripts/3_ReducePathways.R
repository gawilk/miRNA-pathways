# Subsets KEGG pathways to RNASeq data
# Then applies Isomap to each KEGG pathway to get Isomap embedding 

# change to project repo if necessary
load("data/BRCAEsetTPM.RData")
source("source/AnalysisFunctions.R")
library(vegan)
library(parallel)

#==============================================================================
# subset by pathways
#==============================================================================

# subset by pathways
# scale gene vectors to unit variance
GEpaths <- subsetGEbyKEGG(BRCAEsetTPM)
GEpaths <- lapply(GEpaths, scale)
GEpaths <- GEpaths[sapply(GEpaths, ncol) > 5]

#==============================================================================
# ISOMAP
#==============================================================================

BRCAPathIso <- mclapply(GEpaths, function(pathway) {
  knn <- findKisomap(pathway, K = seq(4,20), scale = FALSE, plot = FALSE)
  pathway <- scale(pathway)
  distX <- as.matrix(dist(pathway))
  iso <- isomap(distX, ndim = 6, k = knn$k)
  iso <- c(iso, k = knn$k)
  return(iso)
}, mc.cores = 6, mc.preschedule = FALSE)
save(BRCAPathIso, file = "data/BRCAPathIso.RData")

