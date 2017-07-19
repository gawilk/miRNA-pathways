# Correlates miRNAs and Isomap-reduced pathways class-conditionally
# By either tumor or normal tissue

# change to project repo if necessary
load("data/BRCA_miRNASeq.Rda")
load("data/BRCAPathIso.RData")
source("source/AnalysisFunctions.R")

#==============================================================================
# Some miRNA preprocessing
#==============================================================================

# flip matrix (now samples by miRNAs) and filter out lowly expressed miRNAs
miRNA <- t(BRCA_miRNASeq) 
miRNA <- miRNA[, apply(miRNA, 2, function(X) sum(X > 1) > (0.5 * length(X)))]

#==============================================================================
# Correlate with Isomap-reduced pathways
#==============================================================================

# get 1st coordinate PAS
ISO <- sapply(BRCAPathIso, function(pathway) pathway$points[, 1])

# run correlation and permutation tests
BRCACorlDiffs1e5 <- computemiRNAPSSCorls(miRNA = miRNA, 
	PSS = ISO, resamplings = 1e5)
save(BRCACorlDiffs1e5, file = "data/BRCACorlDiffs1e5.RData")
