# Preprocesses ExpressionSet & normalizes data
# Removes samples that are NOT either tumor (sample code "01")
# or adjacent-normal (sample code "11") prior to normalization 

# change to project repo if necessary
load("data/BRCA_Eset.Rda")
source("source/AnalysisFunctions.R")
library(limma)

#==============================================================================
# normalization
#==============================================================================
                               
#only tumor & normal
BRCAEsetTPM <- BRCA_Eset[, substr(sampleNames(BRCA_Eset), 14, 15) %in% c("01", "11")]

#normalize libraries so all sum to 1
exprs(BRCAEsetTPM) <- apply(exprs(BRCAEsetTPM), 2, function(X) X / sum(X))

#remove low median expression genes
BRCAEsetTPM <- BRCAEsetTPM[apply(exprs(BRCAEsetTPM), 1, median) >= 1e-9, ]

#log2 transform plus a small offset
exprs(BRCAEsetTPM) <- log2((exprs(BRCAEsetTPM) + 1e-10) * 1e6)

#save normalized data in TPM
save(BRCAEsetTPM, file = "data/BRCAEsetTPM.RData")
