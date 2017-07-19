# Creates data objects from raw text files downloaded from TCGA

# change to project repo if necessary
source("source/CombineData.R")

#########################################################################################
# data objects
#########################################################################################

# miRNAHiSeq dataframe 
BRCA_miRNASeq <- CombinemiRNASeqData(file = "rawdata/BRCA/Platform Data/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3")
save(BRCA_miRNASeq, file = "data/BRCA_miRNASeq.Rda")

# RNASeq dataframe (scaled_estimate from RSEM)
BRCA_RNASeq <- CombineRNASeqData(directory = "rawdata/BRCA/Platform Data/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3",
                              output = "*.genes.results", scaled = TRUE, 
                              filemap = "rawdata/BRCA/Platform Data/RNASeqV2/FILE_SAMPLE_MAP.txt")
save(BRCA_RNASeq, file = "data/BRCA_RNASeq.Rda")

# ExpressionSet for RNASeq dataframe
BRCA_Eset <- CreateEset(path = "data/BRCA_RNASeq.RData", 
                      clinical.file = "rawdata/BRCA/Clinical Data/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt",
                      type = "RNA", cancer = "BRCA")
save(BRCA_Eset, file = "data/BRCA_Eset.Rda")



