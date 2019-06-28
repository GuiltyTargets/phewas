
#############################################
############### Load the data ###############
#############################################
home.folder <- "data/MSBB/"

metaData <- read.csv(file=paste(home.folder, "MSBB_clinical.csv", sep=""))
countData <- read.csv(file=paste(home.folder, "MSSM_all_counts_matrix.txt", sep=""),  sep = "\t")
covariates <- read.csv(file=paste(home.folder, "MSBB_RNAseq_covariates.csv", sep=""))

library(plyr)
colData <- join(metaData, covariates, by="individualIdentifier", type="inner")
colData <- colData[colData$BrodmannArea == 'BM22', ]

#############################################
############# Prepare the data ##############
#############################################
rownames(countData) <- gsub('\\..+$', '', countData$feature)
countData <- countData[5:length(rownames(countData)) ,2:length(colnames(countData))]

colData <- colData[complete.cases(colData[ , "NP.1"]),]
colData <- colData[!duplicated(colData$sampleIdentifier), ]

rownames(colData) <- colData$sampleIdentifier
colData <- colData[rownames(colData) %in% colnames(countData), ]

colData$BrodmannArea <- factor(colData$BrodmannArea)
colData$AOD <- factor(colData$AOD)
colData$SEX <- factor(colData$SEX)
colData$RACE <- factor(colData$RACE)
colData$NP.1 <- factor(colData$NP.1)
colData$NP.1 <- relevel(colData$NP.1, ref="1")

countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

#############################################
########### Generate DESeq object ###########
#############################################
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData   = colData,
                              design    = ~ AOD + SEX + NP.1 + RACE)
# Filter out the rows that don't have information on differential gene expression
dds <- dds[rowSums(counts(dds)) > 1, ]

#############################################
###### Differential expression analysis #####
#############################################
dds <- DESeq(dds)
res <- results(dds, contrast=c("NP.1","1","2"))
res
summary(res)

#############################################
############# Write to the file #############
#############################################
library("AnnotationDbi")
library("org.Hs.eg.db")

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="list")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="list")

resDF <- as.data.frame(res
write.table(resDF, file=paste(home.folder, "DifferentialExpression_BM22.csv", sep=""), sep=";", dec=".")

