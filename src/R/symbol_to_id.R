home.folder <- "/home/bit/lacerda/data/STRING/"
home.folder <- "C:/Users/Mauricio/Thesis/data/STRING/"

res <- read.csv(file=paste(home.folder, "string_symbol.edgelist", sep=""), header = FALSE, sep = "\t")
colnames(res) = c('Protein1', 'Protein2', 'assoc')

genes <- unique(append(as.vector(res$Protein1), as.vector(res$Protein2)))

library("AnnotationDbi")
library("org.Hs.eg.db")

symbolToEntrez <- mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
symbolToEntrez <- symbolToEntrez[which(!is.na(symbolToEntrez))]

ensgToEntrez <- mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
ensgToEntrez <- ensgToEntrez[which(!is.na(ensgToEntrez))]

resDF <- append(symbolToEntrez, ensgToEntrez)

write.table(resDF, file=paste(home.folder, "symbol_to_entrez.txt", sep=""), quote = FALSE, col.names = FALSE, sep="\t", dec=".")
