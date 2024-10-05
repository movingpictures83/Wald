library("derfinder")
library("derfinderData")
library("GenomicRanges")
library("knitr")
library("DESeq2")

input <- function(inputfile) {
#####################################################################################
# REGION MATRIX
pheno <<- subset(brainspanPheno, structure_acronym == "AMY")
#p <- pheno[, -which(colnames(pheno) %in% c(
#    "structure_acronym",
#    "structure_name", "file"
#))]
#rownames(p) <<- NULL
regionMat <<- readRDS(inputfile)
#regionMat$chr21$regions
#head(regionMat$chr21$coverageMatrix)
#####################################################################################
}

run <- function() {}

output <- function(outputfile) {

#####################################################################################
# WALD
counts <- round(regionMat$chr21$coverageMatrix)
dse <- DESeqDataSetFromMatrix(counts, pheno, ~ group + gender)
dse <- DESeq(dse, test = "Wald", fitType = "local")
deseq <- regionMat$chr21$regions
mcols(deseq) <- c(mcols(deseq), results(dse))
write.csv(table(deseq$padj < 0.5), outputfile)
#####################################################################################

}
