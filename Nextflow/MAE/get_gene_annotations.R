#' MAE gene annotations
#' Last update: 24-04-2025
#' Argument 1= input path .tsv
#' Argument 2= output path .tsv
#' Argument 3= optional genome build (hg19 or hg38)


library(dplyr)
library(tidyr)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("data.table")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)


args <- commandArgs(trailingOnly = TRUE)

genome_build <- args[3]

# Annotate chr start end.
if (genome_build == "hg19") {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
} else {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
}

res <- read.table(args[1], header=TRUE, sep="\t")

gr <- GRanges(
  seqnames = paste0("chr",res$chr),
  ranges = IRanges(start = res$start, end = res$end)
)

genes <- genes(txdb)

hits <- findOverlaps(gr, genes)

res$entrez <- NA
res$entrez[queryHits(hits)] <- names(genes)[subjectHits(hits)]

res$hgncSymbol = mapIds(org.Hs.eg.db,
                    keys=res$entrez, 
                    column="SYMBOL",
                    keytype="ENTREZID",
                    multiVals="first")

res <- as.data.frame(apply(res,2,as.character))

write.table(res, args[2], sep='\t', append = FALSE, row.names = FALSE, col.names = TRUE)