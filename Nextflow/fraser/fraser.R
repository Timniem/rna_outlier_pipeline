#' FRASER autoencoder 
#' Gagneur-lab FRASER (2.0)
#' Processes start from a samplesheet with SampleID's BAM paths featurecount settings etc. 
#' and creates fraser rds object and results .tsv
#' 28-10-2023
#' Argument 1: Samplesheet
#' Argument 2= input/output folder

library(FRASER)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)

# Setup parallelisation
if(.Platform$OS.type == "unix") {
    register(MulticoreParam(workers=min(4, multicoreWorkers())))
} else {
    register(SnowParam(workers=min(4, multicoreWorkers())))
}

workdir <- args[2]

# Load original sample table

original_settingsTable <- fread(args[1])

fds <- loadFraserDataSet(dir=workdir)
fds <- filterExpressionAndVariability(fds, minDeltaPsi=0, filter=FALSE)
fds <- fds[mcols(fds, type="j")[,"passed"],]

# Hyperparam optim
set.seed(42)
fds <- optimHyperParams(fds, type="jaccard", plot=FALSE)
best_q <- bestQ(fds, type="jaccard")
fds <- FRASER(fds, q=c(jaccard=best_q))

# Using different method of annotation
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(org.Hs.eg.db)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#orgDb <- org.Hs.eg.db
#fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)

fds <- annotateRanges(fds) # previous method

fds <- calculatePadjValues(fds, type="jaccard", geneLevel=TRUE) # geneLevel TRUE -> FALSE

saveFraserDataSet(fds, dir=workdir, name="fraser_out")

register(SerialParam())
res <- as.data.table(results(fds,aggregate=TRUE, all=TRUE))

res <- res[res$sampleID %in% original_settingsTable$sampleID] #filter out samplesheet samples

# Rename for compatibility
names(res)[names(res) == 'seqnames'] <- 'chr'

write.table(res, paste("combined_samples", 'result_table_fraser.tsv', sep='_'), sep='\t', append = FALSE, row.names = FALSE, col.names = TRUE)

# Results per patient
for (sampleid in unique(res$sampleID)){
    sample_out_path = paste(sampleid, 'result_table_fraser.tsv', sep='_')
    write.table(res[res$sampleID == sampleid], sample_out_path, sep='\t', append = FALSE, row.names = FALSE, col.names = TRUE)
}
