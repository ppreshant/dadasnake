log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(BiocParallel)
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}else{
    parallel <- FALSE
    register(SerialParam())
}
library(dada2)

# File parsing
errFfile <- snakemake@input[[1]]
errRfile <- snakemake@input[[2]]

filtF <- snakemake@input[[3]]
filtR <- snakemake@input[[4]]

sampleNameF <- gsub("/","",gsub("filtered/","run.",gsub(".fwd.fastq.gz","",filtF)))
sampleNameR <- gsub("/","",gsub("filtered/","run.",gsub(".rvs.fastq.gz","",filtR)))

if(sampleNameF!=sampleNameR) stop("Forward and reverse files do not match.")

sampleName <- gsub(".+/","",gsub(".fwd.fastq.gz","",filtF))

names(filtF) <- sampleName
names(filtR) <- sampleName

mergefile <- snakemake@output[[1]]

print("merging")
errF <- readRDS(errFfile)
errR <- readRDS(errRfile)

# Sample inference and merger of paired-end reads
print(paste0("make dada object and merge, ",sampleName))
derepF <- derepFastq(filtF)
dadaF <- dada(derepF, err=errF, multithread=snakemake@threads)
derepR <- derepFastq(filtR)
dadaR <- dada(derepR, err=errR, multithread=snakemake@threads)
merger <- mergePairs(dadaF, derepF, dadaR, derepR,
                     minOverlap=snakemake@config[['pair_merging']][['min_overlap']],
                     maxMismatch=snakemake@config[['pair_merging']][['max_mismatch']],
                     justConcatenate=snakemake@config[['pair_merging']][['just_concatenate']],
                     trimOverhang=snakemake@config[['pair_merging']][['trim_overhang']])
saveRDS(merger,mergefile)

print("done")