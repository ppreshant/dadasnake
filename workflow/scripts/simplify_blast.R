# script to pick the top blast hit for each query (OTU)
# Author : Prashant K.
# Date: Jan 2023

# logfile code : copied from report_readNumbers.R <- cutadapt.smk
log <- file(snakemake@log[[1]], open="wt") # this is causing syntax error due to "(".. hmm
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))


# load package
library(tidyverse) # using snakemake loaded env : workflow/envs/tidyverse_env.yaml

# name the blast result file : Grab from snakemake input
blast_result_file <- snakemake@input[[1]] # 'all.seqs.blast.txt' for running through snakemake rule

# load blast result file ----
.dat <- read_tsv(blast_result_file,
  col_names = c('qseqid', 'sseqid', 'evalue', 'mismatch', 'length'))

# rearrange, retain only top hit per OTU
top_hit <- .dat %>%
  nest_by(qseqid) %>% # make one row per query/OTU
  summarize(top_hit = head(data, 1)) %>% # retain only the top hit for each query/row
  unnest(cols = top_hit) # expand the columns of the nesting for top_hit

# Save cleaned result
write_tsv(top_hit, snakemake@output[[1]])
