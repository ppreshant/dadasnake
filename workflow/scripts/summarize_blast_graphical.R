# adhoc analysis of BLAST top hits / OTUs
# TODO : need to make snakemake compatible as in simplify_blast.R ;
# save pngs into analysis/ folder within mapping/

# load package
library(tidyverse) # in snakemake loaded env : workflow/envs/tidyverse_env.yaml

# Load data
.topdat <- read_tsv("top_all.seqs.blast.txt")
.seqtab <- read_tsv("../sequenceTables/all.seqTab.tax.tsv", col_select = -`Row.names`)

# Processing : join with taxonomy and sample reads ; order for plotting
procd <- mutate(.topdat, sseqid = fct_infreq(sseqid) %>% fct_rev) %>% # make factor order for plotting
  left_join(.seqtab, by = c('qseqid' = 'OTU'))

write_tsv(procd, "top_seqs_w_tax.csv", na = '') # save data for visual inspection


# barplot for counts / OTU
ggplot(procd, aes(y = sseqid)) + geom_bar()

# jitterplot for read-lengths
ggplot(.dat, aes(length, sseqid)) + geom_jitter(height = .2, alpha = 0.2)

# jitterplot for mismatches
ggplot(procd, aes(mismatch, sseqid)) + geom_jitter(height = .2, alpha = 0.2)


# summarize read counts per BLAST match organism
summd <- group_by(procd, sseqid) %>% summarize(across(contains("16S"), sum))
write_tsv(summd, "reads_by_matched_org.csv") # save data  
