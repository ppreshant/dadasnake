# script to pick the top blast hit for each query (OTU)
# Author : Prashant K.
# Date: Jan 2023

# logger ----

# logfile code : copied from report_readNumbers.R <- cutadapt.smk
log <- file(snakemake@log[[1]], open="wt") # this is causing syntax error due to "(".. hmm
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))


# load package
library(tidyverse) # using snakemake loaded env : workflow/envs/tidyverse_env.yaml


# load file ----

# load the the counts with taxonomy data : Grab from snakemake input
taxonomy_w_summary <- read_tsv(snakemake@input[[1]])


# collate by organism ----

# collate the matching organism counts (combining ASV identities)
collated_tax_counts <- 
  select(taxonomy_w_summary, -evalue, -alignment_len, -asv_len, -Row.names) %>% 
  replace_na(list(sseqid = 'NA')) %>% 
  
  # collate stuff
  reframe(
    
    mismatch_range = str_c(min(mismatch), max(mismatch), sep = '-'), # show mismatch range
    
    across(where(is.numeric), sum), # sum all the counts (incl. mismatches)
    ASV_count = n(), # count num of ASVs combined
    ASVs = str_replace(ASV, 'ASV_', '') %>% str_c(collapse = ', '), # combine all ASVs (shortened)
    
    .by = sseqid) %>% 
  
  # format properly
  relocate(mismatch_range, mismatch, .after = ASV_count) # move columns to later

# Save the collated counts data by organism
write_tsv(collated_tax_counts, snakemake@output[['collated_summary']])

# plotting ----  

organism_shortener <- 
  c('Shewanella_oneidensis_MR1_16S' = 'S.onei',
    'Vibrio_natriegens_14048_16S' = 'V.nat',
    'Acinetobacter_baylyi_ADP1_16S' = 'A.bayl',
    'Pseudomonas_putida_F1_16S' = 'P.puti',
    'E_coli_MG1655_16S' = 'E.coli',
    'Sinorhizobium_meliloti_1021_16S' = 'S.meli',
    'NA' = 'no match'
  )

# prepare for plotting: shorten organism names, reorder, pivot
collated_clean <- 
  select(collated_tax_counts, -ASV_count, -contains('mismatch'), -ASVs) %>% # select columns for plot
  
  # shorten organism name, reorder
  mutate(sseqid = str_replace_all(sseqid, organism_shortener) %>%  fct_relevel(organism_shortener)) %>% 
  arrange(sseqid) %>% 
  
  # bring all samples and counts into a column each
  pivot_longer(cols = -sseqid, names_to = 'sample', values_to = 'counts')

# to make a simple ggplot, facet by sample. y-axis with species 
# facet order: (all g's in a row); all m's in next row to make easy comparisons

plt1 <- 
  ggplot(collated_clean,
          aes(x = counts, y = sseqid)) + 
      geom_point() +
      
      labs(y = NULL) + # remove y-axis label
      scale_y_discrete(limits = rev) + # reverse y-axis order (more intuitive)
      facet_wrap(facets = vars(sample), nrow = 2)

# save plot
ggsave(snakemake@output[['plot_linear']], plt1)

# plot in logscale
plt1log <- 
  plt1 + 
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

# save plot
ggsave(snakemake@output[['plot_log']], plt1log)
