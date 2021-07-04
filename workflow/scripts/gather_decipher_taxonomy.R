log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

print("reading DECIPHER result")
for(dt in 1:length(snakemake@input)){
  currName <- gsub(".RDS$","",gsub("sequenceTables/tax.decipher.","",snakemake@input[[dt]]))
  decTax <- readRDS(snakemake@input[[dt]])
  decTax <- decTax[,-ncol(decTax)]
  colnames(decTax) <- paste0(colnames(decTax),".decipher.",currName)
  decTax$OTU <- rownames(decTax)
  if(!"lastTab" %in% ls()) lastTab <- decTax else lastTab <- merge(lastTab,decTax,by="OTU", all=T)
  lastTab[is.na(lastTab)] <- ""
}
print("Saving table with all taxonomies")
saveRDS(lastTab,snakemake@output[[1]])

