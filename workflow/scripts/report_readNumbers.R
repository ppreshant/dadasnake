log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

sampleFile <- snakemake@input[[1]]
filesOI <- unique(unlist(snakemake@input)[-1])

print(paste0("reading sample table from ",sampleFile))
sampleTab <- read.delim(sampleFile,stringsAsFactors=F)


# count raw reads
if(snakemake@params[["currentStep"]] == "raw"){
  print("extracting read numbers")
 # print(system2("which",args=c("Rscript"),stdout=T))

  readnums <- sapply(filesOI,function(x){ # Count # of reads in each file
    if(grepl(".gz$",x)){
      as.numeric(system2("zcat",args=c(x,"| wc -l"),stdout=T))/4 # for compressed files (.gz)
    }else{
      as.numeric(unlist(strsplit(system2("wc",args=c("-l",x),stdout=T),split=" "))[1])/4
    }
  })

  # ensure that the prefix directory ends with a single '/'
  prefix <- gsub("[/]{2}$","/",paste0(snakemake@params[["raw_directory"]],"/"))

  # get the raw file names (which should match to the r1_file, r2_file columns in the sample_table.tsv)
  names(readnums) <- gsub(paste0("^[/]{0,1}[^/]*",prefix),"",names(readnums)) # remove the directory path prefix and anything before

  # write the respective read numbers in the table, extracting from the readnums list(/vector?)
  sampleTab$reads_raw_r1 <- sapply(sampleTab$r1_file,function(x) readnums[x])
  sampleTab$reads_raw_r2 <- sapply(sampleTab$r2_file,function(x) readnums[x])
  write.table(sampleTab,snakemake@output[[1]],sep="\t",quote=F,row.names=F)


  # count after primer removal
}else if(snakemake@params[["currentStep"]] == "primers"){
  print("extracting read numbers")
  readnums <- sapply(filesOI,function(x) {
    if(grepl(".gz$",x)){
      as.numeric(system2("zcat",args=c(x,"| wc -l"),stdout=T))/4
    }else{
      as.numeric(unlist(strsplit(system2("wc",args=c("-l",x),stdout=T),split=" "))[1])/4
    }
  })
  prefix <- "preprocessing/"
  suffix <- ".{5}fastq[.gz]*"
  names(readnums) <- gsub(prefix,"",names(readnums))
  runs_libs <- sapply(names(readnums),function(x) unlist(strsplit(x,split="/")))
  runs_libs[2,] <- gsub(suffix,"",runs_libs[2,])
  print(runs_libs)
  sampleTab$reads_primers_fwd <- apply(sampleTab[,c("run","library")],1,
                                       function(x) readnums[which(runs_libs[1,]==x[1]
                                                                  &runs_libs[2,]==x[2]
                                                                  &grepl("fwd.fastq",names(readnums)))])
  sampleTab$reads_primers_rvs <- apply(sampleTab[,c("run","library")],1,
                                       function(x) readnums[which(runs_libs[1,]==x[1]
                                                                  &runs_libs[2,]==x[2]
                                                                  &grepl("rvs.fastq",names(readnums)))])
  write.table(sampleTab,snakemake@output[[1]],sep="\t",quote=F,row.names=F)

  tabPerSample <- aggregate(sampleTab[,c("library","run","r1_file","r2_file")],list(sampleTab$sample),function(x)paste(unique(x),collapse=",",sep=","))
  colnames(tabPerSample)[1] <- "sample"
  numsPerSample <- aggregate(sampleTab[,c("reads_raw_r1","reads_raw_r2","reads_primers_fwd","reads_primers_rvs")],list(sampleTab$sample),sum)
  colnames(numsPerSample)[1] <- "sample"
  perSample <- merge(tabPerSample,numsPerSample,by="sample")
  write.table(perSample,snakemake@output[[2]],sep="\t",quote=F,row.names=F)



  # count reads after filtering
}else if(snakemake@params[["currentStep"]] == "filtered"){
  print("extracting read numbers")
  readnums <- sapply(filesOI,function(x) as.numeric(unlist(strsplit(system2("zcat",args=c(x,"| wc -l"),stdout=T),split=" "))[1])/4)
  prefix <- "filtered/"
  suffix <- ".{5}fastq.gz"
  names(readnums) <- gsub(prefix,"",names(readnums))
  runs_sams <- sapply(names(readnums),function(x) unlist(strsplit(x,split="/")))
  runs_sams[2,] <- gsub(suffix,"",runs_sams[2,])
  sampleTab$reads_filtered_fwd <- apply(sampleTab[,c("run","sample")],1,
                                       function(x) readnums[which(runs_sams[1,]==x[1]
                                                                  &runs_sams[2,]==x[2]
                                                                  &grepl("fwd.fastq",names(readnums)))])
  sampleTab$reads_filtered_rvs <- apply(sampleTab[,c("run","sample")],1,
                                       function(x) readnums[which(runs_sams[1,]==x[1]
                                                                  &runs_sams[2,]==x[2]
                                                                  &grepl("rvs.fastq",names(readnums)))])
  write.table(sampleTab,snakemake@output[[1]],sep="\t",quote=F,row.names=F)
  tabPerSample <- aggregate(sampleTab[,c("library","run","r1_file","r2_file")],list(sampleTab$sample),function(x)paste(unique(x),collapse=",",sep=","))
  colnames(tabPerSample)[1] <- "sample"
  numsPerSample <- aggregate(sampleTab[,c("reads_raw_r1","reads_raw_r2","reads_primers_fwd","reads_primers_rvs")],list(sampleTab$sample),sum)
  colnames(numsPerSample)[1] <- "sample"
  tab2PerSample <- merge(tabPerSample,numsPerSample,by="sample")
  nums2PerSample <- aggregate(sampleTab[,c("reads_filtered_fwd","reads_filtered_rvs")],list(sampleTab$sample),function(x) sum(unique(x)))
  colnames(nums2PerSample)[1] <- "sample"
  perSample <- merge(tab2PerSample,nums2PerSample,by="sample")
  write.table(perSample,snakemake@output[[2]],sep="\t",quote=F,row.names=F)
}else if(snakemake@params[["currentStep"]] == "merged"){
  if(!require(dada2)){
    BiocManager::install("GenomeInfoDbData",update=F,ask=F)
    require(dada2)
  }
  print("extracting read numbers")
  getN <- function(x) sum(getUniques(x))
  if(length(filesOI)>1 | !snakemake@params[['pool']] %in% c("true","pseudo")){
    readnums <- sapply(filesOI,function(x){
      if(file.info(x)$size>0) getN(readRDS(x)) else 0
    })
    runs_sams <- sapply(gsub(".RDS$","",names(readnums)),function(x) unlist(strsplit(x,split="/")))
    sampleTab$reads_merged <- apply(sampleTab[,c("run","sample")],1,
                                  function(x) readnums[which(runs_sams[2,]==x[1]
                                                             &runs_sams[3,]==x[2])])
  }else{
    if(file.info(filesOI)$size>0) readnums <- getN(readRDS(filesOI)) else readnums <- 0
    names(readnums) <- filesOI
    runs_sams <- sapply(gsub(".RDS$","",names(readnums)),function(x) unlist(strsplit(x,split="/")))
    sampleTab$reads_merged <- apply(sampleTab[,c("run","sample")],1,
                                  function(x) readnums[which(runs_sams[1,]==x[1]
                                                             &runs_sams[2,]==x[2])])
  }
  write.table(sampleTab,snakemake@output[[1]],sep="\t",quote=F,row.names=F)
  tabPerSample <- aggregate(sampleTab[,c("library","run","r1_file","r2_file")],list(sampleTab$sample),function(x)paste(unique(x),collapse=",",sep=","))
  colnames(tabPerSample)[1] <- "sample"
  numsPerSample <- aggregate(sampleTab[,c("reads_raw_r1","reads_raw_r2","reads_primers_fwd","reads_primers_rvs")],list(sampleTab$sample),sum)
  colnames(numsPerSample)[1] <- "sample"
  tab2PerSample <- merge(tabPerSample,numsPerSample,by="sample")
  nums2PerSample <- aggregate(sampleTab[,c("reads_filtered_fwd","reads_filtered_rvs","reads_merged")],list(sampleTab$sample),function(x)sum(unique(x)))
  colnames(nums2PerSample)[1] <- "sample"
  perSample <- merge(tab2PerSample,nums2PerSample,by="sample")
  write.table(perSample,snakemake@output[[2]],sep="\t",quote=F,row.names=F)
}else if(snakemake@params[["currentStep"]] == "table"){
  print("extracting read numbers")
  if(file.info(filesOI)$size>0){
    readnums <- rowSums(readRDS(filesOI))
    if(length(readnums)==1 & is.null(names(readnums))) names(readnums) <- sampleTab$sample
    sampleTab$reads_tabled <- sapply(sampleTab$sample,
                                       function(x) if(x %in% names(readnums)) readnums[names(readnums)==x] else 0)
  }else{
    sampleTab$reads_tabled <- 0
  }
  write.table(sampleTab,snakemake@output[[1]],sep="\t",quote=F,row.names=F)
}else if(snakemake@params[["currentStep"]] == "chimera"){
  print("extracting read numbers")
  if(file.info(filesOI[1])$size>0){
    readnums <- rowSums(readRDS(filesOI[1]))
    if(length(readnums)==1 & is.null(names(readnums))) names(readnums) <- sampleTab$sample
    sampleTab$reads_tabled <- sapply(sampleTab$sample,
                                       function(x) if(x %in% names(readnums)) readnums[names(readnums)==x] else 0)
    readnums <- rowSums(readRDS(filesOI[2]))
    if(length(readnums)==1 & is.null(names(readnums))) names(readnums) <- sampleTab$sample
    sampleTab$reads_chimera_checked <- sapply(sampleTab$sample,
                                       function(x) if(x %in% names(readnums)) readnums[names(readnums)==x] else 0)
  }else{
    sampleTab$reads_tabled <- 0
    sampleTab$reads_chimera_checked <- 0
  }
  write.table(sampleTab,snakemake@output[[1]],sep="\t",quote=F,row.names=F)
}else if(snakemake@params[["currentStep"]] == "post"){
  print("extracting read numbers")
  if(file.info(filesOI[1])$size>0){
    tmpOTU <- readRDS(filesOI[1])
    if(length(sampleTab$sample)>1){
      readnums <- colSums(tmpOTU[,colnames(tmpOTU) %in% ifelse(grepl("^[[:digit:]]",sampleTab$sample),paste0("X",sampleTab$sample),sampleTab$sample)])
    }else{
      if(!colnames(tmpOTU)[3] %in% sampleTab$sample) colnames(tmpOTU)[3]<- sampleTab$sample
      readnums <- sum(tmpOTU[,colnames(tmpOTU) %in% sampleTab$sample])
      names(readnums) <- sampleTab$sample
    }
    sampleTab$reads_tax.length_filtered <- sapply(ifelse(grepl("^[[:digit:]]",sampleTab$sample),
                                                         paste0("X",sampleTab$sample),
                                                         sampleTab$sample),
                                       function(x) if(x %in% names(readnums)) readnums[names(readnums)==x] else 0)
  }else{
    sampleTab$reads_tax.length_filtered <- 0
  }
  write.table(sampleTab,snakemake@output[[1]],sep="\t",quote=F,row.names=F)
}
