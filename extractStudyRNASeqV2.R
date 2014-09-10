#!/usr/bin/Rscript --vanilla --slave

args <- commandArgs(trailingOnly=T)
study <- args[1]

#hard code data location
setwd("/home/rj67/group/RNASeqV2")

load(file="candidate_gene_list.RData")

manifest <- read.delim( paste(study, '/FILE_SAMPLE_MAP.txt', sep=""), header =T, strip.white=T, stringsAsFactors = F)
colnames(manifest)[2] <- "barcode"
manifest <- manifest[grep('rsem.genes.normalized_results', manifest$filename),]

file_paths <- list.files(path =paste(study, "/RNASeqV2", sep=""), pattern = "rsem.genes.normalized_results", recursive=T, full.names=T)

all_mrna <- do.call(rbind, lapply( file_paths, function(sample_file){
  df <- read.delim(sample_file, header =T, strip.white = T, stringsAsFactors = F );
 
  #CASP12|120329-> 100506742
  df$gene_id[df$gene_id=="CASP12|120329"] <- "CASP12|100506742" 
  df$gene_id[df$gene_id=="DUX4|22947"] <- "DUX4|100288687" 
  #lose C15orf65,  NUTM2B, SLX1A
  df$Entrez<-sapply(df$gene_id, function(x) strsplit(x, split="|", fixed=T)[[1]][2])
  df <- merge(df, list_goi[c("Entrez.Gene", "Gene")], by.x="Entrez", by.y="Entrez.Gene")
  df$Entrez <- NULL
  df$gene_id <- NULL
 
  path_split <- strsplit(sample_file, split='/', fixed=F)[[1]] ;
  df$study <- study
  df$platform <- path_split[which(path_split=="RNASeqV2")+1]
  df$filename <- rev(path_split)[1]
  return(df)
}))

all_mrna <- merge(all_mrna, manifest, by ="filename")
all_mrna$filename <- NULL
all_mrna$Patient <- substr(all_mrna$barcode, 6, 12)
all_mrna$sample_type <- substr(all_mrna$barcode, 14, 15)
all_mrna$event_uid <- apply(all_mrna[c("Patient", "Gene")], 1, function(x) paste0(x, collapse="-") )

file_out <- paste(study, "_mrna.RData", sep="")
assign(paste(study, "_mrna", sep=""), all_mrna)
save(list=paste(study, "_mrna", sep=""), file=file_out)

