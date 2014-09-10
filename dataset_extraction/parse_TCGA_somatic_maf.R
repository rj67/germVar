#!/usr/bin/Rscript --vanilla --slave

# get a list of maf file names, read them in, concatenate, output
fns <- list.files('June_23', pattern = ".maf", recursive = T, full.names = T)

#fns <- c("all_tumor/genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.2.1.0.maf")
#print(fns)
read_maf<-function(fn){
  print(fn)
  df <- read.delim(fn, header = T, strip.white = T, comment.char = "#", row.names=NULL, sep='\t' )
  
  # change colnames to lower case
  colnames(df) <- tolower(colnames(df))
  
  # some file use chromosome instead of chrom, change it
  colnames(df) <- gsub('chromosome', 'chrom', colnames(df))
  #df$Patient <- sapply(df$tumor_sample_barcode, function(x) gsub('-', '.', substr(x, 1, 12), fixed =T))
  df$Patient <- sapply(df$tumor_sample_barcode, function(x) substr(x, 6, 12))
  df = df[c('hugo_symbol', 'entrez_gene_id','center','ncbi_build','variant_classification','variant_type','Patient', 'chrom','start_position','end_position', 'reference_allele','tumor_seq_allele1','tumor_seq_allele2',
    'tumor_sample_barcode','matched_norm_sample_barcode', 'archive_name')]
  df$tumor_sample_barcode <- factor(df$tumor_sample_barcode)
  df$matched_norm_sample_barcode <- factor(df$matched_norm_sample_barcode)
  df$archive_name <- factor(df$archive_name)
  #colnames(df)[1:6] <- c('Gene', 'Entrez', 'Center', 'Build', 'EFF', 'VariantType')
  # filename as source 
  #df$source = strsplit(fn, split='/', fixed =F)[[1]][2]
  
  return(df)
}
df <- do.call(rbind, lapply(fns, read_maf))
print('dimension')
print( dim(df))
write.table(df, file ='June_23/TCGA_all_somatic.redundent.tsv', row.names= F, quote =F, sep = '\t')

# remove duplicates
df$uid <- apply(df[c("ncbi_build", "chrom", "reference_allele", "tumor_seq_allele1", "tumor_seq_allele2", "Patient")], 1, function(x) paste0(x, collapse="-"))
df <- df[!duplicated(df$uid), ]

print(table(df$ncbi_build))
#fn <- "all_tumor/genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.2.0.0.somatic.maf"

#rownames(df) <- 1:nrow(df)
#print.data.frame(df[1:5,])
#
print('dimension')
print( dim(df))
print('# patient')
print(length(unique(df$Patient)))

write.table(subset(df, ncbi_build==36), file ='June_23/TCGA_all_somatic.hg18.tsv', row.names= F, quote =F, sep = '\t')
write.table(subset(df, ncbi_build!=36), file ='June_23/TCGA_all_somatic.hg19.tsv', row.names= F, quote =F, sep = '\t')
#write.table(subset(df, !ncbi_build %in% c(36, 37)) , file ='June_23/TCGA_all_somatic.bad.txt', row.names= F, quote =F, sep = '\t')
