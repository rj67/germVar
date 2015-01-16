setwd("/Users/snafu/Documents/Project/germVar")

################################################################
#   Get GRCh37 ensembl transcripts
################################################################
load("Results/biomart_ensembl_gff.RData")
ensembl_gff <- droplevels(subset(ensembl_gff, symbol %in% list_goi$Gene))
ensembl_gff$length <- mapply(function(x,y) y-x, ensembl_gff$start, ensembl_gff$end)
#longest_seq <- ensembl_gff %>% arrange(., -length) %>% subset(., !duplicated(symbol))

library("biomaRt")
# up-to-date genome
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensembl <- useMart( biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
# Ensembl Gene ID
gene.ids <- list_goi[["Ensembl.ID"]]

retrieveEnsemblGFF <- function(ensembl, gene.ids){
  # what information to get for each transcript:
  sel.attributes=c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "chromosome_name", "strand", "start_position","end_position", 
                 "transcript_start", "transcript_end", "description", "transcript_biotype")
  # if 'hgnc_symbol' is not a defined attribute for gene symbols in your data species, try to find an equivalent, using commands such as this one:
  #grep("symbol",listAttributes(ensembl)[,1], value=TRUE)

  # retreive information:
  ensembl_gff <- getBM(attributes=sel.attributes, filters="ensembl_gene_id", value=gene.ids, mart=ensembl)

  ## replace attribute names by standardized names
  ensembl_gff <- ensembl_gff %>% plyr::rename(., c(
    "ensembl_gene_id" = "Ensembl.ID",
    "ensembl_transcript_id" = "Transcript",
    "chromosome_name" = "Chrom",
    "hgnc_symbol" = "Gene",
    "transcript_biotype" = "Biotype" 
  ))
  return(ensembl_gff)
}

ensembl_75 <- retrieveEnsemblGFF(ensembl, gene.ids)



################################################################
#     read the latest CCDS release, write to bed file
################################################################

read_CCDS <- function(release){
  # read the main CCDS file  
  CCDS_file <- paste(release, "CCDS.current.txt", sep="/")
  CCDS <- read.delim(CCDS_file)
  table(CCDS$ccds_status)
  # retain only Public entries
  CCDS <- subset(CCDS, ccds_status=="Public")
  # eleven genes not in CCDS
  nrow(subset(list_goi, Gene %in% setdiff(list_goi$Gene, CCDS$gene)))
  # subset CCDS to goi 
  CCDS <- droplevels(subset(CCDS, gene %in% list_goi$Gene))

  # get the ensemble transcript ID
  Seq_file <- paste(release, "CCDS2Sequence.current.txt", sep="/")
  CCDS_seq <- read.delim(Seq_file)
  CCDS_seq <- subset(CCDS_seq, source =="EBI,WTSI" & grepl("ENST", nucleotide_ID))
  colnames(CCDS_seq) <- gsub("X.ccds", "ccds_id", colnames(CCDS_seq))
  # merge with CCDS
  CCDS_seq <- merge(CCDS[c("ccds_id", "gene", "gene_id", "cds_strand")], CCDS_seq)
  
  # get the Exon number
  Exon_file <- paste(release, "CCDS_exons.current.txt", sep="/")
  CCDS_exon <- read.delim(Exon_file)
  # tally Max number of exon merge with CCDS_seq
  CCDS_seq <- plyr::join(CCDS_seq, CCDS_exon %>% group_by(., ccds_id) %>% dplyr::summarise(., ExonTotal = max(exon_ordinal_position)), by="ccds_id")
  
  # get the UniProt Number
  Uniprot_file <- paste(release, "CCDS2UniProtKB.current.txt", sep="/")
  CCDS_uniprot <- read.delim(Uniprot_file)
  colnames(CCDS_uniprot)[1] <- "ccds_id"
  CCDS_uniprot <- subset(CCDS_uniprot, ccds_id %in% CCDS_seq$ccds_id)
  CCDS_uniprot$uid <- interaction(CCDS_uniprot$ccds_id, CCDS_uniprot$UniProtKB, drop=T)
  CCDS_uniprot <- subset(CCDS_uniprot, !duplicated(uid))
  CCDS_seq <- plyr::join(CCDS_seq, CCDS_uniprot[c("ccds_id", "UniProtKB")], by="ccds_id")
  return(CCDS_seq[c("gene", "gene_id", "ccds_id", "nucleotide_ID",  "protein_ID", "UniProtKB","cds_strand", "ExonTotal")])
}

CCDS_17 <- read_CCDS("../dataDump/CCDS/release_17")

CCDS_15 <- read_CCDS("../dataDump/CCDS/release_15")

# genes that are not in CCDS_15 but in CCDS_17
miss_genes <- setdiff(CCDS_17$gene, CCDS_15$gene)
# transcripts that in CCDS_17 and also in ensembl 37
rescued_seq <- subset(subset(CCDS_17, gene %in% miss_genes), nucleotide_ID %in% ensembl_75$Transcript)
rescued_seq <- subset(CCDS_17, gene %in% miss_genes)
CCDS_15 <- rbind(CCDS_15, rescued_seq)

# manually check HG19 gene definition arrive at ENST00000428558 for RECQL4 
# SEMA3B deosn't have a protein coding transcript in GRCh37
# append to list_goi_seq
CCDS_15 <- plyr::rbind.fill(CCDS_15, data.frame(gene=c("RECQL4"), ccds_id =c("none"), nucleotide_ID=c("ENST00000428558"), UniProtKB=c("O94761"), cds_strand =c("-"), ExonTotal=c(22)))

subset(tmp, !nucleotide_ID %in% ensembl_75$Transcript)
problems <- subset(tmp, !nucleotide_ID %in% ensembl_75$Transcript)
# remove all but KMT2B, ENST00000222270
CCDS_15 <- subset(CCDS_15, !(nucleotide_ID %in% problems$nucleotide_ID) | nucleotide_ID=="ENST00000222270")

CCDS_15 <- merge(CCDS_15, ensembl_75[c("Transcript", "Ensembl.ID", "Biotype")], by.x="nucleotide_ID", by.y="Transcript", all.x=T)
length(unique(CCDS_15$gene))
length(unique(CCDS_15$ccds_id))
CCDS_15 <- subset(CCDS_15, is.na(Biotype) | Biotype=="protein_coding")

# dedup CCDS_15
CCDS_15 <- subset(CCDS_15, !duplicated(nucleotide_ID))

# still missing these genes, 
setdiff(miss_genes, CCDS_15$gene)
# take the longest transcript from ensembl
#subset(longest_seq, symbol %in% setdiff(miss_genes, CCDS_15_seq$gene))
# manually check HG19 gene definition arrive at ENST00000428558 for RECQL4 
# SEMA3B deosn't have a protein coding transcript in GRCh37
# append to list_goi_seq
#CCDS_15 <- plyr::rbind.fill(CCDS_15, data.frame(gene=c("RECQL4"), nucleotide_ID=c("ENST00000428558"), UniProtKB=c("O94761"), cds_strand =c("-"), ExonTotal=c(22)))


write(subset(CCDS_15, !duplicated(ccds_id)$nucleotide_ID), file="Output/list_goi_CCDS_transcript.txt")

save(CCDS_15, file="Results/list_goi_CCDS_transcript.RData")

remove(rescue_seq)


# get the exon number
if(F){
CCDS_exon <- read.delim(file="../dataDump/CCDS/CCDS_exons.current.release17.txt", stringsAsFactors=F)
CCDS <- CCDS_exon %>%
  group_by(ccds_id) %>%
  summarise(Num_exon = max(exon_ordinal_position)) %>%
  merge(., CCDS) 



CCDS <- merge(CCDS, CCDS_enst[c("ccds_id", "nucleotide_ID")])
  group_by(gene) %>%
  summarise(Num_ccds = length(unique(ccds_id)))

#table(CCDS$ccds_status)

#View(subset(CCDS, ccds_status!="Public"))
#setdiff(CCDS$gene, subset(CCDS, ccds_status=="Public")$gene)

# remove Y chromosome
#CCDS <- subset(CCDS, X.chromosome !="Y")
# get exon CDS
CCDS$cds_locations <- gsub("[", "", fixed=T, CCDS$cds_locations)
CCDS$cds_locations <- gsub("]", "", fixed=T, CCDS$cds_locations)
CCDS$cds_locations <- gsub(" ", "", fixed=T, CCDS$cds_locations)

tmp <- do.call(rbind, apply(CCDS[c("ccds_id", "cds_locations")], 1, function(x) {ranges = strsplit(x[2], split=",", fixed=T)[[1]] ;
                      data.frame(ccds_id = rep(x[1], length(ranges)), names=paste(x[1], seq(1, length(ranges)), sep="."),  ranges=as.character(ranges))}))
CCDS <- merge(CCDS, tmp, by="ccds_id")
#CCDS$ranges <- levels(CCDS$ranges)[CCDS$ranges]
CCDS$start <- sapply(CCDS$ranges, function(x) strsplit(x, split="-", fixed=T)[[1]][1])
CCDS$end <- sapply(CCDS$ranges, function(x) strsplit(x, split="-", fixed=T)[[1]][2])
CCDS$seqnames <- paste("chr", CCDS$X.chromosome, sep="")
# output bed for liftover
write.table(CCDS[c("seqnames", "start", "end", "names")], sep="\t", quote=F, row.names=F, col.names=F, file="Output/CCDS_goi.exons.hg38.bed")

################################################################
#     lift over CCDS bed file, read back in
################################################################
#read lifted over
liftover <- read.delim(file="../dataDump/CCDS/CCDS_goi.exons.hg38Tohg19.bed", stringsAsFactors=F, header=F)
colnames(liftover) <- c("seqnames", "start", "end", "names")
CCDS$start <- NULL
CCDS$end <- NULL
CCDS <- merge(CCDS, liftover[c("start", "end", "names")], by="names", all.x=T)
#compute the number of transcript per gene
CCDS <- merge(CCDS, dplyr::summarise(group_by(CCDS, gene), N_tran = length(unique(ccds_id))), by="gene")
CCDS <- merge(CCDS, dplyr::summarise(group_by(CCDS, ccds_id), N_exon = length(unique(ranges))), by="ccds_id")
subset(list_goi, Gene %in% setdiff(list_goi$Gene, CCDS$gene))
#mapply(queryCCDS, as.list, tmp[c("Gene", "POS", "REF", "ALT")])

library(IRanges)
CCDS_IRanges <- sapply(unique(CCDS$gene), function(x){
   tmp <- subset(CCDS, gene==x & !is.na(start));
   return(IRanges(start=tmp$start, end=tmp$end, names=tmp$ccds_id))
})

CCDS <- CCDS[!duplicated(CCDS$ccds_id), ]
save(CCDS, CCDS_IRanges, file="Results/CCDS_summary.RData")
}
