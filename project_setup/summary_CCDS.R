setwd("/Users/snafu/Documents/Project/germVar")
################################################################
#     read the latest CCDS release, write to bed file
################################################################
release <- "../dataDump/CCDS/release_15"
read_CCDS <- function(release){
  CCDS_file <- paste(release, "CCDS.current.txt", sep="/")
  CCDS <- read.delim(CCDS_file)
  table(CCDS$ccds_status)
  # retain only Public entries
  CCDS <- subset(CCDS, ccds_status=="Public")
  # eleven genes not in CCDS
  nrow(subset(list_goi, Gene %in% setdiff(list_goi$Gene, CCDS$gene)))
  # subset CCDS to goi 
  CCDS <- droplevels(subset(CCDS, gene %in% list_goi$Gene))
  return(CCDS)
}


read_CCDS_seq <- function(release){
  # get the ensemble transcript ID
  Seq_file <- paste(release, "CCDS2Sequence.current.txt", sep="/")
  CCDS_seq <- read.delim(Seq_file)
  CCDS_seq <- subset(CCDS_seq, source =="EBI,WTSI" & grepl("ENST", nucleotide_ID))
  colnames(CCDS_seq) <- gsub("X.ccds", "ccds_id", colnames(CCDS_seq))
  return(CCDS_seq[c("nucleotide_ID", "ccds_id", "protein_ID")])
}

CCDS_17 <- read_CCDS("../dataDump/CCDS/release_17")
CCDS_17_seq <- read_CCDS_seq("../dataDump/CCDS/release_17")
CCDS_17_seq <- droplevels(subset(CCDS_17_seq, ccds_id %in% CCDS_17$ccds_id))
CCDS_17_seq <- merge(CCDS_17[c("ccds_id", "gene")], CCDS_17_seq)

CCDS_15 <- read_CCDS("../dataDump/CCDS/release_15")
CCDS_15_seq <- read_CCDS_seq("../dataDump/CCDS/release_15")
CCDS_15_seq <- droplevels(subset(CCDS_15_seq, ccds_id %in% CCDS_15$ccds_id))
CCDS_15_seq <- merge(CCDS_15[c("ccds_id", "gene")], CCDS_15_seq)

miss_genes <- setdiff(CCDS_17$gene, CCDS_15$gene)
rescued_seq <- subset(subset(CCDS_17_seq, gene %in% miss_genes), !nucleotide_ID %in% ensembl_gff$name)


ensembl_gff <- droplevels(subset(ensembl_gff, symbol %in% list_goi$Gene))
ensembl_gff$length <- mapply(function(x,y) y-x, ensembl_gff$start, ensembl_gff$end)
longest_transcript <- dplyr::summarise(group_by(ensembl_gff, symbol), max_length=max(length))

remove(rescue_genes,)


# get the exon number
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
