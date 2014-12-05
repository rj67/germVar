CCDS <- read.delim(file="../dataDump/CCDS/CCDS.current.release17.txt", stringsAsFactors=F)
#setdiff(list_goi$Gene, CCDS$gene)

# subset CCDS to goi
CCDS <- subset(CCDS, gene %in% list_goi$Gene)
#table(CCDS$ccds_status)

# retain only Public entries
CCDS <- subset(CCDS, ccds_status=="Public")
#View(subset(CCDS, ccds_status!="Public"))
#setdiff(CCDS$gene, subset(CCDS, ccds_status=="Public")$gene)

# remove Y chromosome
CCDS <- subset(CCDS, X.chromosome !="Y")
# get exon CDS
CCDS$cds_locations <- gsub("[", "", fixed=T, CCDS$cds_locations)
CCDS$cds_locations <- gsub("]", "", fixed=T, CCDS$cds_locations)
CCDS$cds_locations <- gsub(" ", "", fixed=T, CCDS$cds_locations)

tmp <- do.call(rbind, apply(CCDS[c("ccds_id", "cds_locations")], 1, function(x) {ranges = strsplit(x[2], split=",", fixed=T)[[1]] ;
                      data.frame(ccds_id = rep(x[1], length(ranges)), names=paste(x[1], seq(1, length(ranges)), sep="."),  ranges=as.character(ranges))}))
CCDS <- merge(CCDS, tmp, by="ccds_id")
CCDS$ranges <- levels(CCDS$ranges)[CCDS$ranges]
CCDS$start <- sapply(CCDS$ranges, function(x) strsplit(x, split="-", fixed=T)[[1]][1])
CCDS$end <- sapply(CCDS$ranges, function(x) strsplit(x, split="-", fixed=T)[[1]][2])
CCDS$seqnames <- paste("chr", CCDS$X.chromosome, sep="")
# output bed for liftover
write.table(CCDS[c("seqnames", "start", "end", "names")], sep="\t", quote=F, row.names=F, col.names=F, file="output/CCDS_goi.exons.hg38.bed")

#read lifted over
liftover <- read.delim(file="../dataDump/CCDS/CCDS_goi.exons.hg38Tohg19.bed", stringsAsFactors=F, header=F)
colnames(liftover) <- c("seqnames", "start", "end", "names")
CCDS$start <- NULL
CCDS$end <- NULL
CCDS <- merge(CCDS, liftover[c("start", "end", "names")], by="names", all.x=T)
#compute the number of transcript per gene
CCDS <- merge(CCDS, dplyr::summarise(group_by(CCDS, gene), N_tran = length(unique(ccds_id))), by="gene")

#mapply(queryCCDS, as.list, tmp[c("Gene", "POS", "REF", "ALT")])

require(IRanges)
CCDS_IRanges <- sapply(unique(CCDS$gene), function(x){
   tmp <- subset(CCDS, gene==x & !is.na(start));
   return(IRanges(start=tmp$start, end=tmp$end, names=tmp$ccds_id))
})
