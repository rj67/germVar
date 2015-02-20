ASCAT_stat <- read.table("Results/ASCAT_out/ploidy_purity.txt", fill=T)
colnames(ASCAT_stat) <- c("Patient", "ploidy", "purity")
ASCAT_stat <- subset(ASCAT_stat, !is.na(purity))

seq_lengths<- data.frame(CHROM = c(seq(1:22), "x"), length=seqlengths(seqinfo(nonsyn_GT))[1:23])

parse_ascat <- function(Patient){
  ascat_file <- paste("Results/ASCAT_out/", Patient, ".segments_withSize.txt", sep="") 
  ascat_df <- read.csv(ascat_file, strip.white = T)
  # some file have the "Different chromosome boundaries!" error message, remove them
  ascat_df <- subset(ascat_df, !grepl("Different", Segment..))
  ascat_grange <- GRanges( seqnames = Rle(ascat_df$chromosome), 
                           ranges = IRanges(start = ascat_df$startBP, end = ascat_df$endBP, names = ascat_df$Segment..))
  return(list(df=ascat_df, grange=ascat_grange))
}

ASCAT <-lapply(ASCAT_stat$Patient, parse_ascat) 
A

tmp<-subset(nsSNP_GT, uid %in% subset(nsSNP_rv, pred_patho=="tier1")$uid)
tmp<-plyr::join(tmp, nsSNP_vr[c("uid", "aa_uid")])
tmp$uid <- as.character(levels(tmp$uid)[tmp$uid])
tmp2<-mapply(query_ascat, tmp$Patient, tmp$uid)
tmp2<-as.data.frame(t(tmp2))

df<-ASCAT[[25]]$df
plot_ASCAT <- function(df){
  library(grid)
  # reorder chrmomosome level
  df$chromosome <- factor(df$chromosome, levels=c(seq(1, 22), "X"))
  # for small segments, extend both end to achieve at least size_thresh 
  size_thresh <- 500000
  df$startBP[df$size<size_thresh] <- df$startBP[df$size<size_thresh] - (size_thresh -df$size[df$size<size_thresh])/2
  df$endBP[df$size<size_thresh] <- df$endBP[df$size<size_thresh] + (size_thresh -df$size[df$size<size_thresh])/2
  # cap nA at threshold
  nA_thresh <- 5
  # offset nA nB for plotting
  nAB_offset <- 0.09
  df$nA[df$nA>nA_thresh] <- nA_thresh
  df$nA <- df$nA - nAB_offset
  df$nB <- df$nB + nAB_offset
  p <- ggplot(df) + geom_segment(aes(x=startBP, xend=endBP, y=nA, yend=nA), size=5, color="#d7191c") 
  p <- p + geom_segment(aes(x=startBP, xend=endBP, y=nB, yend=nB), size=5, color="#2c7bb6") 
  p <- p + facet_grid(.~chromosome, scales="free_x", space="free") 
  p <- p + theme_few() + ylab("") + xlab("") + theme(axis.text.x=element_blank()) 
  p <- p + scale_y_continuous(breaks=seq(0, 5)) + scale_x_continuous(breaks=NULL)
  p
}
ggsave(filename="tmp.png", width=10, height=5)
