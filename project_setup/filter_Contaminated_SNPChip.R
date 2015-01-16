setwd("/Users/snafu/Documents/Project/germVar")
# read in all tcga germline sample SNP Chip file, flag suspicious
studies <- list.files(path = "../dataDump/TCGA/SNPARRAY", recursive=F, full.names=F)
print(studies)

######################### ##################################################
read_study <- function(study){
  #study <-"SKCM" 
  mani_path <- paste("../dataDump/TCGA/SNPARRAY/", study, "/FILE_SAMPLE_MAP.txt", sep="")
  manifest <- read.delim(mani_path, header =T, strip.white=T, stringsAsFactors = F)
  colnames(manifest)[2] <- "barcode"
  # only look at hg19 after common cnv removed
  manifest <- manifest[grep("nocnv_hg19", manifest$filename),]
  manifest$Sample <- substr(manifest$barcode, 1, 20)
  manifest$Tissue <- substr(manifest$Sample, 14, 15)
  # only look at normal sample
  manifest <- subset(manifest, Tissue %in% c("10", "11"))
  manifest$study <- study

  cnv_trans <- as.data.frame(do.call(rbind, lapply( manifest$filename, function(sample_file){
    #sample_file <- manifest$filename[163]
    file_path <- paste("../dataDump/TCGA/SNPARRAY/", study, "/normal_hg19/", sample_file, sep="");
    df <- read.delim(file_path, header =T, strip.white = T );
    df <- subset(df, abs(Segment_Mean) >= 0.10 & Num_Probes >= 8)
    df$Probe_cut <- cut(df$Num_Probes, breaks= 5^seq(1, 7)) 
    #df$Mean_cut <- cut(abs(df$Segment_Mean), breaks= c(0.10, 0.3, 0.6, 10))
    df$Mean_cut <- cut(abs(df$Segment_Mean), breaks= c(0.10, 0.3, 10))
    cuts <- table(df$Probe_cut, df$Mean_cut)  
    #cuts
    cuts<-t(apply(cuts, 1, function(x) rev(cumsum(rev(x)))))
    cuts<-(apply(cuts, 2, function(x) rev(cumsum(rev(x)))))
    total_probes <- sapply(c(0.1, 0.3), function(thresh) { 
      num_probes <- with(df, (sum(Num_Probes[abs(Segment_Mean) >= thresh])));
      log_num <- ifelse(num_probes==0, 0, log(num_probes));
      return(log_num)})
    return(c(sqrt(as.numeric(cuts)), total_probes))
  })))
  rownames(cnv_trans) <- manifest$filename
  return(list("manifest"=manifest, "cnv_trans"=cnv_trans))
}
######################### ##################################################

# actually read the profiles
list_cnv<-lapply(studies, read_study)
#form feature matrix
cnv_trans<-do.call(rbind, lapply(list_cnv, function(x) x[["cnv_trans"]]))
print(dim(cnv_trans))
# remove empty columns
cnv_trans<- cnv_trans[colSums(cnv_trans) > 0]
# remove empty rows
cnv_trans <- cnv_trans[cnv_trans$V1>0,]
print(dim(cnv_trans))
print(ht(cnv_trans))

# form file info matrix
manifest<-do.call(rbind, lapply(list_cnv, function(x) x[["manifest"]]))
manifest$Specimen <- substr(manifest$Sample, 6, 16)
#(ht(manifest))
print(table(manifest$Tissue))

# solid normal matrix
solid_cnv <- cnv_trans[rownames(cnv_trans) %in% subset(manifest, Tissue=="11")$filename, ]

#blood_cnv <- cnv_trans[rownames(cnv_trans) %in% subset(manifest, Tissue=="10")$filename, ]

library(rrcov)
# fit the solid tissue
solid_pca<-PcaHubert(solid_cnv, scale=T, mcd=T, trace=T, alpha=0.75, k=1)
#screeplot(solid_pca)
# two dimensional plot of PCA
plot(solid_pca)
#extract info
solid_df<-as.data.frame( attributes(solid_pca)[c("od","sd","flag")])
solid_df$filename<- rownames(solid_df)
solid_df<-merge(solid_df, manifest)

print("cutoff.od")
print(attributes(solid_pca)$cutoff.od)
print("cutoff.sd")
print(attributes(solid_pca)$cutoff.sd)
print('filtered out samples')
print(sum(!solid_df$flag))

write.table(solid_df, quote=F, row.names=F, sep = "\t", file="Output/SNPChip_solid_flagged.tsv" )
#solid_flagged<-arrange(subset(solid_df, od>attributes(solid_pca)$cutoff.od|sd>attributes(solid_pca)$cutoff.sd), study)
# include anything where V4>0
#misses<-setdiff(rownames(solid_cnv[solid_cnv$V4>0,]), solid_flagged$filename)
#solid_flagged <- rbind(solid_flagged, solid_df[solid_df$filename %in% misses,] )

# fit the blood normal
#blood_pca<-PcaHubert(blood_cnv, scale=T, mcd=T, trace=T, alpha=0.75, k=1)
#screeplot(blood_pca)
#plot(blood_pca)
#blood_df<-data.frame(od=attributes(blood_pca)[["od"]], sd=attributes(blood_pca)[["sd"]])
#blood_df$filename<- rownames(blood_df)
#blood_df<-merge(blood_df, manifest)
#blood_flagged<-arrange(subset(blood_df, od>5|sd>3), study)
#only flag the ones with long deletion
# fit the blood normal

#misses<-setdiff(rownames(solid_cnv[solid_cnv$V4>0,]), solid_flagged$filename)

#write.table(rbind(solid_flagged, blood_flagged)$filename, quote=F, col.names=F, row.names=F, file="output/cnv_flagged_names.list" )

#save(solid_flagged, solid_df, solid_cnv, blood_flagged, blood_df, blood_cnv, file="output/filter_SNPChip.RData")
rm(studies, list_cnv, cnv_trans, manifest, solid_cnv, solid_pca, solid_df)
