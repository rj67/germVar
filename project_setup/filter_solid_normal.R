setwd("/Users/snafu/Documents/Project/germVar")
# read in all tcga germline sample SNP Chip file, flag suspicious
studies <- list.files(path = "../dataDump/TCGA/SNPARRAY", recursive=F, full.names=F)
studies <- studies[!grepl("tar", studies)]
print(studies)

######################### ##################################################
summarise_cnv_study <- function(study){
  #study <-"SKCM" 
  print(study)
  manifest <- read.delim( paste("../dataDump/TCGA/SNPARRAY/", study, "/FILE_SAMPLE_MAP.txt", sep=""), header =T, strip.white=T, stringsAsFactors = F)
  colnames(manifest)[2] <- "barcode"
  # only look at hg19 after common cnv removed
  manifest %<>% subset(., grepl("nocnv_hg19", filename)) %>%
                mutate(., Sample = substr(barcode, 1, 20), Tissue = substr(Sample, 14, 15), Specimen = substr(Sample, 6, 16)) %>%
                subset(., Tissue %in% c("11")) # only look at solid normal sample
  #print(manifest)
  if (nrow(manifest) > 0){
    manifest$disease <-  study
    solid_cnv <- do.call(rbind, lapply(manifest$filename, summarise_cnv_file, study=study))
    colnames(solid_cnv) <- paste("V", seq(1, ncol(solid_cnv)), sep="")
    solid_cnv <- cbind(manifest, as.data.frame(solid_cnv))
    return(solid_cnv)
  }
}


summarise_cnv_file <- function(sample_file, study){
  #sample_file <- manifest$filename[163]
  file_path <- paste("../dataDump/TCGA/SNPARRAY/", study, "/", sample_file, sep="");
  df <- read.delim(file_path, header =T, strip.white = T );
  #df <- subset(df, abs(Segment_Mean) >= 0.10 & Num_Probes >= 8)
  df <- subset(df, abs(Segment_Mean) >= 0.10 )
  df$Probe_cut <- cut(df$Num_Probes, breaks= c(0, 5^seq(1, 7)))
  #df$Probe_cut <- cut(df$Num_Probes, breaks= 10^seq(1, 5))
  df$Mean_cut <- cut(abs(df$Segment_Mean), breaks= c(0.1, 0.5, 10))
  cuts <- table(df$Probe_cut, df$Mean_cut)  
    #cuts
  cuts<-t(apply(cuts, 1, function(x) rev(cumsum(rev(x)))))
  cuts<-(apply(cuts, 2, function(x) rev(cumsum(rev(x)))))
  #total_probes <- sapply(c(0.1, 0.3), function(thresh) { 
  total_probes <- sapply(c(0.1, 0.5), function(thresh) { 
    num_probes <- with(df, (sum(Num_Probes[abs(Segment_Mean) >= thresh])));
    log_num <- ifelse(num_probes==0, 0, (num_probes));
    return(log_num)})
  return(c(sqrt(as.numeric(cuts)), total_probes))
}

copy_file <- function(sample_file, study){
   orig_path <- paste("../dataDump/TCGA/SNPARRAY/", study, "/", sample_file, sep="")
   dest_path <- "../dataDump/TCGA/solid_suspect"
   command <- paste("cp",orig_path, dest_path, sep=" ")
   system(command)
}

######################### ##################################################

# actually read the profiles
all_cnv<-do.call(rbind, lapply(studies, summarise_cnv_study))

# plot both solid/blood normal
ggplot(all_cnv, aes(V15)) + geom_bar(fill="darkgrey") + facet_grid(Tissue~., scale="free_y") + scale_x_log10() + geom_vline(xintercept=2000, color="red")

group_by(all_cnv, Tissue) %>% dplyr::summarise(., log_med = median(log10(V15)), log_sd = IQR(log10(V15))/1.349) 
# copy the solid suspect files
solid_cnv <- subset(all_cnv, Tissue==11)
solid_cnv$suspect <- solid_cnv$V15>2000
with(subset(solid_cnv, suspect ), mapply(copy_file, filename, disease))

# Output the solid normals
p <- ggplot(solid_cnv, aes(V15)) + geom_bar(fill="darkgrey") + theme_few() + scale_x_log10() + geom_vline(xintercept=2000, color="red")
p <- p + xlab('Number of probes with Segment Mean > 0.1') + ylab('Number of solid normal samples')
p 
ggsave(file="solid_normal_probe_hist.png", width=4, height=4)

write.csv(subset(solid_cnv, suspect)[c("filename", "barcode", "Sample", "disease", "V15")], quote=F, row.names=F, file="Output/SNPChip_solid_flagged.csv" )
save(solid_cnv, file="Results/")
rm(studies, all_cnv, copy_file, summarise_cnv_file, summarise_cnv_study)

#blood_cnv <- subset(all_cnv, Tissue==10)
#blood_cnv$suspect <- blood_cnv$V15>10000

#form feature matrix
# remove empty columns
#cnv_trans<- cnv_trans[colSums(cnv_trans) > 0]
# 

#   
# library(rrcov)
# # fit the solid tissue
# solid_pca<-PcaHubert(dplyr::select(solid_cnv, starts_with("V")), scale=T, mcd=T, trace=T, alpha=0.75, k=1)
# #screeplot(solid_pca)
# # two dimensional plot of PCA
# plot(solid_pca)
# 
# #extract info
# solid_cnv<-cbind(solid_cnv, as.data.frame( attributes(solid_pca)[c("od","sd","flag")]))
# table(solid_cnv$flag)
# 
# 
# print("cutoff.od")
# print(attributes(solid_pca)$cutoff.od)
# print("cutoff.sd")
# print(attributes(solid_pca)$cutoff.sd)
# print('filtered out samples')
# print(sum(!solid_df$flag))
# 
# write.csv(solid_cnv, quote=F, row.names=F, file="Output/SNPChip_solid_flagged.csv" )
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
#rm(studies, list_cnv, cnv_trans, manifest, solid_cnv, solid_pca, solid_df)
