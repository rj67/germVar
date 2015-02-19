# RNASeq V2
require(data.table)
mrna_files <- list.files(path = "./Results/RNASeqV2", pattern="RData", recursive=F, full.names=T, )
RNASeq <- rbindlist( lapply(mrna_files, function(x){
  load(file=x);
  mrna_object <- gsub(".RData", "", strsplit(x, split="/")[[1]][4], fixed=T) ;
  call_set <- get(mrna_object);
  call_set <- data.table(call_set, key="Gene,Patient");
  rm(mrna_object);
  return(call_set)
}))
RNASeq_norm<-subset(RNASeq, sample_type==11)
RNASeq <- subset(RNASeq, sample_type!=11)
RNASeq <- RNASeq[!duplicated(RNASeq$event_uid), ]

# RNASeq old
mrna_files <- list.files(path = "./Results/RNASeq", pattern="RData", recursive=F, full.names=T, )
RNASeq_old <- rbindlist( lapply(mrna_files, function(x){
  load(file=x);
  mrna_object <- gsub(".RData", "", strsplit(x, split="/")[[1]][4], fixed=T) ;
  call_set <- get(mrna_object);
  call_set <- data.table(call_set, key="Gene,Patient");
  rm(mrna_object);
  return(call_set)
}))
# take the largest isoform
RNASeq_old <- subset(RNASeq_old, sample_type=="01")

RNASeq_old <- merge(RNASeq_old[!duplicated(RNASeq_old$event_uid),],  dplyr::summarise(group_by(RNASeq_old, event_uid, isof), normalized_count = max(RPKM)), by="event_uid")
# merge the RNASeq_old
RNASeq<-rbind(RNASeq, RNASeq_old[, colnames(RNASeq), with = F])


# across board Gene stat
Gene_stat <- dplyr::summarise(group_by(RNASeq, Gene), med = median(sqrt(normalized_count)), sd2 = IQR(sqrt(normalized_count))/1.349)
no_express <- subset(Gene_stat, med==0 & sd2==0)$Gene
Gene_stat <- subset(Gene_stat, ! Gene %in% no_express)

# remove genes that don't express
RNASeq <- subset(RNASeq, !Gene %in% no_express)

Gene_study_stat <- dplyr::summarise(group_by(RNASeq, study, platform, Gene), med = median(sqrt(normalized_count)), sd = IQR(sqrt(normalized_count))/1.349)
Gene_noexpr_stat <-dplyr::summarise(group_by(subset(Gene_study_stat, sd==0), Gene), noexpr_tissue=length(unique(study)))
subset(Gene_noexpr_stat, noexpr_tissue>13)

# substitue the cross board sd for when there is no sd
Gene_study_stat <- merge(Gene_study_stat, Gene_stat[, c("Gene", "sd2"), with=F], by= "Gene")
Gene_study_stat$sd <- apply(Gene_study_stat[,c("sd", "sd2"), with=F], 1, function(x) ifelse(x[1]!=0, x[1], x[2]))

# merge the med and sd back into RNASeq
RNASeq <- merge(RNASeq, Gene_study_stat[, c("study", "platform", "Gene", "med", "sd"), with=F], by=c("study", "platform", "Gene"))
RNASeq$mrnaz <- with(RNASeq, (sqrt(normalized_count)-med)/sd)

RNASeq <- subset(RNASeq, Patient %in% all_tcga$Patient)
save(RNASeq, no_express, file="Results/RNASeq_summary.RData")
