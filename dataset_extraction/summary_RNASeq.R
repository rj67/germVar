library(ff)
library(ffbase)
# read in the house keeping gene list
House_Keep <- read.delim(file="input/NHKS_homo_sapien.txt", skip=3)
setdiff(House_Keep$GID, table_HGNC$Ensembl.ID)
House_Keep <- subset(table_HGNC, Ensembl.ID %in% House_Keep$GID)[c("Approved.Symbol", "Entrez.Gene", "Approved.Name")]
colnames(House_Keep)[1] <- "Gene"
save(House_Keep, file="Results/house_keep.RData")

# RNASeq V2
library(data.table)
mrna_files <- list.files(path = "./Results/RNASeqV2", pattern="RData", recursive=F, full.names=T, )
RNASeq <- lapply(mrna_files, function(x){
  load(file=x);
  mrna_object <- gsub(".RData", "", strsplit(x, split="/")[[1]][4], fixed=T) ;
  de_object <- gsub("mrna", "de", mrna_object) ;
  mrna <- get(mrna_object);
  mrna <- data.table(mrna);
  de <- get(de_object);
  study <- gsub("_mrna", "", mrna_object);
  mrna$study <- study;
  de$study <- study;
  return(list(RNASeq=mrna, de=de))
})
#RNASeq_norm<-subset(RNASeq, sample_type==11)
RNASeq <- subset(RNASeq, sample_type!=11)
RNASeq <- RNASeq[!duplicated(RNASeq$event_uid), ]
tmp1 <- do.call(rbind, lapply(RNASeq, function(x) x[[1]]))
tmp2 <- do.call(rbind, lapply(RNASeq, function(x) x[[2]]))

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

tmp <- do.call(rbind, lapply(RNASeq, function(x) x[[2]]))

RNASeq_old <- merge(RNASeq_old[!duplicated(RNASeq_old$event_uid),],  dplyr::summarise(group_by(RNASeq_old, event_uid, isof), normalized_count = max(RPKM)), by="event_uid")
# merge the RNASeq_old
RNASeq<-rbind(RNASeq, RNASeq_old[, colnames(RNASeq), with = F])
# remove event_uid
RNASeq$event_uid <- NULL
# only take profiled patients
RNASeq <- subset(RNASeq, Patient %in% all_tcga$Patient)

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

RNASeq$study <- factor(RNASeq$study)
RNASeq$platform <- factor(RNASeq$platform)
RNASeq$Patient <- factor(RNASeq$Patient)
RNASeq$sample_type <- factor(RNASeq$sample_type)
RNASeq$barcode <- NULL

RNASeq <- as.ffdf(RNASeq)
ffsave(RNASeq, file="DataSet_Results/RNASeq_summary.ff")

#save(RNASeq, no_express, file="Results/RNASeq_summary.RData")



load("Results/RNASeqV2/ACC_mrna.RData")
tmp<-subset(ACC_mrna, HK & sample_type=="01")
library(reshape2)
tmp2 <- dcast(tmp[c("Gene", "Patient", "normalized_count")], Gene~Patient)
hist(apply(as.matrix(tmp2[-1]), 1, function(x) sd(x, na.rm=T)/mean(x, na.rm=T)))

library(princomp)
fit <- prcomp(t(as.matrix(tmp2[-1])), cor=TRUE)
summary(fit)
biplot(fit)

loadings(fit)
plot(fit,type="lines")
