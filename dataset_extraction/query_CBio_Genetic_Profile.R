#get all the RNASeq data

# Cancer Genome Data Server
require('cgdsr')
# creat a connection to the server
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
# list all available cancer studies
getCancerStudies(mycgds)[, 1:2]
studies <- c("kirc", "coadread", "luad", "thca", "prad", "gbm", "ucec", "lgg", "cesc", "ov", "kirp", "brca", "skcm", "stad", 
             "paad", "blca", "lihc", "sarc", "hnsc", "acc", "kich", "ucs")

require(data.table)

CBio_query <- rbindlist( lapply(studies, function(study) {
  require(reshape2)
  print(study)
  study_id <- paste(study, "_tcga", sep="");
  case_list_id = paste(study_id, "_all", sep="");
  # log2CNA
  genetic_profile_id = paste(study_id, "_log2CNA", sep="");
  portion1 <- getProfileData(mycgds, list_goi$Gene[1:1000], genetic_profile_id, case_list_id);
  portion2 <- getProfileData(mycgds, list_goi$Gene[1001:nrow(list_goi)], genetic_profile_id, case_list_id);
  returns <- cbind(portion1, portion2)
  # make sure all rownames of return are the same length
  if (length(table(sapply(rownames(returns), nchar)))!=1) stop("rownames longer than 12")
  returns$Patient <- substr(gsub(".", "-", rownames(returns), fixed=T), 6, 12)
  returns <- melt(returns, id.vars= "Patient", variable.name="Gene", value.name="log2CNA")
  #returns$uid <- interaction(returns$Patient, returns$Gene, drop=T)
  
  # gistic
  genetic_profile_id = paste(study_id, "_gistic", sep="");
  portion1 <- getProfileData(mycgds, list_goi$Gene[1:1000], genetic_profile_id, case_list_id);
  portion2 <- getProfileData(mycgds, list_goi$Gene[1001:nrow(list_goi)], genetic_profile_id, case_list_id);
  returns2 <- cbind(portion1, portion2)
  if (length(table(sapply(rownames(returns2), nchar)))!=1) stop("rownames longer than 12")
  returns2$Patient <- substr(gsub(".", "-", rownames(returns2), fixed=T), 6, 12)
  returns2 <- melt(returns2, id.vars= "Patient", variable.name="Gene", value.name="gistic")
  #returns2$uid <- interaction(returns2$Patient, returns2$Gene, drop=T)
  
  # mrnaz
  genetic_profile_id = paste(study_id, "_rna_seq_v2_mrna_median_Zscores", sep="");
  portion1 <- getProfileData(mycgds, list_goi$Gene[1:1000], genetic_profile_id, case_list_id);
  portion2 <- getProfileData(mycgds, list_goi$Gene[1001:nrow(list_goi)], genetic_profile_id, case_list_id);
  returns3 <- cbind(portion1, portion2);
  if (length(table(sapply(rownames(returns3), nchar)))!=1) stop("rownames longer than 12");
  returns3$Patient <- substr(gsub(".", "-", rownames(returns3), fixed=T), 6, 12);
  returns3 <- melt(returns3, id.vars= "Patient", variable.name="Gene", value.name="mrnaz");
  #returns3$uid <- interaction(returns3$Patient, returns3$Gene, drop=T)
  
  # merge all  
  returns <- cbind(returns, gistic=returns2$gistic, mrnaz=returns3$mrnaz)
  returns$study <- study
  return(returns)
}))
CBio_query <- subset(CBio_query, !(is.na(log2CNA) & is.na(gistic) & is.na(mrnaz) ))
CBio_query$event_uid <- apply(CBio_query[,c("Patient", "Gene"),with=F], 1, function(x) paste0(x, collapse="-") )
CBio_query$gistic2 <- sign(CBio_query$gistic)

#CBio_query <- subset(CBio_query, Patient %in% all_tcga$Patient)

save(CBio_query, file="Results/CBio_genetic_profile.RData")
rm(studies)
