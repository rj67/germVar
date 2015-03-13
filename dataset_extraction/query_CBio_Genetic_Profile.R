setwd("/Users/snafu/Documents/Project/germVar")

#get all the RNASeq data

# Cancer Genome Data Server
require('cgdsr')
# creat a connection to the server
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
# list all available cancer studies
getCancerStudies(mycgds)[, 1:2]
studies <- c("kirc", "coadread", "luad", "lusc", "thca", "prad", "gbm", "ucec", "lgg", "cesc", "ov", "kirp", "brca", "skcm", "stad", 
             "paad", "blca", "lihc", "sarc", "hnsc", "acc", "kich", "ucs", "pcpg", "esca")

library(data.table)

# Only Protein data, 18 studies and 99 Genes
CBio_rppa <- rbindlist( lapply(studies, query_rppa))
CBio_rppa <- subset(CBio_rppa, !is.nan(rppa))
CBio_rppa<-as.data.frame(CBio_rppa)
CBio_rppa$Patient <- factor(CBio_rppa$Patient)
CBio_rppa$disease <- factor(CBio_rppa$disease)
CBio_rppa <- as.ffdf(CBio_rppa)
ffsave(CBio_rppa, file="DataSet_Results/CBio_rppa_profile.ff")

# The other stuff, 
CBio_query <- rbindlist( lapply(studies, query_study))
CBio_query <- subset(CBio_query, !(is.na(log2CNA) & is.na(gistic) & is.na(mrnaz) ))
CBio_query<-as.data.frame(CBio_query)
CBio_query$Patient <- factor(CBio_query$Patient)
CBio_query$disease <- factor(CBio_query$disease)

library(ff)
library(ffbase)
CBio_query <- as.ffdf(CBio_query)
ffsave(CBio_query, file="DataSet_Results/CBio_genetic_profile.ff")
rm(studies)


query_rppa <- function(study) {
  library(reshape2)
  print(study)
  study_id <- paste(study, "_tcga", sep="");
  case_list_id = paste(study_id, "_all", sep="");
  
  # check genetic profile data
  getGeneticProfiles(mycgds, "cesc_tcga")
  
  # RPPA
  genetic_profile_id <- paste(study_id, "_RPPA_protein_level", sep="")
  if( genetic_profile_id %in% getGeneticProfiles(mycgds, study_id)$genetic_profile_id ){
    returns <- query_data(genetic_profile_id, case_list_id, "rppa")
    returns$disease <- toupper(study)
    return(returns)
  }
}

query_study <- function(study) {
  library(reshape2)
  print(study)
  study_id <- paste(study, "_tcga", sep="");
  case_list_id = paste(study_id, "_all", sep="");
  
  # check genetic profile data
  #getGeneticProfiles(mycgds, "cesc_tcga")
  
  # log2CNA
  returns <- query_data(paste(study_id, "_log2CNA", sep=""), case_list_id, "log2CNA")
  # gistic copy number call
  returns2 <- query_data(paste(study_id, "_gistic", sep=""), case_list_id, "gistic")
  # RNASeq V2 zscore
  #if not v2, then retry
  if(paste(study_id, "_rna_seq_v2_mrna_median_Zscores", sep="") %in% getGeneticProfiles(mycgds, study_id)$genetic_profile_id){
    returns3 <- query_data(paste(study_id, "_rna_seq_v2_mrna_median_Zscores", sep=""), case_list_id, "mrnaz")
  }else if(paste(study_id, "_rna_seq_mrna_median_Zscores", sep="") %in% getGeneticProfiles(mycgds, study_id)$genetic_profile_id){
    returns3 <- query_data(paste(study_id, "_rna_seq_mrna_median_Zscores", sep=""), case_list_id, "mrnaz")
  }else{
    return3 <- NA
  }
  # merge all  
  if(class(return3)=="data.frame"){
    returns <- cbind(returns, gistic=returns2$gistic, mrnaz=returns3$mrnaz)
  }else{
    returns <- cbind(returns, gistic=returns2$gistic, mrnaz=rep(NA, nrow(returns2)))
  }
  returns$disease <- toupper(study)
  return(returns)
}

query_data <- function(genetic_profile_id, case_list_id, value_name){
  portion1 <- getProfileData(mycgds, list_goi$Gene[1:500], genetic_profile_id, case_list_id);
  portion2 <- getProfileData(mycgds, list_goi$Gene[501:1000], genetic_profile_id, case_list_id);
  portion3 <- getProfileData(mycgds, list_goi$Gene[1001:1500], genetic_profile_id, case_list_id);
  portion4 <- getProfileData(mycgds, list_goi$Gene[1501:nrow(list_goi)], genetic_profile_id, case_list_id);
  returns <- cbind(portion1, portion2, portion3, portion4)
  # make sure all rownames of return are the same length
  if (length(table(sapply(rownames(returns), nchar)))!=1) stop("rownames longer than 12")
  returns$Patient <- substr(gsub(".", "-", rownames(returns), fixed=T), 6, 12)
  returns <- melt(returns, id.vars= "Patient", variable.name="Gene", value.name=value_name)
  return(returns)
}


