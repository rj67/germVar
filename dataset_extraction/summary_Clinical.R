# read in all tcga study clinical information
clin_files <- list.files(path = "../dataDump/TCGA/Clinical_bulk", recursive=T, pattern="clinical_patient", full.names=T)
#filename <- clin_files[15]

library(reshape2)
get_clinical <- function(filename){
  study <- gsub(".txt", "", strsplit(filename, split = "_patient_")[[1]][2], fixed=T)
  columnames <- scan(filename, what = "character", skip = 1, nlines = 1)
  clin<-read.delim(filename, header=F, skip = 3, stringsAsFactors=F, na.strings="[Not Available]")
  colnames(clin) <- columnames
  clin$study <- factor(toupper(study))
  to_select <- c("bcr_patient_barcode", "bcr_patient_uuid", "study", "gender", "race", "ethnicity",
                      "person_neoplasm_cancer_status", "vital_status", "pathologic_stage", "clinical_stage", "primary_therapy_outcome_success", "new_tumor_event_after_initial_treatment",
                        "days_to_birth", "age_at_initial_pathologic_diagnosis", "days_to_death", "days_to_last_followup" )
  clin <- clin[intersect(to_select, colnames(clin))]
  return(clin)
}

all_clin <- do.call(plyr::rbind.fill, lapply(clin_files, get_clinical))
library(DescTools)
#Desc(all_clin)
colnames(all_clin) <- gsub("age_at_initial_pathologic_diagnosis", "age", colnames(all_clin), fixed=T)
all_clin$Patient<-sapply(all_clin$bcr_patient_barcode, function(x) substr(x, 6, 12))

# figure the ethnicitiy
all_clin$EA <- with(all_clin, !is.na(race) & race=="WHITE" & !is.na(ethnicity) & ethnicity=="NOT HISPANIC OR LATINO") 

#barchart of patient age

library(RColorBrewer)
#all_clin <- droplevels(subset(all_clin, !is.na(age)))
p <- ggplot(aes(reorder(study, age, quantile, 0.5, na.rm=T), age), data = all_clin)
p <- p + geom_boxplot(outlier.shape = NA, guide="none") + geom_jitter( aes(color = study, alpha=0.6))
p <- p + xlab("TCGA study")
p <- p + scale_color_manual(values=colorRampPalette(brewer.pal(8, "Dark2"))(nlevels(all_clin$study)))
p <- p + theme( panel.background = element_rect(fill='white',colour='black'), axis.text= element_text(angle=90), legend.position="none")
p
ggsave(filename="All_patient_age_boxplot.pdf",p, width=8, height=6)
            #panel.grid.major = element_blank(),
            #panel.grid.minor = element_blank())
#t test of people with or without HR germline mutation
#t.test(BIWU.clinic$age[BIWU.clinic$GERMMUT==0], BIWU.clinic$age[BIWU.clinic$GERMMUT==1] )


#all_tcga <- plyr::join(all_tcga, all_clin[c("Patient", "race", "ethnicity")], by="Patient")
# figure out the ethnicity 
# cauca1 is a strict definition of caucasian, cauca2 is looser definition
#all_tcga$cauca1 <- with(all_tcga, !is.na(race) & race=="WHITE" & !is.na(ethnicity) & ethnicity=="NOT HISPANIC OR LATINO") 
#all_tcga$cauca2 <- with(all_tcga, (is.na(race) | race %in% c("[Not Evaluated]", "[Unknown]", "WHITE" )) & (is.na(ethnicity) | ethnicity!="HISPANIC OR LATINO")) 

save(all_clin, file="Results/all_tcga_clinical.RData")
rm(clin_files, get_clinical)
