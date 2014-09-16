
setwd("/Users/snafu/Documents/Project/germVar")

#all_tcga <- read.csv('input/all_TCGA_Aug_11.csv', header =T, stringsAsFactors = F)
all_tcga <- read.csv('input/all_TCGA.csv', header =T, stringsAsFactors = F)
print("raw cgquery")
print(dim(all_tcga))

# only illumina, hg19, center!=CGHUB
all_tcga <- subset(all_tcga, refassem_abbr =="19" & !(platform %in% c("ABI_SOLID", "LS454")) & center !="CGHUB")
# remove cell lines, weird sample type 12, 13
all_tcga <- subset(all_tcga, !sample_type %in% c("12", "13", "50"))
# remove LAML, DLBC, LUSC
#all_tcga <- subset(all_tcga, !disease %in% c("LAML", "DLBC", "LUSC"))
# remove filenames that contain HOLD_QC
all_tcga <- droplevels(all_tcga[-grep("HOLD", all_tcga$filename),])
all_tcga$Patient <- substr(all_tcga$sample_id, 6, 12)
all_tcga$Specimen <- substr(all_tcga$sample_id, 6, 16)
all_tcga$ToN <- sapply(all_tcga$sample_type, function(x) ifelse(x %in% c(10, 11), "N", "T"))
all_tcga$Year <- as.numeric(substr(all_tcga$upload_date, 1, 4))
all_tcga$Month <- as.numeric(substr(all_tcga$upload_date, 6, 7))
all_tcga$Day <- as.numeric(substr(all_tcga$upload_date, 9, 10))
print("after filter")
print(dim(all_tcga))

#bam_paths <- read.table("input/all_exist.list", strip.white=T, stringsAsFactors=F)
bam_paths <- read.table("input/all_bam_uid_SM.list", strip.white=T, stringsAsFactors=F)
colnames(bam_paths) <- c("analysis", "SM", "file_path")
bam_paths$Patient2 <- substr(bam_paths$SM, 6, 12)

# merge tcga dataframe with all existing bam files
all_tcga <- merge( all_tcga, bam_paths, by ="analysis")
# remove ~16 samples where the cgquery info dont match up with SM
print(c("number of samples where cgquery and SM dont match", sum(with(all_tcga, Patient!=Patient2))))
all_tcga <- subset(all_tcga, Patient==Patient2)

# sort by date and remove duplicates
all_tcga$exist <- all_tcga$analysis %in% bam_paths$analysis
all_tcga <- arrange(all_tcga, disease, Patient, ToN, legacy_sample_id, !exist, -Year, -Month, -Day)
print(c("number of duplicates", sum(duplicated(all_tcga$legacy_sample_id))))

all_tcga <- all_tcga[!duplicated(all_tcga$legacy_sample_id), ]

#filter out contaminated solid normals
solid_df <- read.table(file="output/SNPChip_solid_flagged.tsv" , stringsAsFactors=F, header=T)
print(c("number of contaminated samples", nrow(subset(all_tcga, Specimen %in% subset(solid_df, !flag)$Specimen))))
all_tcga <- subset(all_tcga, !Specimen %in% subset(solid_df, !flag)$Specimen)

# manually remove bam files that are actually aligned to HG18
all_tcga <- subset(all_tcga, !analysis %in% c("fc393f8e-5242-456c-bbcb-7edcadd5968a", "b9620e09-fec8-4533-9495-ab9d3720f9b9",
"349eeb2b-8f73-4d15-97cb-c2089921f6eb", "d2f523b7-6895-4e54-9b37-ba97d0985b55",
"721af97c-dff3-487e-97f1-259bbfc3c104", "5695946b-f5ea-458a-85b8-aea7bbd2fceb",
"71b9ff39-7d3c-459b-baf0-1beaa284629c", "6e82197f-1142-48ab-a933-b255f9d5c514",
"7d3a79d5-08f3-42a0-9995-318f894d2874", "6246ae12-0e5a-43b7-bd8e-48bb5e48da0b",
"b5a71869-acbe-4e81-a09b-05b609313521", "6a96f390-47b0-4412-8cf6-ea442ff21c8f",
"7ca06096-03fd-4ae7-998c-e3a3fdc8226a", "3b210e8d-3071-42bd-8b60-d8e23f99c8c8",
"45324dd9-ee4d-441f-a5c0-205267213175", "b596f377-c7f3-4bed-bb0c-e8c32a4deee1",
"31dcd6e7-c89f-4a2d-ae42-9de37de10a2b"))


#check each patient T/N status
ToN_stat <- dplyr::summarise(group_by(all_tcga, Patient), status = paste0(unique(ToN),collapse=""))
print(table(ToN_stat$status))
# remove samples where only Tumor sample is available
all_tcga <- merge(all_tcga, ToN_stat, by="Patient")
all_tcga <- subset(all_tcga, status !="T")

write.table(subset(all_tcga, sample_type %in% c(10, 11))[c("analysis", "file_path")], file="output/all_norm_bam_uid.list", sep="\t", row.names=F, col.names=F, quote=F)
write.csv(all_tcga,  file="output/all_tcga.csv", row.names=F, quote=F)
# tally
tally_tcga <- dplyr::summarise(group_by(all_tcga, disease, center, refassem_abbr), 
  num_patient=length(unique(participant)), num_sample=length(unique(sample_id)), num_analysis=length(unique(analysis)))
print(tally_tcga)

#dups <- subset(all_tcga, disease=="LUAD" & duplicated(sample_id))$sample_id
#View(subset(all_tcga, sample_id %in% dups))
#for (study in unique(all_tcga$disease)){
#  write(paste("/cbio/cslab/home/rj67/vcf/norm_redux/norm.", subset(all_tcga, disease==study&ToN=="N")$analysis, ".raw.gvcf", sep=""), file=paste0(c("output/gvcf_lists/", study, "_gvcf.list"), collapse=""))
#}
# add 4 new cancers
#all_new <- subset(all_tcga, disease %in% c("ESCA", "ACC", "PCPG", "UCS"))

# split off the 4 new cancers, look into cancer I've downloaded
#all_old <- subset(all_tcga, !disease %in% c("ESCA", "ACC", "PCPG", "UCS"))

# previous to June, existing bam files 
#exist_germ <- read.table('input/all_hg19_uid_SM.list', stringsAsFactors=F)
#exist_soma <- read.table('input/all_soma_uid_SM.list', stringsAsFactors=F)

#all_exist_uid<- c(exist_germ$V1, exist_soma$V1)
#all_exist_sample <- subset(all_tcga, analysis %in% all_exist_uid)$sample_id

#all_absent_uid <- subset(all_tcga, !sample_id %in% all_exist_sample)$analysis

#all_tcga_nr <- subset(all_tcga, analysis %in% c(all_exist_uid, all_absent_uid))
#all_tcga_nr <- subset(all_tcga, analysis %in% c(all_exist_uid, all_absent_uid))


rm(bam_paths, tally_tcga, solid_df, ToN_stat)
