################################################################
#     get three sets all together
################################################################
load("Results/TCGA_recurrent_variants.RData")
load("Results/COSMIC_recurrent_variants.RData")
nonsyn_var <- do.call(rbind, list(recur_cosm[c("CHROM", "POS", "REF", "ALT", "Gene", "var_uid", "aa_uid", "site")], 
                    recur_tcga[c("CHROM", "POS", "REF", "ALT", "Gene", "var_uid", "aa_uid", "site")],
                    clinvar_txt[c("CHROM", "POS", "REF", "ALT", "Gene", "var_uid", "aa_uid", "site")] ))
nonsyn_var <- nonsyn_var[!duplicated(nonsyn_var$var_uid),]
nonsyn_var <- merge(nonsyn_var, recur_tcga[c("var_uid", "tcga_vcount", "tcga_scount")], by="var_uid", all.x=T)
nonsyn_var <- merge(nonsyn_var, recur_cosm[c("var_uid", "cosm_vcount", "cosm_scount")], by="var_uid", all.x=T)
nonsyn_var <- merge(nonsyn_var, clinvar_txt[c("var_uid", "RCVaccession", "PhenotypeIDs", "Origin","ClinicalSignificance", "OtherIDs")], by="var_uid", all.x=T)
#writeBED(nonsyn_var, "output/nonsyn_var.bed")

#manual literature curation

#generation of the csv 
#write.csv(subset(nonsyn_calls, var_uid %in% VAR$var_uid)@VAR[colnames(tmp)], file="Results/out.csv", row.names=F)
tmp<- read.csv(file="Results/nonsyn_variant_pathogenicity_curation.csv", header=T, stringsAsFactors=F)
nonsyn_var <- merge(nonsyn_var, tmp[c("var_uid","Pathogenicity","Cancer", "Note")], by="var_uid", all.x=T)

#tmp<-clin_var[c("Gene", "AA_Change", "RCVaccession", "ClinicalSignificance", "OtherIDs")]
#tmp<-merge(table_HGNC[c("Approved.Symbol", "UniProt.ID")], tmp, by.y="Gene", by.x="Approved.Symbol")
#colnames(tmp)[1] <- "Gene"
#write.table(tmp, file="output/clinvar_cancer_gene_nonsyn_variants.tsv", sep="\t", quote=F, row.names=F, col.names=T)
################################################################
#       nonsyn Variant INFO
################################################################
library(VariantAnnotation)
# variant info
# nonsyn_VAR <- expand(readVcf(file="./Results/norm_nonsyn.VAR.vcf.gz", "hg19"))
# 
# # preprocess
# nonsyn_VAR <- nonsyn_VAR %>% labelUid  %>%   
#               subset(., Coding == "CODING" & VARTYPE=="SNP") %>%
#               subset(., BIOTYPE == "protein_coding" ) %>%
#               fixSnpEffGene %>%
#               subset(., Gene %in% list_goi$Gene) %>%
#               subset(., Impact=="MODERATE") %>%
#               subset(., Gene==SYMBOL)   %>% 
#               labelVarUid %>% 
#               fixVEPCCDS %>% 
#               labelVEPUid %>% 
#               convSIFT %>% 
#               labelAAUid %>% 
#               labelSite %>% 
#               reduceCOL 

nonsyn_VAR <- 

# preprocess
nsSNP_VAR <- expand(readVcf(file="./Results/norm_nsSNP.VAR.vcf.gz", "hg19")) %>% labelUid %>%   
  subset(., !is.na(Coding) & Coding == "CODING" & VARTYPE=="SNP") %>%
  fixSnpEffGene %>%
  subset(., Gene %in% list_goi$Gene) %>%
  subset(., Impact=="MODERATE") %>%
  labelVarUid %>% 
  labelSnpEffUid %>%
  labelAAUid %>% 
  labelSite %>% 
  reduceCOL %>%
  mergeVCF(., list_gff_long[c("Transcript", "strand")], by="Transcript") %>%
  annoESPX2kG(.)



#### get the mutation info
## intersect with clinvar/COSMIC/TCGA mutations
nsSNP_VAR <- subset(nsSNP_VAR, site %in% c(cosm_nonsyn$site, clinvar_txt$site ))
## COSMIC 
#use merge to get cosm_vcount and cosm_scount
nsSNP_VAR <- mergeVCF(nsSNP_VAR, cosm_nonsyn[c("var_uid", "cosm_vcount")], by="var_uid")
nsSNP_VAR <- mergeVCF(nsSNP_VAR, subset(cosm_nonsyn, !duplicated(site))[c("site", "cosm_scount")], by="site")
## Clinvar
# flag whether site is mutated in Clinvar
info(header(nsSNP_VAR))["clin_site",] <- list("0" ,"Flag", "Whether site in clinvar list")
info(nsSNP_VAR)$clin_site <- info(nsSNP_VAR)$site %in% clinvar_txt$site
# get the ClinVar annotation by merging on aa_uid 
nsSNP_VAR <- mergeVCF(nsSNP_VAR, subset(clinvar_txt, !duplicated(aa_uid))[c("aa_uid",  "RCVaccession", "PhenotypeIDs", "ClinicalSignificance", "OtherIDs")], by="aa_uid") 


# label whether exact variant match
#info(header(nonsyn_VAR))["var_match",] <- list("0" ,"Flag", "Whether variant in curated list")
#info(nonsyn_VAR)$var_match <- info(nonsyn_VAR)$var_uid %in% nonsyn_var$var_uid
# label whether aachange match
#info(header(nonsyn_VAR))["aa_match",] <- list("0" ,"Flag", "Whether amino acid change in curated list")
#info(nonsyn_VAR)$aa_match <- info(nonsyn_VAR)$aa_uid %in% nonsyn_var$aa_uid
# only keep rare variants if not aa_match
#nonsyn_VAR <- subset(nonsyn_VAR, aa_match | AF<=0.02)
# fix the origin column
#info(nonsyn_VAR)$Origin[grep("germline", info(nonsyn_VAR)$Origin)] <- "germline" #anything containing germline is germline
#info(nonsyn_VAR)$Origin[!info(nonsyn_VAR)$Origin %in% c("germline", "somatic")] <- NA # anythin not germline or somatic, designate NA first
#info(nonsyn_VAR)$Origin[is.na(info(nonsyn_VAR)$Origin) & (!is.na(info(nonsyn_VAR)$tcga_scount) | !is.na(info(nonsyn_VAR)$cosm_scount))] <- "somatic" # not labeled in ClinVar, labeled in cosm/tcga are somatic
#info(nonsyn_VAR)$Origin[is.na(info(nonsyn_VAR)$Origin)] <- "other" # the remaining are "other"
info(nonsyn_VAR)$Origin <- c("Clinvar","HSMutation")[as.numeric(with(info(nonsyn_VAR), !is.na(tcga_scount)|!is.na(cosm_scount))) +1 ]

# classify the ClinVar
info(header(nsSNP_VAR))["Clinvar",] <- list("1" ,"String", "Classify clinvar annotation")
info(nsSNP_VAR)$Clinvar <- sapply(info(nsSNP_VAR)$ClinicalSignificance, function(x) ifelse(is.na(x),"" ,"Patho/Risk"))
info(nsSNP_VAR)$Clinvar[grepl("uncertain",info(nsSNP_VAR)$ClinicalSignificance,  ignore.case=T  )|
                           grepl("benign",info(nsSNP_VAR)$ClinicalSignificance,  ignore.case=T  ) ]<-"Conflicting"
info(nsSNP_VAR)$Clinvar[is.na(info(nsSNP_VAR)$Clinvar)] <- "Unknown"
info(nsSNP_VAR)$Clinvar <- factor(info(nsSNP_VAR)$Clinvar, levels=c("Patho/Risk", "Conflicting", "Unknown"))

######################
# bioinformatic prediction
######################

# output for FATHMM variant annotation
# view(info(nonsyn_VAR)[c("ENSP", "vep_aa")])
# fathmm <- read.delim("Results/nonsyn_VAR.annotated.fathmm.disease.tsv", header=T)
# fathmm <- fathmm %>% plyr::rename(., replace = c("Protein.ID"="ENSP",  "Substitution"="AAChange",	"Prediction"="fathmm_pred", "Score"= "fathmm_score"))
# fathmm$fathmm_pred[fathmm$fathmm_pred==""] <- "no weights"
# #fathmm$fathmm_pred <- factor(fathmm$fathmm_pred, levels=c("CANCER", "PASSENGER/OTHER", "no weights"))
# fathmm$fathmm_pred <- factor(fathmm$fathmm_pred, levels=c("DAMAGING", "TOLERATED", "no weights"))
# info(header(nonsyn_VAR))["fathmm_pred",] <- list("1" ,"String", "fathmm prediction")
# info(header(nonsyn_VAR))["fathmm_score",] <- list("1" ,"Float", "fathmm score")
# info(nonsyn_VAR)$fathmm_pred <- fathmm$fathmm_pred
# info(nonsyn_VAR)$fathmm_score <- fathmm$fathmm_score

# ### unique variant set, different VEP alter variants are ranked by fathmm prediction
# uniq_vep <- as.data.frame(info(nonsyn_VAR)) %>% arrange(., fathmm_pred) %>% subset(., !duplicated(var_uid)) %>% dplyr::select(., matches("vep_uid"))
# nonsyn_VARu <- subset(nonsyn_VAR, vep_uid %in% uniq_vep$vep_uid)
# # annotate ESP X2kG AF
# nonsyn_VARu <- annoESPX2kG(nonsyn_VARu)
# 
# View(subset(info(nonsyn_VAR), var_uid %in%  subset(tmp, out!=1)$var_uid)[c("var_uid", "vep_aa", "SIFT", "SIFT_score", "PolyPhen","Cscore", "fathmm_pred", "fathmm_score")])
# 
# View(subset(info(nonsyn_VARu), pred_patho=="Deleterious")[c("var_uid", "vep_aa", "SIFT", "SIFT_score", "PolyPhen","Cscore", "fathmm_pred", "fathmm_score")])

# antoher one, separates into 3 tiers
#info(header(nonsyn_VARu))["pred_patho3",] <- list("1" ,"Integer", "Bioinformatic prediction")
#info(nonsyn_VARu)$pred_patho3 <- as.numeric(info(nonsyn_VARu)$fathmm_pred=="DAMAGING") + as.numeric(info(nonsyn_VARu)$Cscore>17)

### CADD
# output for CADD variant annotation
writeVcf(nsSNP_VAR, filename="Output/nsSNP_VAR_CADD.txt")
# read back CADD score and merge into info
CADD <- read.delim("Results/nsSNP_VAR.annotated.CADD.tsv", header=T, skip=1)
CADD$uid <- apply(CADD[c("X.Chrom", "Pos", "Ref", "Alt")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
# 2 variants are not annotated nsSNP in CADD,  this include , CLDN23-Pro232Ser, MEF2B-Arg17Trp  
CADD <- CADD %>% subset(., GeneName %in% list_goi$Gene)
nsSNP_VAR <- mergeVCF(nsSNP_VAR, CADD[c("uid", "RawScore", "PHRED")], by="uid")
# set NA CADD score to 0.001,CLDN23-Pro232Ser, MEF2B-Arg17Tr
info(nsSNP_VAR)$PHRED[is.na(info(nsSNP_VAR)$PHRED)] <- 0.001


# output for Mutation Assessor variant annotation
writeMA(nsSNP_VAR, filename="Output/nsSNP_VAR_MA.txt")
MutAss <- read.csv("Results/nsSNP_VAR.annotated.MutationAssessor.csv", header=T)
#MutAss$Func..Impact[MutAss$Func..Impact==""] <- "neutral"
## Mutation Assessor has different gene names
# Gene  AAChange.p	MA_Gene	AA.variant
# AGO1	Arg88His	EIF2C1	R88H
# BRINP3	Arg598Trp	FAM5C	R598W
# ACKR3	Arg161Cys	CXCR7	R161C
# UVSSA	Ala48Thr	KIAA1530	A48T
# KDM7A	Arg644Cys	JHDM1D	R644C
# KMT2C	Ala4490Thr	MLL3	A4433T
# NUTM2A	Ala160Thr	FAM22A	A88T
# KDM4E	Arg100Cys	KDM4DL	R100C
# KMT2A	Leu2368Val	MLL	L2365V
# KMT2D	Arg4420Trp	MLL2	R4420W
# ERCC5	Ala874Thr	BIVM-ERCC5	A874T
# NUTM1	Ser165Phe	C15orf55	S137F
# SNX29	Val357Leu	RUNDC2A	V204L
# MEF2B	Arg17Trp	MEF2BNB-MEF2B	R17W
# CRLF2	Arg215Trp	CRLF2 /// CRLF2	
# P2RY8	Thr331Met	P2RY8 /// P2RY8	T331M
# AMER1	Ala29Thr	FAM123B	A29T
info(nsSNP_VAR)[MutAss$Gene!=info(nsSNP_VAR)$Gene,]
MutAss[MutAss$Gene!=info(nonsyn_VARu)$Gene,]
info(header(nsSNP_VAR))["ma_pred",] <- list("1" ,"String", "MutAss prediction")
info(nsSNP_VAR)$ma_pred <- MutAss$Func..Impact
# combine MutationAssessor and CADD prediction
info(header(nsSNP_VAR))["pred_patho",] <- list("1" ,"String", "Bioinformatic prediction")
info(nsSNP_VAR)$pred_patho <- sapply((info(nsSNP_VAR)$ma_pred %in% c("medium","high")) & (info(nsSNP_VAR)$PHRED>=20), function(x) ifelse(x, "Deleterious", "Benign"))
info(nsSNP_VAR)$pred_patho <- factor(info(nsSNP_VAR)$pred_patho, levels=c("Deleterious", "Benign"))

################################################################
#       nonsyn Variant GT
################################################################
# read in genotype by chromosome
readChrGT <- function(chr){
  print(chr)
  print("read in GT", quote=F)
  dummy_GT <- (readVcf(file=paste0("./Results/norm_nsSNP.", chr ,".GT.vcf.gz", collapse=""), "hg19"))
  print("filter GT", quote=F)
  dummy_GT <- dummy_GT  %>% labelUid %>% subset(., uid %in% info(nsSNP_VAR)$uid) 
  print("sumamrise GT", quote=F)
  return(summaryVCFGT(dummy_GT))
}
nsSNP_tally <- lapply( c(as.character(seq(1, 22)), "X"), readChrGT )
# merge the variant tally info
nsSNP_VAR <- mergeVCF(nsSNP_VAR, do.call(rbind, lapply(nsSNP_tally, function(x) x$tally)))
# sum up the Genotype
nsSNP_GT <-  do.call(rbind, lapply(nsSNP_tally, function(x) x$GT))
# remove mis call entry
nsSNP_VAR <- subset(nsSNP_VAR, !is.na(EAF))

# for filtering
nsSNP_VAR <- labelDomiAllele(nsSNP_VAR)
# quality filter 
nsSNP_VAR <- hardFilter(nsSNP_VAR)

# 6 commom variants have AB_med between 0.2~0.25, only CDK11B found in ExAC, other are suspect
# CDK11B-1571841-A-C
# NOTCH2-120611960-C-T
# MST1-49726070-G-A
# PARP4-25021323-A-G
# PKD1-2164211-G-A
# CDC27-45214519-A-G

# park common variants away
nsSNP_VAR_comm <- subset(nsSNP_VAR, EAF>=0.01 & (X2kG_AF>=0.01 | ESP_AF_>=0.01))
# now remove common variants from futher consideration
nsSNP_VAR <- subset(nsSNP_VAR, !var_uid %in% info(nsSNP_VAR_comm)$var_uid)

# remove two variants that are not in database but common in our data set, not sure if FP
nonsyn_VARu <- subset(nonsyn_VARu, !var_uid %in% c("RSBN1L-77379331-T-G", "KDM4E-94759020-G-A"))
nonsyn_VARu <- subset(nonsyn_VARu, pass)

# write FP input
all_uid <- nonsyn_calls@VAR$var_uid
chunks <- split(all_uid, ceiling(seq_along(all_uid)/100))
for (i in seq(1,3)){
  print(i)
  writeFPfilter(subset(nonsyn_calls, var_uid %in% chunks[[i]]), paste0(c("output/nonsyn_fpinp.", i, ".tsv"), collapse=""))
}

################################################################
#       nonsyn both tumor and normal
################################################################
VAR <- read.delim("Results/nonsyn_fpout.VAR.tsv", strip.white=T, stringsAsFactors=F, na.strings = "")
VAR <- subset(VAR, MQRankSum >-12)
VAR <- merge(VAR, nonsyn_calls@VAR[c("var_uid", "Gene", "AAChange")])
GT <- read.delim("Results/nonsyn_fpout.GT.tsv", strip.white=T, stringsAsFactors=F, na.strings = "")
GT <- merge(GT, VAR[c("ALT_idx", "Gene", "AAChange","var_uid")])
GT$is_var <- with(GT, GT %in% c("0/1", "1/1"))
GT$SAMPLE_ADs <- sapply(GT$AD, function(x) as.numeric(strsplit(x, split=",")[[1]]))
GT$SAMPLE_DP <- sapply(GT$SAMPLE_ADs, sum)
GT$SAMPLE_DP[is.na(GT$SAMPLE_DP)] <- 0
GT$SAMPLE_ALT_AD <- mapply(function(idx, ADs) {if(is.na(idx)){return(ADs[2])}else{return(ADs[idx+1])}}, GT$ALT_idx, GT$SAMPLE_ADs)
GT$SAMPLE_ALT_AD[is.na(GT$SAMPLE_ALT_AD)] <- 0
GT$SAMPLE_AB <-  mapply(function(ALT_AD, DP) {if(DP==0){return(0)}else{return(ALT_AD/DP)}}, GT$SAMPLE_ALT_AD, GT$SAMPLE_DP)    
GT <- merge(GT, all_tcga, by.x="SAMPLE", by.y="SM")
GT$mut_uid <- apply(GT[c("Patient", "var_uid")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
GT$event_uid <- apply(GT[c("Patient", "Gene")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))


#ggplot(aes(NAC), data=VAR) + geom_bar()+ facet_grid(~ patho)+ theme_set(theme_bw())
  
mut_stat <-  dplyr::summarise(group_by(GT, Gene, var_uid, Patient,event_uid, mut_uid), N_het = any(GT=="0/1"& (SAMPLE_DP-SAMPLE_ALT_AD) >= 3 & ToN=="N"), 
                              T_hom = any(SAMPLE_AB>=0.75 & SAMPLE_DP >= 20), LOH = N_het & T_hom, PIT = any(is_var & ToN=="T") | all(ToN=="N"))

VAR <- merge(VAR, dplyr::summarise(group_by(GT, var_uid), NAC = length(unique( mut_uid [is_var & ToN =="N"])), TAC = length(unique( mut_uid [is_var & ToN =="T"])) ,
                                   NAB_med = median(SAMPLE_AB[ToN=="N"]), TAB_med = median(SAMPLE_AB[ToN=="T"])))
VAR <- merge(VAR, dplyr::summarise(group_by(mut_stat, var_uid), LOH = sum(LOH)))

#VAR$patho <- VAR$var_uid %in% patho_uid
nonsyn_calls<- merge(nonsyn_calls, VAR[c("var_uid", "NAC", "TAC", "NAB_med", "TAB_med", "LOH")], by="var_uid")
nonsyn_miss_calls <- subset(nonsyn_calls, TAC==0)
nonsyn_calls <- subset(nonsyn_calls, TAC>0)
# remove two variants that are suspect
nonsyn_calls <- subset(nonsyn_calls, !var_uid %in% c("HNF1A-121437407-A-G", "JAK2-5073770-G-T"))

# missing in tumor
View(arrange(subset(GT, mut_uid %in% subset(mut_stat, !PIT)$mut_uid & ToN=="N" & var_uid %in% patho_uid), 
             Patient)[c("Patient", "Gene", "AAChange", "GT", "AD", "SAMPLE_ALT_AD", "SAMPLE_AB", "SAMPLE", "var_uid", "analysis", "disease", "center")])
# remove missing in tumor variants and sample
mut_stat <- subset(mut_stat, PIT)
#GT <- subset(GT, mut_uid %in% mut_stat$mut_uid)
VAR <- subset(VAR, var_uid %in% mut_stat$var_uid)


nonsyn_calls <- merge(nonsyn_calls, VAR[c("var_uid", "NAC", "TAC", "patho")])
#patho_calls<-subset(nonsyn_calls, !patho)

#look at the frequent variants
View(subset(nonsyn_calls, AC>15)@VAR[c("var_uid","AAChange","AC","NAC","TAC","AF", "ESP_AF", "X1kG_AF", "X1kG_LDAF")])

####################################
# compare AF###
####################################

tmp <- cbind(nonsyn_calls@VAR, t(apply(nonsyn_calls@VAR[c("NAC", "AN", "X2kG_AC", "X2kG_AN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "X2kG"))})))
tmp <- cbind(tmp, t(apply(nonsyn_calls@VAR[c("NAC", "AN", "ESP_EA_AC", "ESP_EA_AN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "ESP"))})))
tmp$pass <- with(tmp, ESP.sign*X2kG.sign>0)
View(arrange(subset(tmp, NAC>=15)[c("var_uid", "Gene", "AAChange", "Pathogenicity","NAC", "AN","ESP_EA_AC", "ESP_EA_AN", "ESP_EA_AF", "X2kG_AC", "X2kG_AN", "X2kG_AF", 
                  "ESP.pval", "ESP.sign", "ESP.conf.lo", "ESP.conf.hi", "X2kG.pval", "X2kG.sign","X2kG.conf.lo", "X2kG.conf.hi")],  Gene))

View(arrange(subset(tmp, NAC<15 & (ESP.sign==-2 | X2kG.sign==-2))[c("var_uid", "Gene", "AAChange", "patho","NAC", "AN","ESP_EA_AC", "ESP_EA_AN", "ESP_AF", "X2kG_AC", "X2kG_AN", "X2kG_AF", 
                                   "ESP.pval", "ESP.sign", "ESP.conf.lo", "ESP.conf.hi", "X2kG.pval", "X2kG.sign","X2kG.conf.lo", "X2kG.conf.hi", "pass")], !patho, Gene))

# create a subset of pathogenic variants associated with cancer
patho_calls <- subset(nonsyn_calls, patho & !Cancer=="No")
patho_calls <- subset(patho_calls, NAC < 15 | var_uid %in% c("RET-43613908-A-T", "TRAF3-103342015-C-T"))

####################################
##look at disease distribution
####################################

apply( subset(nonsyn_calls@VAR, NAC>15&patho), 1, function(x){
  print(x["var_uid"]);
  to_test <- subset(nonsyn_calls@GT, var_uid==x["var_uid"]);
  bg <- subset(CBio_query, Gene==x["Gene"])
  print(chisqTestCallSet(to_test, "gistic", bg))
})



patient_hit<-dplyr::summarise(group_by(patho_calls@GT, Patient), hit = length(unique(var_uid)))
subset(patient_hit, hit>1)
mut_stat<- merge(mut_stat, CBio_query[c("event_uid", "log2CNA", "gistic", "mrnaz")], by="event_uid", all.x=T) 


#VAR <- merge(VAR, nonsyn_calls@VAR[c("var_uid", setdiff(colnames(nonsyn_calls@VAR), colnames(VAR)))], by="var_uid")
#tmp<-new("VariantCallSet", GT = GT, VAR = VAR)
#tmp <- subset(tmp, Pathogenicity=="Yes" & Cancer!="No")

patho_uid <- subset(nonsyn_calls@VAR, (Pathogenicity=="Yes" & Cancer!="No"))$var_uid

# add somatic expression copy number  info
load("Results/CBio_genetic_profile.RData")
nonsyn_calls@GT<- merge(nonsyn_calls@GT, CBio_query[c("event_uid", "log2CNA", "gistic", "mrnaz")], all.x=T) 

patho_calls <- subset(nonsyn_calls, Pathogenicity=="Yes" & Cancer!="No")

#notes <- read.csv(file="Results/recur_calls_pathogenicity.csv")
#nonsyn_calls <- merge(nonsyn_calls, notes[c("var_uid", "Pathogenicity", "Note")], all.x=T)
#tmp <- nonsyn_calls@VAR[c("var_uid", "Gene", "AAChange", "AC", "AF", "ESP_AC", "ESP_AF", "X1kG_AF", "tcga_count", "cosm_count", "ClinicalSignificance", "OtherIDs", "Pathogenicity", "Note")]
#sapply(tmp, function(x), sum(is.na()))
#manual literature
table(tmp$Pathogenicity)
#tcga
sum(!is.na(tmp$tcga_count))
sum(!is.na(tmp$cosm_count))
sum(!is.na(tmp$cosm_count)|!is.na(tmp$tcga_count))
sum(!is.na(tmp$ClinicalSignificance))
view(tmp ) 

writeVCF(nonsyn_calls, "output/nonsyn_check_anno.vcf")


tmp2<-do(group_by(tmp@GT, var_uid), pval = chisqTestCallSet(., "disease"))
