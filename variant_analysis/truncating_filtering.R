
trunc_VAR <- expand(readVcf(file="./Results/norm_trunc.VAR.vcf.gz", "hg19"))

# preprocess
trunc_VAR <- trunc_VAR %>% labelUid  %>%   
  subset(., Coding == "CODING" ) %>%
  subset(., BIOTYPE == "protein_coding" ) %>%
  fixSnpEffGene %>%
  subset(., Gene %in% list_goi$Gene) %>%
  subset(., Impact=="HIGH" & EFF!="STOP_LOST") %>% 
  subset(., Gene==SYMBOL)  %>% 
  labelVarUid %>% 
  fixVEPCCDS %>% 
  labelVEPUid(., "trunc") %>% 
  reduceCOL(., "trunc") 

# remove 6 genes that dont have CCDS
no_ccds<-as.data.frame(info(trunc_VAR)[c("var_uid", "ccds_id")]) %>% group_by(var_uid) %>% dplyr::summarise(., num_ccds = length(unique(ccds_id[!is.na(ccds_id)]))) %>% subset(., num_ccds==0)
table(subset(info(trunc_VAR), var_uid %in% no_ccds$var_uid)$Gene)
trunc_VAR <- trunc_VAR %>% subset(., !is.na(ccds_id)) 

# set LoF NA equal unknown
info(trunc_VAR)$LoF[is.na(info(trunc_VAR)$LoF)] <- "unknown"
info(trunc_VAR)$LoF <- factor(info(trunc_VAR)$LoF, levels=c("HC", "LC", "unknown"))

# make sure same vep_uid has same LoF flag
as.data.frame(info(trunc_VAR)[c("vep_uid", "LoF")]) %>% group_by(vep_uid) %>% dplyr::summarise(., num_lof = length(unique(LoF))) %>% subset(., num_lof>1)

# rank vep_uid according to LoF, only take one vep_uid per var_uid
uniq_vep <- as.data.frame(info(trunc_VAR)) %>% arrange(., LoF) %>% subset(., !duplicated(var_uid)) %>% dplyr::select(., matches("vep_uid"))
# get the uniq var_uid set
trunc_VARu <- subset(trunc_VAR, vep_uid %in% uniq_vep$vep_uid) %>% subset(., !duplicated(var_uid))

#trouble with annoESP and uid
trunc_VARu <- subset(trunc_VARu, uid!="19-57334157-G-A")
trunc_VARu <- annoESPX2kG(trunc_VARu)

# tally the total high confidence variants in 1000G
comm_genes <- as.data.frame(info(trunc_VARu)) %>% subset(., LoF=="HC") %>% group_by(., Gene) %>% dplyr::summarise(., tot_X2kG = sum(X2kG_AF)) %>% subset(., tot_X2kG>=0.02) 
comm_genes

# remove these genes
trunc_VARu <- subset(trunc_VARu, !Gene %in% comm_genes$Gene)


apply(trunc_calls@VAR[c("Gene", "POS", "REF", "ALT")], 1, function(x) do.call(queryCCDS, as.list(x)))


################################################################
#       Trunc Variant GT
################################################################
# Genotype
chrs <- c(as.character(seq(1, 22)), "X")
# read in genotype by chromosome
chr_tally <- list()
for (chr in chrs){ 
  print(chr)
  print("read in GT", quote=F)
  trunc_GT <- (readVcf(file=paste0("./Results/norm_trunc.GT.", chr ,".vcf.gz", collapse=""), "hg19"))
  #update GT 
  print("filter GT", quote=F)
  trunc_GT <- trunc_GT  %>% labelUid %>% subset(., uid %in% info(trunc_VARu)$uid) 
  print("sumamrise GT", quote=F)
  chr_tally[[chr]] <- summaryVCFGT(trunc_GT)
}
trunc_tally <- list()
trunc_tally$tally <- do.call(rbind, lapply(chr_tally, function(x) x$tally))
# six variants are missing
subset(info(trunc_VARu), !uid%in% trunc_tally$tally$uid)
trunc_VARu <- subset(trunc_VARu, uid %in% trunc_tally$tally$uid)
# get the nonsyn_GT variant summary, EAC, EAN, sample AD and AB, etc
# update nonsyn_VARu with the tallied info
trunc_VARu <- mergeVCF(trunc_VARu, trunc_tally$tally, by="uid")

# for filtering
trunc_VARu <- labelDomiAllele(trunc_VARu)
# quality filter 
trunc_VARu <- hardFilter(trunc_VARu)

# 
trunc_VARu <- subset(trunc_VARu, pass & LoF!="LC")


#######################################################################################################################################################################
#### Previous code
#######################################################################################################################################################################



load("Results/norm_trunc.merged.annovar.RData")
comment(call_set)
trunc_calls <- prepRawCallSet(call_set)
trunc_calls <- subset(trunc_calls, Impact=="HIGH")
trunc_calls@VAR$CDS <- apply(trunc_calls@VAR[c("Gene", "POS", "REF", "ALT")], 1, function(x) do.call(queryCCDS, as.list(x)))
trunc_calls <- subset(trunc_calls, CDS !=0)
trunc_calls@VAR <- annoESPX2kG(trunc_calls@VAR)

#####################################################################
#### remove suspicious samples 
#####################################################################
# tally the number of truncating variants per sample
sample_tally <- data.frame(table(subset(trunc_calls@GT, SAMPLE_ALT_AD >= 3 & SAMPLE_AB>=0.15)$SAMPLE))
#sample_tally <- data.frame(table(subset(trunc_calls@GT)$SAMPLE))
sample_tally <- arrange(merge(sample_tally, all_tcga, by.x="Var1", by.y="SM"), -Freq)
plot(subset(sample_tally, Freq>15)$Freq)
ggplot(subset(sample_tally, filesize>2&analyte!="X"), aes(Freq))   + geom_bar(stat="bin")+scale_y_sqrt()+facet_grid(.~analyte)
ggsave("trunc_var_per_sample_tally_histogram.png", width=6, height=5)
abline(60, 0)
View(subset(sample_tally, Freq>75)) 
trunc_calls@GT <- merge(trunc_calls@GT, sample_tally[c("Var1", "Freq")], all.x=T, by.x="SAMPLE", by.y="Var1")
# remove the top samples with more than 75 truncating variants
clean_GT <- droplevels(subset(trunc_calls@GT, !SAMPLE %in% subset(sample_tally, Freq>25)$Var1))
# number of filtered out variants
length(setdiff(trunc_calls@GT$var_uid, clean_GT$var_uid))
# corresponding VAR slot
filtered_out <- subset(trunc_calls@VAR, !var_uid %in% trunc_calls@GT$var_uid)
trunc_calls@VAR <-subset(trunc_calls@VAR, var_uid %in% trunc_calls@GT$var_uid)

#####################################################################
#### remove some variants, common, near tail 
#####################################################################
# remove reference genome annotation error
subset(trunc_calls, AF>=0.99)@VAR
trunc_calls <- subset(trunc_calls, AF<0.99)
# remove STOP_LOST
trunc_calls <- subset(trunc_calls, EFF!="STOP_LOST")
# remove variants near end of protein
trunc_calls <-subset(trunc_calls, (is.na(AA_pos) | abs(AA_pos-AALength)/AALength > 0.05 | (Pfam==1 & AftDom==0)))
# remove genes that have common truncating variants in 1000G, note 1000G doesn't have splice site 
View(subset(trunc_calls,(X1kG_LDAF>=0.005|X2kG_AF>=0.005)  &  QD>4  & (abs(AA_pos-AALength)/AALength > 0.20 | (Pfam==1 & AftDom==0)))@VAR)

# show all the homozygous TMPRSS11A samples
table(subset(subset(trunc_calls, var_uid=="TMPRSS11A-68829109-C-T")@GT, SAMPLE_AB>0.8)$disease)
# show all the homozygous PRODH samples
table(subset(subset(trunc_calls, var_uid=="PRODH-18912677-C-T")@GT, SAMPLE_AB>0.8)$disease)
# remove genes that have frequent germline truncation in 1000G
common_loh <- unique(subset(trunc_calls,(X1kG_LDAF>=0.05|X2kG_AF>=0.05)  &  QD>4  & (abs(AA_pos-AALength)/AALength > 0.20 | (Pfam==1 & AftDom==0)))@VAR$Gene)
lowfq_loh <- unique(subset(trunc_calls,(X2kG_AF>=0.005)  &  QD>4  & (abs(AA_pos-AALength)/AALength > 0.20 | (Pfam==1 & AftDom==0)))@VAR$Gene)
lowfq_loh <- setdiff(lowfq_loh, common_loh)
lowfq_calls <- subset(trunc_calls, Gene%in% lowfq_loh & ( X2kG_AF>=0.005)  &  QD>4  & (abs(AA_pos-AALength)/AALength > 0.20 | (Pfam==1 & AftDom==0)))

# test copy number change
apply( lowfq_calls@VAR, 1, function(x){
  print(x["var_uid"]);
  to_test <- subset(lowfq_calls@GT, var_uid==x["var_uid"]);
  if(x["Gene"] %in% CBio_query$Gene){
    bg <- subset(CBio_query, Gene==x["Gene"])
    print(wilcoxTestCallSet(to_test, "gistic2", bg))
  }else{
    print("not in CBioquery")
  }
})
# test disease
lapply( lowfq_calls@VAR$var_uid , function(x){ print(x); to_test <- subset(lowfq_calls@GT, var_uid==x);
    print(chisqTestCallSet(to_test, "disease"))
})
# test allele frequency
apply(lowfq_calls@VAR[c("AC", "AN", "ESP_AC", "ESP_AN")], 1, function(x) {names(x) <-NULL;do.call(varBurden, c(as.list(x),"CNTR", "greater"))})
#lowfq_loh <- lowfq_loh[lowfq_loh!="POLM"]
trunc_freq <- subset(trunc_calls, Gene %in% c(common_loh, lowfq_loh) )
trunc_calls <- subset(trunc_calls, !Gene %in% c(common_loh, lowfq_loh) )

################################################################
#       trunc Variant GT
################################################################
trunc_GT <- (readVcf(file="./Results/norm_trunc.GT.vcf.gz", "hg19"))
# update both GT and VAR
trunc_GT <- trunc_GT  %>% labelUid %>% subset(., uid %in% info(trunc_VARu)$uid) %>% calcAC
trunc_VARu <- mergeVCF(trunc_VARu, as.data.frame(info(trunc_GT)[c("uid", "EAC", "EAN", "EAF", "CAC1", "CAN1", "CAC2", "CAN2")]), by="uid")
#  
trunc_VARu <- labelDomiAllele(trunc_VARu)
# 6 variants have NA for QD, not sure why, filter out
trunc_VARu <- subset(trunc_VARu, !is.na(QD))
tmp<- preFilter(trunc_VARu)


#####################################################################
#### read in hg19 annotated 
#####################################################################
writeVCF(trunc_freq, "output/freq_check_anno.vcf")
# run snpEFF with hg19 on above VCF
load("Results/trunc_check_anno.snpEff.RData")
trunc_hg19 <- subset(call_set, Gene %in% list_goi$Gene)
View(subset(trunc_calls@VAR, !var_uid %in% subset(trunc_hg19, Impact=="HIGH")$var_uid))
# remove trunc_calls variants that are not labeled HIGH in trunc_hg19
#trunc_calls <- subset(trunc_calls, var_uid %in% subset(trunc_hg19, Impact=="HIGH")$var_uid)


#####################################################################
#### hard filter
#####################################################################
trunc_lowq <- hardFilter(trunc_calls, F)
trunc_filt <- hardFilter(trunc_calls, T)
#temporarily remove splice site
trunc_filt <- subset(trunc_filt, !is.na(AA_pos))
# remove some frequent ones
# GEN1, RAD52, CBLC for nearing the end
# CDC27 for weird, not in 1000G or ESP
trunc_filt <- subset(trunc_filt, !var_uid %in% c("CBLC-45296846-A-AC","RAD52-1023218-G-T", "GEN1-17962993-CAAGTT-C", "CDC27-45214558-G-A", "CDC27-45234406-CA-C"))

#####################################################################
#### joint analysis Tumor/Normal
#####################################################################
# write FP input
all_uid <- trunc_filt@VAR$var_uid
chunks <- split(all_uid, ceiling(seq_along(all_uid)/200))
for (i in seq(1,21)){
  print(i)
  writeFPfilter(subset(trunc_filt, var_uid %in% chunks[[i]]), paste0(c("output/trunc_fpinp.", i, ".tsv"), collapse=""))
}
View()

require(plyr)
VAR_files <- list.files(path = "./Results/trunc_fpout", pattern="VAR", recursive=F, full.names=T, )
GT_files <- list.files(path = "./Results/trunc_fpout", pattern="GT", recursive=F, full.names=T, )
VAR <- do.call(rbind.fill, lapply(VAR_files, function(x) read.delim(x, strip.white=T, stringsAsFactors=F, na.strings = "")))
VAR <- merge(VAR, trunc_filt@VAR[c("var_uid", "Gene")])
GT <- do.call(rbind, lapply(GT_files, function(x) read.delim(x, strip.white=T, stringsAsFactors=F, na.strings = "")))
GT <- merge(GT, VAR[c("ALT_idx", "Gene", "var_uid")])

GT$is_var <- with(GT, GT %in% c("./1", "0/1", "1/1"))
GT$SAMPLE_ADs <- sapply(GT$AD, function(x) as.numeric(strsplit(x, split=",")[[1]]))
GT$SAMPLE_DP <- sapply(GT$SAMPLE_ADs, sum)
GT$SAMPLE_DP[is.na(GT$SAMPLE_DP)] <- 0
GT$SAMPLE_ALT_AD <- mapply(function(idx, ADs) {if(is.na(idx)){return(ADs[2])}else{return(ADs[idx+1])}}, GT$ALT_idx, GT$SAMPLE_ADs)
GT$SAMPLE_ALT_AD[is.na(GT$SAMPLE_ALT_AD)] <- 0
GT$SAMPLE_AB <-  mapply(function(ALT_AD, DP) {if(DP==0){return(0)}else{return(ALT_AD/DP)}}, GT$SAMPLE_ALT_AD, GT$SAMPLE_DP)    
GT <- merge(GT, all_tcga, by.x="SAMPLE", by.y="SM")
GT$mut_uid <- apply(GT[c("Patient", "var_uid")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
GT$event_uid <- apply(GT[c("Patient", "Gene")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
VAR <- merge(VAR, dplyr::summarise(group_by(GT, var_uid), NAC = length(unique( mut_uid [is_var & ToN =="N"])), TAC = length(unique( mut_uid [is_var & ToN =="T"])) ,
                                                          NACq = length(unique( mut_uid [is_var & ToN =="N" & SAMPLE_ALT_AD>=3 & SAMPLE_AB >= 0.15 ])), 
                                                          TACq = length(unique( mut_uid [is_var & ToN =="T" & SAMPLE_ALT_AD>=3 & SAMPLE_AB >= 0.15 ])) ,
                                                          NALT_AD_med = median(SAMPLE_ALT_AD[ToN=="N"]), TALT_AD_med = median(SAMPLE_ALT_AD[ToN=="T"]),
                                                          NAB_med = median(SAMPLE_AB[ToN=="N"]), TAB_med = median(SAMPLE_AB[ToN=="T"]), Nonly = all(ToN=="N")))
VAR <- annoESPX2kG(VAR)
VAR <- subset(VAR, TAC!=0 | Nonly)
VAR <- subset(VAR, NAB_med > 0)
VAR <- subset(VAR, QD > 2)
VAR <- subset(VAR, QD > 4 | AC<20)
VAR <-subset(VAR, FS < 60 | VARTYPE!= "SNP")
VAR <-subset(VAR, FS < 200) 
VAR <- subset(VAR, ReadPosRankSum> -8 | VARTYPE!="SNP")

GT <- subset(GT, var_uid %in% VAR$var_uid)
mut_stat <-  merge(GT[!duplicated(GT$mut_uid),][c("Gene", "var_uid", "Patient", "event_uid", "mut_uid")], dplyr::summarise(group_by(GT, mut_uid), N_het = any(GT=="0/1"& (SAMPLE_DP-SAMPLE_ALT_AD) >= 3 & ToN=="N"), 
                              T_hom = any(SAMPLE_AB>=0.75 & SAMPLE_DP >= 20), LOH = N_het & T_hom, PIT = any(is_var & ToN=="T") | all(ToN=="N")), by="mut_uid")

trunc_filt <- subset(trunc_filt, var_uid %in% VAR$var_uid)
View(subset(trunc_filt, !var_uid %in% VAR$var_uid))

#####################################################################
#### wilcox test
#####################################################################

apply( subset(trunc_filt@VAR, AC>15), 1, function(x){
  print(x["var_uid"]);
  to_test <- subset(GT, var_uid==x["var_uid"]);
  bg <- subset(CBio_query, Gene==x["Gene"]);
  if(nrow(bg)==0){
    print("no such Gene")
  }else{
    print(wilcoxTestCallSet(to_test, "mrnaz", bg))
  }
})

####################################################################
### old
# remove frequent germline truncations in ESP
freq_loh <- c(freq_loh, unique(subset(trunc_filt, !is.na(AA_pos) & !is.na(ESP_AF) & ESP_AF>=0.01  &  AF>=0.01 & 
                                      (abs(AA_pos-AALength)/AALength > 0.2 | (Pfam==1 & AftDom==0)))@VAR$Gene))
trunc_filt <- subset(trunc_filt, !Gene %in% freq_loh) 

trunc_filt@VAR$length <- with(trunc_filt@VAR, abs(nchar(REF)-nchar(ALT)))

View(subset(trunc_filt, length>20)@GT)
#UNC5B very frequent long indel

# 200 most famous cancer gene
# remove some variants near the end
# remove ECT2L, wierd oncogene, high FP
# remove CREBBP, wierd compound mutations
trunc_200 <- subset(trunc_filt, cgc200 & !(var_uid%in% c("CBLC-45296846-A-AC")) & !(Gene %in% c("ECT2L", "CREBBP")))
# clean GT
trunc_200@GT <- subset(trunc_200@GT, SAMPLE_AB>0.15)
#number of Gene
length(unique(trunc_200@VAR$Gene))
#number of Patient
length(unique(trunc_200@GT$Patient))

trunc_200@GT<- merge(trunc_200@GT, CBio_query[c("event_uid", "log2CNA", "gistic", "mrnaz")], all.x=T) 

# compare copy number gain/loss with all patients
tmp<-subset(CBio_query, Gene %in% unique(trunc_200@VAR$Gene))
table(tmp$gistic)
table(trunc_200@GT$gistic)
# patient with truncating variants more likely to have CNV loss
chisq.test(matrix(c(201, 385, 113, 146373, 507784, 140777), nrow=2), simulate.p.value =F)

# look for multiple hit same patient
patient_hit<-dplyr::summarise(group_by(trunc_200@GT, Patient), hit = length(unique(Gene)))
subset(patient_hit, hit>1)
# look for multiple hit same patient/gene
event_hit<-dplyr::summarise(group_by(trunc_200@GT, event_uid), hit = length(unique(var_uid)))
# 26 patient double hit
subset(patient_hit, hit>1)
                    
# intersect with somatic mutation, 13 event, all TS
intersect(trunc_200@GT$event_uid, tcga_mut$event_uid)


                    