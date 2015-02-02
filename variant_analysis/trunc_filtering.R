# read vcf using VariantAnnotation interface
trunc_VAR <- expand(readVcf(file="./Results/norm_trunc.VAR.vcf.gz", "hg19"))

trunc_VAR <- trunc_VAR %>% labelUid  %>%   
  subset(., !is.na(Transcript) ) %>%
  subset(., Coding == "CODING" ) %>%
  subset(., is.na(BIOTYPE) | BIOTYPE == "protein_coding" ) %>%
  fixSnpEffGene %>%
  subset(., Gene %in% list_goi$Gene) %>%
  subset(., Impact=="HIGH" ) %>%
  subset(., is.na(SYMBOL) | Gene==SYMBOL) %>% 
  labelVarUid %>% 
  labelSnpEffUid  %>% 
  labelVEPUid(., "trunc") %>% 
  reduceCOL(., "trunc") 


### pick longest unique CCDS transcripts for SnpEFF annotation
trunc_VARu <- pickSnpEffTranscript(trunc_VAR)

# fix LOF_P_Transcripts NAs
info(trunc_VARu)$LOF_P_Transcripts <- replaZero(info(trunc_VARu)$LOF_P_Transcripts)
info(trunc_VARu)$NMD_P_Transcripts <- replaZero(info(trunc_VARu)$NMD_P_Transcripts)


### Simplify VEP Consequence field
trunc_VARu <- trunc_VARu %>% simplifyVEPConseq

### Remove VEP annotation fields that aren't in CCDS, excluding variants that don't have any CCDS annotation
trunc_VARu <- trunc_VARu %>% filterVEPTranscript(., CCDS_15[["nucleotide_ID"]])

### Rank vep_uid according to LoF, only take one vep_uid per var_uid
trunc_VARu <- trunc_VARu %>% pickVEPTranscript

# check vairant effects, remove all variants that contain stop_lost since it's near the end
table(info(trunc_VARu)$EFF)
trunc_VARu <- subset(trunc_VARu, !grepl("stop_lost", info(trunc_VARu)$EFF))

# consolidate SnpEff LoF prediction
info(header(trunc_VARu))["LoF_snpeff",] <- list("1" ,"Flag", "LoF_snpeff")
info(trunc_VARu)$LoF_snpeff <- with(as.data.frame(info(trunc_VARu)[c("NMD_P_Transcripts", "LOF_N")]), NMD_P_Transcripts>0 | LOF_N>0) 

# show contingency table of SnpEff and VEP annotated LoF
with(as.data.frame(info(trunc_VARu)[c("LoF_snpeff", "LoF")]), table(LoF_snpeff, LoF))
info(header(trunc_VARu))["tier1",] <- list("1" ,"Flag", "tier1")
info(header(trunc_VARu))["tier2",] <- list("1" ,"Flag", "tier2")
info(trunc_VARu)$tier1 <- with(as.data.frame(info(trunc_VARu)[c("LoF_snpeff", "LoF")]), LoF_snpeff & LoF!= "LC")
info(trunc_VARu)$tier2 <- with(as.data.frame(info(trunc_VARu)[c("LoF_snpeff", "LoF")]), (LoF_snpeff & LoF== "LC" )| (!LoF_snpeff & LoF=="HC" ))
#remove LCs
#trunc_VAR <- subset(trunc_VAR, LoF!="LC")

#Annotate ESP and X1kG population frequency with annoESP and uid
trunc_VARu <- annoESPX2kG(trunc_VARu)

#manually set POLN 2074702-TG-T  and GEN 17962993-CAAGTT-C to LC
#info(trunc_VARu)["2-17962993-CAAGTT-C","LoF"] <- "LC"
#info(trunc_VARu)["4-2074702-TG-T","LoF"] <- "LC"
# exclude CASP8-202122956-T-C, since it hits only 1/6 CCDS Transcripts.
# exclude GEN1-17962993-CAAGTT-C, last exon and near 10% end.
# "BRCA1-41276044-ACT-A" should be a tier1 , for some reason SnpEff didn't say LoF.

# remove AF>0.99
trunc_VARu <- subset(trunc_VARu, AF<0.99)

# tally the total high confidence variants in 1000G
comm_genes <- as.data.frame(info(trunc_VARu)) %>% subset(., tier1) %>% group_by(., Gene) %>% dplyr::summarise(., tot_X2kG = sum(X2kG_AF)) %>% subset(., tot_X2kG>=0.02) 
comm_genes

# remove these genes
trunc_VARu <- subset(trunc_VARu, !Gene %in% comm_genes$Gene | Gene %in% c("GEN1", "ECT2L", "CASP8"))

# update the GRange with Gene info so as to query CCDS coding region
#dummy_grange <- as.data.frame(rowData(trunc_VARu))
#dummy_grange$Gene <- info(trunc_VARu)$Gene
#dummy_grange$start <- as.numeric(dummy_grange$start)
#dummy_grange$end <- as.numeric(dummy_grange$end)
#info(header(trunc_VARu))["CCDS_hits",] <- list("1" ,"Integer", "CCDS_hits")
#info(trunc_VARu)$CCDS_hits <- with(dummy_grange, mapply(queryCCDS, Gene, start, end))



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
  trunc_GT <- (readVcf(file=paste0("./Results/norm_trunc.", chr ,".GT.vcf.gz", collapse=""), "hg19"))
  #update GT 
  print("filter GT", quote=F)
  trunc_GT <- trunc_GT  %>% labelUid %>% subset(., uid %in% info(trunc_VARu)$uid) 
  print("sumamrise GT", quote=F)
  chr_tally[[chr]] <- summaryVCFGT(trunc_GT)
}
trunc_tally <- list()
trunc_tally$tally <- do.call(rbind, lapply(chr_tally, function(x) x$tally))
# three variants are missing
subset(info(trunc_VARu), !uid%in% trunc_tally$tally$uid)
trunc_VARu <- subset(trunc_VARu, uid %in% trunc_tally$tally$uid)
# get the nonsyn_GT variant summary, EAC, EAN, sample AD and AB, etc
# update nonsyn_VARu with the tallied info
trunc_VARu <- mergeVCF(trunc_VARu, trunc_tally$tally, by="uid")

# for filtering
trunc_VARu <- labelDomiAllele(trunc_VARu)
# quality filter 
trunc_VARu <- hardFilter(trunc_VARu)

# set variants that are in 1kG pass
info(trunc_VARu)$pass[with(as.data.frame(info(trunc_VARu)[c("pass", "X2kG_AC")]), !pass & X2kG_AC>0)] <- T

trunc_VARu <- subset(trunc_VARu, pass)


#####################################################################
#### joint analysis Tumor/Normal
#####################################################################

require(plyr)
# joint call variant file
VAR_files <- list.files(path = "./Results/trunc_fpout", pattern="VAR", recursive=F, full.names=T, )
trunc_VARj <- do.call(plyr::rbind.fill, lapply(VAR_files, function(x) read.delim(x, strip.white=T, stringsAsFactors=F, na.strings = "")))
# joint call GT file
GT_files <- list.files(path = "./Results/trunc_fpout", pattern="GT", recursive=F, full.names=T, )
trunc_GTj <- do.call(rbind, lapply(GT_files, function(x) read.delim(x, strip.white=T, stringsAsFactors=F, na.strings = "")))
trunc_GTj <- plyr::join(trunc_GTj, trunc_VARj[c("ALT_idx", "var_uid")], by="var_uid")
                        
joint_tally <- summaryJointGT(trunc_GTj)

# update trunc_VARj
trunc_VARj <-  plyr::join(trunc_VARj, joint_tally$tally, by="var_uid")

# further filter trunc_VARj
trunc_VARj <-trunc_VARj %>% subset(., TAC!=0 & !Nonly ) %>%
             subset(., QD > 2) %>%
             subset(., QD > 4 | NAC<10 ) %>% # rescue one variant since present in X2kG
             subset(., FS < 60 | VARTYPE!= "SNP") %>%
             subset(., FS < 200) 

# four PEG3 variants seem to form complex variants.
#VAR <- subset(VAR, ReadPosRankSum> -8 | VARTYPE!="SNP") # filters out two PEG3 variants PEG3-57335984-T-A, PEG3-57335976-C-T

#update trunc_VARu
trunc_VARu <- subset(trunc_VARu, var_uid %in% trunc_VARj$var_uid)
trunc_VARu <- mergeVCF(trunc_VARu, joint_tally$tally, by="var_uid")

# visual inspect splice_site_acceptor variants
# "CAV1-116166578-G-A"     "PRKDC-48805815-A-AC"    "TP53INP1-95952448-TC-T" "SLC39A4-145641477-T-C"  "RAD9B-110969377-A-C"    "TINF2-24709362-CTG-C"   "MYBBP1A-4442271-C-A"    "TP53-7578555-C-CT"  
# non LoF variants include, "TP53-7578555-C-CT" , "TINF2-24709362-CTG-C", ?"PRKDC-48805815-A-AC"

# visual inspect splice_site_donor variants
# "KDM1A-23397852-GGT-G"     "PAX3-223066130-C-G"       "CADM2-85985007-TGTAA-T"   "MSH5-31721216-TG-T"       "MAPK15-144799962-GGTGA-G" "KANK1-742405-CG-C"        "APTX-33001564-C-A"        "FAN1-31203019-G-T"       
#" CASC5-40939272-GGTAAA-G"  "FANCI-89849425-TGTGA-T"   "GLTSCR2-48253543-AGT-A"   "LZTR1-21337378-GGTGA-G"
# none LoF include, "LZTR1-21337378-GGTGA-G", "FANCI-89849425-TGTGA-T", "CASC5-40939272-GGTAAA-G", "MAPK15-144799962-GGTGA-G", "CADM2-85985007-TGTAA-T", "KDM1A-23397852-GGT-G"
trunc_VARu <- subset(trunc_VARu, !var_uid %in% c("TP53-7578555-C-CT" , "TINF2-24709362-CTG-C", "LZTR1-21337378-GGTGA-G", "FANCI-89849425-TGTGA-T", 
                                                "CASC5-40939272-GGTAAA-G", "MAPK15-144799962-GGTGA-G", "CADM2-85985007-TGTAA-T", "KDM1A-23397852-GGT-G"))



tmp2 <- dplyr::summarise(group_by(tmp, var_uid), N_LOH = sum(LOH), N_PIT = sum(PIT))
var_tally <- merge(var_tally, tmp2, by="var_uid")





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


                    