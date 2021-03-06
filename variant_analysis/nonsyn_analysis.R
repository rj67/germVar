nsSNP_VAR <- mergeVCF(nsSNP_VAR, list_goi[c("Gene", "cgc200", "HER", "Cat")], by="Gene")


##########################################
# Low frequency variants in Clinvar
##########################################

lf_clinvar <- as.data.frame(subset(info(nsSNP_VAR),  EAC>=8 & Clinvar=="Patho/Risk" & Gene %in% subset(list_goi, cgc500)$Gene))
lf_clinvar <- arrange(subset(lf_clinvar[c("Gene", "Transcript","AAChange.p", "EAC", "PHRED", "ma_pred", "pred_patho", "RCVaccession")]), pred_patho, Gene, PHRED)
view(subset(lf_clinvar))


##########################################
# LOH
##########################################
LOH_patient <- read.table("Results/ASCAT_out/ploidy_purity.txt")
colnames(LOH_patient) <- c("Patient", "ploidy", "purity")

ascat<- lapply(LOH_patient$Patient, parse_ascat) 
names(ascat) <- LOH_patient$Patient

parse_ascat <- function(Patient){
  ascat_file <- paste("Results/ASCAT_out/", Patient, ".segments_withSize.txt", sep="") 
  ascat_df <- read.csv(ascat_file, strip.white = T)
  ascat_grange <- GRanges( seqnames = Rle(ascat_df$chromosome), 
           ranges = IRanges(start = ascat_df$startBP, end = ascat_df$endBP, names = ascat_df$Segment..))
  return(list(df=ascat_df, grange=ascat_grange))
}
  
query_ascat <- function(Patient, uid){
  mysnp <- GRanges(seqnames=Rle(strsplit(uid, split="-")[[1]][1]), 
                   ranges = IRanges(start = as.integer(strsplit(uid, split="-")[[1]][2]), 
                                      end = as.integer(strsplit(uid, split="-")[[1]][2]), names= "1" ))
  overlap <- findOverlaps(ascat[[Patient]]$grange, mysnp)
  return(subset(ascat[[Patient]]$df, Segment..==queryHits(overlap))[c("chromosome", "nA", "nB", "size", "startBP", "endBP")])
  #print(mysnp)
  #print(subset(ascat[[Patient]]$df, Segment..==queryHits(overlap)))
}
       
  findOverlaps(ascat_grange, mysnp)
  seqinfo = Seqinfo(seqnames=seqinfo(rowData(nonsyn_GT))@seqnames,
                   seqlengths=seqinfo(rowData(nonsyn_GT))@seqlengths))
  SeqInfo(tmp) <- seqinfo(rowData(nonsyn_GT))
}

tmp <- droplevels(subset(nonsyn_tally$GT, Patient %in% LOH_patient$Patient))
tmp <- droplevels(subset(tmp,  uid %in% info(nonsyn_VARu)$uid))
tmp <- merge(tmp, as.data.frame(info(nonsyn_VARu)[c("uid", "Gene", "AAChange", "EAC", "pred_patho", "Clinvar")]), by="uid")
tmp$uid <- as.character(levels(tmp$uid)[tmp$uid])

loh_query <- do.call(rbind,mapply(query_ascat, tmp$Patient, tmp$uid, SIMPLIFY=F))
tmp<-cbind(tmp, loh_query)
# get the seqlength
seqlengths <- as.data.frame(seqinfo(nonsyn_GT))
seqlengths$chromosome <- rownames(seqlengths)
tmp<-merge(tmp, seqlengths[c("chromosome", "seqlengths")], by="chromosome")
tmp$fraction <- tmp$size/tmp$seqlengths

###########################################################################################
# low ferquency variants
nsSNP_lf <- as.data.frame(subset(info(nsSNP_VAR), pass & EAF>=0.001)) %>% plyr::join(., list_goi, by="Gene") 
# rare variants
nsSNP_rv <- as.data.frame(subset(info(nsSNP_VAR), pass & EAF < 0.001)) %>% plyr::join(., list_goi, by="Gene") 

with(nsSNP_lf, table(pred_patho, Clinvar))
with(nsSNP_lf, table(Clinvar))
# manual literature curation
#tmp<- read.csv(file="Results/nonsyn_variant_pathogenicity_curation.csv", header=T, stringsAsFactors=F)
#nsSNP_VAR_lf <- plyr::join(nsSNP_VAR_lf, tmp[c("var_uid","Pathogenicity","Cancer", "Note")], by="var_uid")

# all the predicted deleterious CGC
#view(arrange(subset(nsSNP_VAR_lf,  pred_patho!="Benign" & mem=="CGC"),!HER)[c("Gene", "AAChange", "pred_patho","Origin", "EAC")])

# all the predicted benignCGC
#View(arrange(subset(nsSNP_VAR_lf,  pred_patho!="Benign" & mem=="CGC"),!HER)[c("Gene", "AAChange", "pred_patho","Origin", "EAC", "Pathogenicity", "ClinicalSignificance", "Note")])


View(subset(nsSNP_VAR_lf,  Gene=="BUB1B")[c("uid","Gene", "AAChange", "SIFT", "fathmm_pred","Cscore", "ma_pred","pred_patho","ClinicalSignificance","RCVaccession", "OtherIDs", "EAC", "CAC1", "ESP_AC", "ESP_fAC", "ESP_EA_AC","Cancer", "Note")])

##########################################
# compare AF with ESP&X2kG, 
##########################################
#lf_VAR <- rbind(nsSNP_VAR_lf[c("uid","Gene","AAChange.p", )], 
#                as.data.frame(subset(info(trunc_VARu), EAC>=10)[c("uid", "Gene", "AAChange.p", "CAC1", "CAN1", "ESP_EA_AC", "ESP_EA_AN", "EAC", "EAN", "X2kG_AC")]))
# EA1 with ESP_fAC
nsSNP_fish <- nsSNP_lf[c("uid", "Gene", "AAChange.p", "CAC1", "CAN1", "ESP_EA_AC", "ESP_EA_AN", "EAC", "EAN", "X2kG_AC", "cosm_scount", "Clinvar", "PHRED", "ma_pred")] 
nsSNP_fish <- cbind(nsSNP_fish, do.call(rbind, apply(nsSNP_lf[c("CAC1", "CAN1", "ESP_EA_AC", "ESP_EA_AN")], 1, function(x) varBurden(x, "EA1" ))))
nsSNP_fish$EA1.t.padj <- p.adjust(nsSNP_fish$EA1.t.pval, method="BH") 
nsSNP_fish$EA1.f.padj <- p.adjust(nsSNP_fish$EA1.f.pval, method="BH") 

# compare with X2kG
nsSNP_fish$X2kG_AN <- 5008
nsSNP_fish <- cbind(nsSNP_fish, do.call(rbind, apply(nsSNP_fish[c("EAC", "EAN", "X2kG_AC", "X2kG_AN")], 1, function(x) varBurden(x, "X2kG" ))))
nsSNP_fish$X2kG.t.padj <- p.adjust(nsSNP_fish$X2kG.t.pval, method="BH") 
nsSNP_fish$X2kG.f.padj <- p.adjust(nsSNP_fish$X2kG.f.pval, method="BH") 
 
nsSNP_fish <- plyr::join(nsSNP_fish, nsSNP_lf[c("uid", "Gene", "AAChange.p")])

View(subset(nsSNP_fish, uid %in% c("22-29121087-A-G", "16-2168022-C-A", "8-90990521-T-C")))
to_plot <- subset(nsSNP_lf, uid %in% c("22-29121087-A-G", "16-2168022-C-A", "8-90990521-T-C"))
to_plot <- mutate(to_plot, TCGA_AF = 100*EAF, TCGA_EA_AF = CAC1/CAN1*100,  X2kG_AF = 100*X2kG_AF, ESP_EA_AF = ESP_EA_AC/ESP_EA_AN*100)[c("Gene", "TCGA_AF", "TCGA_EA_AF","X2kG_AF", "ESP_EA_AF")]
to_plot <- melt(to_plot)
to_plot$Gene <- factor(to_plot$Gene, levels=c("CHEK2", "PKD1", "NBN"))
p <- ggplot(to_plot, aes(Gene, value, fill=variable)) + theme_few() + geom_bar(stat="identity", position="dodge", alpha=0.9, width=0.8)
p <- p + xlab("") + ylab("Allele frequency(%)") 
p <- p + scale_x_discrete(labels = c( "CHEK2"="CHEK2/I200T", "PKD1"="PKD1/R324L","NBN"="NBN/I171V"))
p <- p + scale_fill_brewer(palette="RdYlBu", guide = guide_legend(title=NULL), labels=c("TCGA", "TCGA_EA", "1000G", "ESP_EA")) + theme(legend.position=c(0.8, 0.8))
p
ggsave(filename="nsSNP_lf_fisher_AF.png", width=4, height=4)

#lf_table <- subset(nsSNP_VAR_lf, (EA1.padj <0.15 & X2kG.padj<0.2))[c( "EAF", "CAC1", "CAN1", "ESP_EA_AC", "ESP_EA_AN", "X2kG_AF", "Gene", "var_uid", "AAChange", "EA1.padj", "X2kG.padj")]
#lf_table <- lf_table %>% mutate( EAF = signif(EAF, 2), CAF1 = signif(CAC1/CAN1, 2), ESP_EA_AF = signif(ESP_EA_AC/ESP_EA_AN, 2), X2kG_AF = signif(X2kG_AF, 2), EA1.padj = signif(EA1.padj, 2), X2kG.padj = signif(X2kG.padj, 2))
#lf_table <- arrange(lf_table[c("Gene", "AAChange", "EAF", "CAF1", "ESP_EA_AF", "X2kG_AF", "EA1.padj", "X2kG.padj", "var_uid")], EA1.padj)
#lf_table[c("EAF", "CAF1", "ESP_EA_AF", "X2kG_AF")] <- sapply(lf_table[c("EAF", "CAF1", "ESP_EA_AF", "X2kG_AF")], function(x) x*100)

#view(subset(nsSNP_VAR_lf, EA1.padj <0.2 & X2kG.padj<0.2  )[c("AC","EAC", "EAN", "EAF", "CAC1", "CAN1","ESP_AC", "ESP_AN", "ESP_fAC", "ESP_fAN", "ESP_EA_AC", "ESP_EA_AN", "X2kG_AF", "Gene", "uid", "AAChange", 
#                                            "ESPf.pval", "EA1.pval","X2kG.padj", "Clinvar", "SIFT", "Cscore", "fathmm_pred", "RCVaccession", "OtherIDs"  )])

##########################################
# expression change
##########################################

nsSNP_mrnaz <- group_by(subset(nsSNP_GT, uid %in% nsSNP_lf$uid), uid) %>% do(mrnaz_test(.))
nsSNP_mrnaz$padj <- p.adjust(nsSNP_mrnaz$pval , method="BH") 
nsSNP_mrnaz <- plyr::join(nsSNP_mrnaz, nsSNP_lf[c("uid", "Gene", "AAChange.p")])

nsSNP_cna <- group_by(subset(nsSNP_GT, uid %in% nsSNP_lf$uid), uid) %>% do(cna_test(.))
nsSNP_cna$padj <- p.adjust(nsSNP_cna$pval , method="BH") 

# plot 
# save GEN1 RNASeq for nonsyn_analysis.Rmd
ffload("DataSet_Results/RNASeq_summary.ff")
GEN1_RNA <- getRNASeq("GEN1")
GEN1_RNA$VAR <- GEN1_RNA$Patient %in% subset(nsSNP_GT, uid=="2-17962005-C-G")$Patient
OGG1_RNA <- getRNASeq("OGG1")
OGG1_RNA$VAR <- OGG1_RNA$Patient %in% subset(nsSNP_GT, uid=="3-9792107-G-A")$Patient
to_plot<-rbind(GEN1_RNA[c("Gene", "Patient", "mrnaz", "VAR")], OGG1_RNA[c("Gene", "Patient", "mrnaz", "VAR")])
p <- ggplot(aes(Gene, mrnaz), data=to_plot) + geom_boxplot(outlier.shape=NA, width=0.7) + theme_few()
p <- p + geom_jitter(data=subset(to_plot, VAR), aes(color=Gene), position = position_jitter(width = .25), size=3, alpha=0.8) 
p <- p + scale_color_manual(values=c("#f46d43", "#74add1"), labels=c("GEN1/S509W", "OGG1/R46Q"))
p <- p + xlab("") + ylab("Tumor mRNA level(z-score)") + scale_y_continuous(limits=c(-4., 6.5))
p <- p + theme(legend.position=c(0.75, 0.85), legend.title = element_blank())
p
ggsave(filename="nsSNP_GEN1_OGG1_mrnaz_boxplot.png",height=4., width=4)

##########################################
# average age
##########################################
nsSNP_age <- group_by(subset(nsSNP_GT, uid %in% nsSNP_lf$uid), uid) %>% do(age_test(.))
nsSNP_age$padj <- p.adjust(nsSNP_age$pval , method="BH") 


pats <- getGTPat("17-56355397-G-A", nonsyn_GT)
bg <- subset(all_clin, Patient %in% pats$All_pat)
fg <-subset(bg, Patient%in%pats$Var_pat)
p <- ggplot(aes(reorder(study, age, median), age), data=subset(bg, study %in% fg$study)) +geom_boxplot(outlier.shape = NA) + theme_few()
p <- p + geom_jitter( position = position_jitter(width = .2), aes(color = study, size=3), data=fg)
p <- p + scale_color_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(nlevels(fg$study)))
p <- p + theme(axis.text= element_text(angle=90), legend.position="none") + xlab("Study") + ylab("age")
p

##########################################
# Rare variants
##########################################

# missing in tumor
tmp<- read.csv(file="Results/nonsyn_patho_missing_in_tumor.csv", header=T, stringsAsFactors=F)
tmp <- dplyr::summarise(group_by(tmp, var_uid), MAC = length(unique(Patient)))
tmp <- subset(tmp, var_uid %in% nsSNP_rv$var_uid)
nsSNP_rv <- plyr::join(nsSNP_rv, tmp, by="var_uid")
subset(nsSNP_rv, MAC>=EAC*2/3)[c("var_uid", "pred_patho", "EAC", "MAC")]
# remove MAC==EAC
nsSNP_rv <- subset(nsSNP_rv, MAC<EAC*2/3 | is.na(MAC))
nsSNP_rv <- subset(nsSNP_rv, !Gene %in% c("SF3B1", "U2AF1"))

# 
table(subset(nsSNP_rv, pred_patho=="tier1" & driver )$Cat)
table(subset(nsSNP_rv, pred_patho=="tier1" & driver & cosm_acount>=5 )$Cat)
table(subset(nsSNP_rv, pred_patho=="tier1" & driver & cosm_acount< 5 )$Cat)
to_plot <- matrix(c(4, 24, 21, 74), 2)
colnames(to_plot) <- c("")
rownames(to_plot) <- NULL
mosaicplot(to_plot, color="white")

# manual literature curation
#tmp<- read.csv(file="Results/nonsyn_variant_pathogenicity_curation.csv", header=T, stringsAsFactors=F)
#nonsyn_VARu_rr <- plyr::join(nonsyn_VARu_rr, tmp[c("var_uid","Pathogenicity","Cancer", "Note")], by="var_uid")

# all rare variants
gene_stat <- nsSNP_rv  %>% group_by(., Gene) %>% dplyr::summarise(., N_pat=sum(EAC)) %>% subset(N_pat>=8)
nsSNP_rv_mrnaz <- nsSNP_GT %>% subset(., uid %in% nsSNP_rv$uid) %>% subset(., Gene %in% gene_stat$Gene) %>% group_by(., Gene) %>% do(mrnaz_test(.))
# only tier1/2
gene_stat <- nsSNP_rv %>% subset(., pred_patho!="tier3") %>% group_by(., Gene) %>% dplyr::summarise(., N_pat=sum(EAC)) %>% subset(N_pat>=8)
nsSNP_rv_mrnaz <- nsSNP_GT %>% subset(., uid %in% subset(nsSNP_rv, pred_patho=="tier1")$uid) %>% subset(., Gene %in% gene_stat$Gene) %>% group_by(., Gene) %>% do(mrnaz_test(.))

nsSNP_rr_mrnaz$padj <- p.adjust(nsSNP_rr_mrnaz$pval , method="BH") 

nsSNP_rv_cna <- nsSNP_GT %>% subset(., uid %in% subset(nsSNP_rv, pred_patho!="tier3")$uid) %>% subset(., Gene %in% gene_stat$Gene) %>% group_by(., Gene) %>% do(cna_test(.))
nsSNP_rv_cna <- nsSNP_GT %>%   subset(., uid %in% nsSNP_rv$uid) %>% subset(., Gene %in% gene_stat$Gene) %>% group_by(., Gene) %>% do(cna_test(.))


nsSNP_mrnaz <- plyr::join(nsSNP_mrnaz, nsSNP_lf[c("uid", "Gene", "AAChange.p")])

driver_RNA <- lapply(unique(subset(nsSNP_rv, driver)$Gene), getRNASeq)
driver_RNA <- do.call(rbind, driver_RNA)
tmp<-merge(nsSNP_GT, nsSNP_rv[c("uid", "driver", "Cat", "pred_patho")], by="uid")
tmp<-merge(tmp, driver_RNA, by=c("Patient", "Gene"))

driver_CBio <- lapply(unique(subset(nsSNP_rv, driver)$Gene), getCBio)
driver_CBio <- do.call(rbind, driver_CBio)
tmp<-merge(nsSNP_GT, nsSNP_rv[c("uid", "driver", "Cat", "pred_patho")], by="uid")
tmp<-merge(tmp, driver_CBio, by=c("Patient", "Gene"))


# Oncogene
OG_pat <- lapply(subset(nonsyn_VARu_rr,  mem=="CGC" &  pred_patho!="Benign" & Cat=="OG" & (ESP_AC+X2kG_AC<=1))$uid, function(x) {Patient <- getVarPat(x, nonsyn_GT); return(data.frame(uid=rep(x, length(Patient)), Patient=Patient))})
OG_pat <-  merge(do.call(rbind, OG_pat), nonsyn_VARu_rr[c("Gene", "aa_uid", "uid")], by="uid")
genes <- unique(OG_pat$Gene)
OG <- as.data.frame(subset(RNASeq, Gene %in% genes))
OG <- merge(OG_pat, OG, by=c("Patient", "Gene"), all.x=T)

# Tumor suppressor Gene
TSG_pat <- lapply(subset(nonsyn_VARu_rr,  mem=="CGC" &  pred_patho!="Benign" & Cat=="TSG" & (ESP_AC+X2kG_AC<=1))$uid, function(x) {Patient <- getVarPat(x, nonsyn_GT); return(data.frame(uid=rep(x, length(Patient)), Patient=Patient))})
TSG_pat <-  merge(do.call(rbind, TSG_pat), nonsyn_VARu_rr[c("Gene", "aa_uid", "uid")], by="uid")
genes <- unique(TSG_pat$Gene)
TSG <- as.data.frame(subset(CBio_query, Gene %in% genes))
TSG <- merge(TSG_pat, TSG, by=c("Patient", "Gene"), all.x=T)
View(arrange(TSG, Gene))

tmp<-rbind(TSG_pat, OG_pat)

# ret study
RET <- as.data.frame(subset(RNASeq, Gene =="RET"))
subset(RET, Patient%in% subset(OG_pat, Gene=="RET")$Patient)

# test RNASeq expression sd between studies for CGC genes
genes <- unique(nonsyn_VARu_rr$Gene)
all_rna <- as.data.frame(subset(RNASeq, Gene %in% genes))
all_rna$id <- with(all_rna, interaction(study, Gene, drop=T))
all_rna <- subset(all_rna, !duplicated(id))
gene_sd <- arrange(dplyr::summarise(group_by(all_rna, Gene), bet_sd = range(sd)[2]-median(sd)), bet_sd)

# tally per gene
sort(table(nonsyn_VARu_rr$Gene))
# many TP53 variants, 
view(subset(nonsyn_VARu_rr,  Gene=="TP53")[c("uid","var_uid","Gene", "AAChange", "SIFT", "Cscore", "fathmm_pred","ma_pred", "pred_patho", "Pathogenicity","ClinicalSignificance", "EAC", "ESP_AC", "ESP_fAC", "RCVaccession", "OtherIDs","Cancer", "Note")])

#MLH1
View(subset(nonsyn_VARu_rr,  mem=="CGC" & Gene!="TP53" & (Clinvar=="Patho/Risk" | pred_patho!="Benign"))[c("uid","var_uid","Gene", "AAChange", "SIFT", "Cscore", "fathmm_pred","ma_pred", "pred_patho", "Pathogenicity","ClinicalSignificance", "EAC", "ESP_AC", "ESP_fAC", "RCVaccession", "OtherIDs","Cancer", "Note")])

View(subset(nonsyn_VARu_rr,  mem=="CGC" & Gene!="TP53"  & pred_patho!="Benign" & Cat%in%c("OG","TSG") & (ESP_AC+X2kG_AC<=1))[c("uid","var_uid","Gene", "AAChange", "SIFT", "Cscore", "fathmm_pred","ma_pred", "pred_patho", "Pathogenicity","ClinicalSignificance", "EAC", "ESP_AC", "ESP_fAC","X2kG_AC", "RCVaccession", "OtherIDs","Cancer", "Note")])

view(arrange(subset(nonsyn_VARu_rr,  mem=="CGC" & Gene!="TP53"  & pred_patho!="Benign" & Cat%in%c("OG","TSG") & (ESP_AC+X2kG_AC<=1))[c("Gene", "AAChange", "SIFT", "Cscore", "fathmm_pred","ma_pred", "pred_patho", "EAC", "ESP_AC", "X2kG_AC", "Clinvar", "Cat")], Cat, Gene))

View(subset(nonsyn_VARu_rr,  Pathogenicity=="Yes"& mem=="CGC")[c("uid","Gene", "AAChange", "SIFT", "Cscore", "fathmm_pred","ClinicalSignificance", "EAC", "ESP_AC", "ESP_fAC", "Cancer", "Note")])

View(subset(nonsyn_VARu_comm, pred_patho=="Deleterious")[c("uid","var_uid","Gene", "AAChange", "SIFT", "Cscore", "fathmm_pred","ma_pred", "pred_patho", "Pathogenicity","ClinicalSignificance", "EAC", "ESP_AC", "ESP_fAC", "RCVaccession", "OtherIDs","Cancer", "Note")])



View(subset(nonsyn_VARu_rr, mem=="CGC" & HER)[c("AC","EAC", "EAN", "EAF", "CAC1", "CAN1","ESP_AC", "ESP_AN", "ESP_fAC", "ESP_fAN", "ESP_EA_AC", "ESP_EA_AN", "X2kG_AF", "Gene", "uid", "AAChange", 
   "Clinvar", "SIFT", "Cscore", "fathmm_pred", "RCVaccession", "OtherIDs"  )])




testVar <- function(uid, Genex, vcf, data_set, factor="CNA", test="oneWay"){
  print(uid)
  print(Genex)
  pats <- getGTPat(uid, vcf)
  print(sapply(pats, length))
  bg <- as.data.frame(subset(get(data_set), Gene==Genex))
  print(dim(bg))
  bg <- subset(bg, Patient %in% pats$All_pat)
  print(dim(bg))
  return.name <- paste(factor, "pval", sep=".")
  if(test=="oneWay"){
    pval <- oneWayTestCallSet(pats$Var_pat, factor, bg, F, F)
  }else if(test=="chisq"){
    pval <- chisqTestCallSet(pats$Var_pat, factor, bg, F)
  }
  return( pval )
}


testPat <- function(uid, vcf, data_set, factor="CNA", test="oneWay"){
  print(uid)
  pats <- getGTPat(uid, vcf)
  print(sapply(pats, length))
  bg <- as.data.frame(get(data_set) )
  print(dim(bg))
  bg <- subset(bg, Patient %in% pats$All_pat)
  print(dim(bg))
  return.name <- paste(factor, "pval", sep=".")
  if(test=="oneWay"){
    pval <- oneWayTestCallSet(pats$Var_pat, factor, bg, F, F)
  }else if(test=="chisq"){
    pval <- chisqTestCallSet(pats$Var_pat, factor, bg, F)
  }
  return( pval )
}

all_clin2<-all_clin %>%group_by(., study) %>% do(scale_age(.))
scale_age<-function(df){ age<-scale(df$age, center=T, scale=F) ; return(data.frame(Patient=df$Patient, age=age))}

testVar("THBD", "20-23028659-G-A")



tmp<-subset(mut_stat, var_uid %in% subset(nonsyn_calls, patho & ESP_AC==0 & X2kG_AC==0)@VAR$var_uid)
tmp<-as.data.frame(table(tmp$Gene))
colnames(tmp)[1] <- "Gene"
tmp <- merge(tmp, list_goi[c("Gene", "TS", "ONCO")] )
tmp <- arrange(tmp, TS, ONCO, -Freq)
view(tmp)

View(subset(nonsyn_calls, patho & ESP_AC!=0 & X2kG_AC!=0))


tmp<-subset(mut_stat, var_uid %in% subset(nonsyn_calls, patho & NAC<=3)@VAR$var_uid)
length(unique(tmp$Gene))
tmp<-merge(tmp, CBio_query[c("event_uid","gistic2","mrnaz")])
length(unique(tmp$Gene))
table(tmp$gistic2)



tmp2 <- subset(CBio_query, Gene%in% tmp$Gene)
table(tmp2$gistic2)
tmp3 <- rbind(table(tmp$gistic2), table(tmp2$gistic2)/sum(table(tmp2$gistic2))*sum(table(tmp$gistic2)))
rownames(tmp3) <- c("Gene", "Variant")
mosaicplot(tmp3, main="", color=c("pink", "lightgreen", "lightblue"), cex.axis=1.5)

tmp<-(subset(nonsyn_calls, patho & AC>=20)@VAR[c("var_uid","Gene","AAChange", "AC","NAC","TAC",
                                                 "X2kG_AC", "ESP_AC","ESP_AN")])
tmp$ESP.ratio <- signif(with(tmp, NAC/ESP_AC*ESP_AN/15000),2)
tmp$X2kG.ratio <- signif(with(tmp, NAC/X2kG_AC*5008/15000),2)
tmp$geomean<- with(tmp, sqrt(ESP.ratio*X2kG.ratio))
tmp<-arrange(tmp, -geomean)
view(tmp[c("Gene", "AAChange", "NAC", "ESP.ratio", "X2kG.ratio")])

mapply( function(x, y){
  print(x);
  to_test <- subset(nonsyn_common_calls@GT, var_uid==x);
  library(ffbase);
  bg <- subset.ffdf(CBio_query, Gene==y);
  #bg <- subset.ff(CBio_query, Gene==unique(x["Gene"]));
  print(dim(bg))
  #print(chisqTestCallSet(to_test, "mrnaz", bg))
}, nonsyn_common_calls@VAR[1:5,]$var_uid, nonsyn_common_calls@VAR[1:5,]$Gene)
