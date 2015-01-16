
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
#############
# low ferquency variants
nonsyn_VARu_lf <- as.data.frame(subset(info(nonsyn_VARu), pass & EAC>=10))
nonsyn_VARu_lf <- plyr::join(nonsyn_VARu_lf, list_goi, by="Gene") 

lf_genes <- subset(list_goi, Gene%in% nonsyn_VARu_lf$Gene) 
# manual literature curation
tmp<- read.csv(file="Results/nonsyn_variant_pathogenicity_curation.csv", header=T, stringsAsFactors=F)
nonsyn_VARu_lf <- plyr::join(nonsyn_VARu_lf, tmp[c("var_uid","Pathogenicity","Cancer", "Note")], by="var_uid")

# all the predicted deleterious CGC
view(arrange(subset(nonsyn_VARu_lf,  pred_patho!="Benign" & mem=="CGC"),!HER)[c("Gene", "AAChange", "pred_patho","Origin", "EAC")])

# all the predicted benignCGC
View(arrange(subset(nonsyn_VARu_lf,  pred_patho!="Benign" & mem=="CGC"),!HER)[c("Gene", "AAChange", "pred_patho","Origin", "EAC", "Pathogenicity", "ClinicalSignificance", "Note")])


View(subset(nonsyn_VARu_lf,  Gene=="BUB1B")[c("uid","Gene", "AAChange", "SIFT", "fathmm_pred","Cscore", "ma_pred","pred_patho","ClinicalSignificance","RCVaccession", "OtherIDs", "EAC", "CAC1", "ESP_AC", "ESP_fAC", "ESP_EA_AC","Cancer", "Note")])

##########################################
# compare AF with ESP
##########################################
# with ESP_fAC
nonsyn_VARu_lf <- cbind(nonsyn_VARu_lf, as.data.frame(t(apply(nonsyn_VARu_lf[c("EAC", "EAN", "ESP_fAC", "ESP_fAN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "ESPf", "t.test", "greater"))}))))
nonsyn_VARu_lf$ESPf.padj <- p.adjust(nonsyn_VARu_lf$ESPf.pval, method="BH") 
# EA1 with ESP_fAC
nonsyn_VARu_lf <- cbind(nonsyn_VARu_lf, as.data.frame(t(apply(nonsyn_VARu_lf[c("CAC1", "CAN1", "ESP_EA_AC", "ESP_EA_AN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "EA1", "t.test", "greater"))}))))
nonsyn_VARu_lf$EA1.padj <- p.adjust(nonsyn_VARu_lf$EA1.pval, method="BH") 
# EA2 with ESP_fAC
nonsyn_VARu_lf <- cbind(nonsyn_VARu_lf, as.data.frame(t(apply(nonsyn_VARu_lf[c("CAC2", "CAN2", "ESP_EA_AC", "ESP_EA_AN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "EA2", "t.test", "greater"))}))))
nonsyn_VARu_lf$EA2.padj <- p.adjust(nonsyn_VARu_lf$EA2.pval, method="BH") 

# compare with X2kG
nonsyn_VARu_lf$X2kG_AN <- 5008
nonsyn_VARu_lf <- cbind(nonsyn_VARu_lf, as.data.frame(t(apply(nonsyn_VARu_lf[c("EAC", "EAN", "X2kG_AC", "X2kG_AN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "X2kG", "t.test", "greater"))}))))
nonsyn_VARu_lf$X2kG.padj <- p.adjust(nonsyn_VARu_lf$X2kG.pval, method="BH") 
 
lf_table <- subset(nonsyn_VARu_lf, (EA1.padj <0.15 & X2kG.padj<0.2))[c( "EAF", "CAC1", "CAN1", "ESP_EA_AC", "ESP_EA_AN", "X2kG_AF", "Gene", "var_uid", "AAChange", "EA1.padj", "X2kG.padj")]
lf_table <- lf_table %>% mutate( EAF = signif(EAF, 2), CAF1 = signif(CAC1/CAN1, 2), ESP_EA_AF = signif(ESP_EA_AC/ESP_EA_AN, 2), X2kG_AF = signif(X2kG_AF, 2), EA1.padj = signif(EA1.padj, 2), X2kG.padj = signif(X2kG.padj, 2))
lf_table <- arrange(lf_table[c("Gene", "AAChange", "EAF", "CAF1", "ESP_EA_AF", "X2kG_AF", "EA1.padj", "X2kG.padj", "var_uid")], EA1.padj)

view(lf_table)
#view(subset(nonsyn_VARu_lf, EA1.padj <0.2 & X2kG.padj<0.2  )[c("AC","EAC", "EAN", "EAF", "CAC1", "CAN1","ESP_AC", "ESP_AN", "ESP_fAC", "ESP_fAN", "ESP_EA_AC", "ESP_EA_AN", "X2kG_AF", "Gene", "uid", "AAChange", 
#                                            "ESPf.pval", "EA1.pval","X2kG.padj", "Clinvar", "SIFT", "Cscore", "fathmm_pred", "RCVaccession", "OtherIDs"  )])

##########################################
# expression change
##########################################

# test RNASeq mrnaz
nonsyn_VARu_lf$mrnaz.pval <- apply(nonsyn_VARu_lf[c("uid", "Gene")], 1, function(x) {names(x) <- NULL;do.call(testVar, c(as.list(x), nonsyn_GT, "RNASeq", "mrnaz", "oneWay"))})
nonsyn_VARu_lf$mrnaz.padj <- p.adjust(nonsyn_VARu_lf$mrnaz.pval , method="BH") 
# test log2CNA as continous variable
nonsyn_VARu_lf$CNA.pval <- apply(nonsyn_VARu_lf[c("uid", "Gene")], 1, function(x) {names(x) <- NULL;do.call(testVar, c(as.list(x), nonsyn_GT, "CBio_query", "log2CNA", "oneWay"))})
nonsyn_VARu_lf$CNA.padj <- p.adjust(nonsyn_VARu_lf$CNA.pval , method="BH") 
# test disease distribution
nonsyn_VARu_lf$disease.pval <- apply(nonsyn_VARu_lf[c("uid", "Gene")], 1, function(x) {names(x) <- NULL;do.call(testVar, c(as.list(x), nonsyn_GT, "all_tcga", "disease", "chisq"))})
nonsyn_VARu_lf$CNA.padj <- p.adjust(nonsyn_VARu_lf$CNA.pval , method="BH") 
subset(nonsyn_VARu_lf, mrnaz.pval<0.01 )
# save GEN1 RNASeq for nonsyn_analysis.Rmd
GEN1_RNA <- as.data.frame(subset(RNASeq, Gene=="GEN1"))

# age
nonsyn_VARu_lf$age.pval <- sapply(nonsyn_VARu_lf$uid, function(x) testPat(x, nonsyn_GT, "all_clin2", "age", "oneWay"))
nonsyn_VARu_lf$age.padj <- p.adjust(nonsyn_VARu_lf$age.pval , method="BH") 


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

nonsyn_VARu_rr <- as.data.frame(subset(info(nonsyn_VARu), pass & EAC<10))
# missing in tumor
tmp<- read.csv(file="Results/nonsyn_patho_missing_in_tumor.csv", header=T, stringsAsFactors=F)
tmp$aa_uid <- paste(tmp$Gene,  tmp$AAChange, sep="-")
tmp <- dplyr::summarise(group_by(tmp, aa_uid), MAC = length(unique(Patient)))
tmp <- subset(tmp, aa_uid %in% nonsyn_VARu_rr$aa_uid)
nonsyn_VARu_rr <- plyr::join(nonsyn_VARu_rr, tmp, by="aa_uid")
subset(nonsyn_VARu_rr, EAC<MAC)[c("uid","aa_uid", "pred_patho", "EAC", "MAC")]
# remove MAC==EAC
nonsyn_VARu_rr <- subset(nonsyn_VARu_rr, EAC!=MAC | is.na(MAC))
# Gene annotation
nonsyn_VARu_rr <- plyr::join(nonsyn_VARu_rr, list_goi, by="Gene") 
nonsyn_VARu_rr$patho <- with(nonsyn_VARu_rr, interaction(Clinvar, pred_patho, drop=T))
# manual literature curation
tmp<- read.csv(file="Results/nonsyn_variant_pathogenicity_curation.csv", header=T, stringsAsFactors=F)
nonsyn_VARu_rr <- plyr::join(nonsyn_VARu_rr, tmp[c("var_uid","Pathogenicity","Cancer", "Note")], by="var_uid")

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



# test CBio_query
tmp<- read.csv(file="Results/nonsyn_variant_pathogenicity_curation.csv", header=T, stringsAsFactors=F)
to_test <- merge(to_test, tmp[c("var_uid","Pathogenicity","Cancer", "Note")], by="var_uid", all.x=T)
pvals2 <- mapply(testVar, to_test$uid, to_test$Gene)
to_test <- cbind(to_test, as.data.frame(t(pvals)))
to_test$cna_padj <- p.adjust(to_test$cna_pval, method="BH") 
to_test$mrnaz2_padj <- p.adjust(to_test$mrnaz2_pval, method="BH") 

View(subset(to_test, cna_pval<0.1 | mrnaz2_pval<0.05)[c("uid", "aa_uid", "ClinicalSignificance","SIFT", "Cscore","cna_pval","log2_pval", "mrnaz_pval", "mrnaz2_pval")])


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
