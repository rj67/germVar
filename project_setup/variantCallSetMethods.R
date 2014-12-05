
setClass("VariantCallSet", slots=c(GT="data.frame", VAR="data.frame"))

setMethod("show", signature(object="VariantCallSet"),
          function(object) { print(head(object@VAR)); print(head(object@GT))}
)
setMethod("subset", "VariantCallSet",
          function(x, ...) {
            print(c("before nrow ", nrow(x@VAR)), quote=F);
            x@VAR <- subset(x@VAR, ...);
            print(c("after nrow ", nrow(x@VAR)), quote=F);
            x@GT <- subset(x@GT, var_uid %in% x@VAR$var_uid);
            return(x)
          } 
)
setMethod("merge", "VariantCallSet",
          function(x, y, ...) {
            print(c("before nrow ", nrow(x@VAR)), quote=F);
            x@VAR <- merge(x@VAR, y, ...);
            print(c("after nrow ", nrow(x@VAR)), quote=F);
            x@GT <- subset(x@GT, var_uid %in% x@VAR$var_uid);
            return(x)
          } 
)
setMethod("View", "VariantCallSet",
          function(x, ...) {
            View(x@VAR)
          } 
)
setMethod("View", "VCF",
          function(vcf, ...) {
            View(info(vcf))
          } 
)
setMethod("subset", "VCF",
          function(vcf, ...) {
            print(c("before nrow ", nrow(info(vcf))), quote=F);
            vcf <- vcf[with(info(vcf), ...),];
            print(c("after nrow ", nrow(info(vcf))), quote=F);
            return(vcf)
          } 
)
labelUid <- function(vcf){
            if(inherits(vcf, "CollapsedVCF")){
              ALT <- sapply(rowData(vcf)$ALT, function(x) as.character(x[[1]]))
            }else{
              ALT <- as.character(rowData(vcf)$ALT)
            };
            uid <- apply( cbind(as.character(seqnames(rowData(vcf))), start(ranges(rowData(vcf))), as.character(rowData(vcf)$REF), ALT), 1, function(x) paste0(x, collapse="-"));
            names(ranges(rowData(vcf))) <- uid;
            info(header(vcf))["uid",] <- list("1" ,"String", "uid")
            info(vcf)$uid <- uid; 
            return(vcf)
}

labelVarUid <- function(vcf){
  if (inherits(vcf, "VCF")){
    info(header(vcf))["var_uid",] <- list("1" ,"String", "var_uid")
    info(vcf)$var_uid <- apply(cbind(as.character(info(vcf)$Gene), start(ranges(rowData(vcf))), as.character(rowData(vcf)$REF), as.character(rowData(vcf)$ALT)), 1, function(x) paste0(x, collapse="-")) ;
  }else{
    if (!all(c("Gene", "POS", "REF", "ALT") %in% colnames(vcf))) stop("vcf missing columns")
    vcf$var_uid <- apply(vcf[c("Gene", "POS", "REF", "ALT")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  }
  return(vcf)
}

labelVEPUid <- function(vcf){
  info(header(vcf))["vep_uid",] <- list("1" ,"String", "vep_uid")
  AA_change <- sapply(as.character(info(vcf)$HGVSp), function(x) gsub("p.", "", strsplit(x, split=":")[[1]][2], fixed=T ))
  info(vcf)$vep_uid <- apply(cbind(as.character(info(vcf)$var_uid), AA_change), 1, function(x) paste0(x, collapse="-")) ;
  return(vcf)
}

# Gene-AAChange(by SnpEff)
labelAAUid <- function(vcf){
  if (inherits(vcf, "VCF")){  
    info(header(vcf))["aa_uid",] <- list("1" ,"String", "aa_uid")
    info(vcf)$aa_uid <- apply(cbind(as.character(info(vcf)$Gene), as.character(info(vcf)$AAChange)), 1, function(x) paste0(x, collapse="-")) ;
  }else{
    if (!all(c("Gene", "AAChange") %in% colnames(vcf))) stop("labelAAUid: vcf missing columns")
    vcf$aa_uid <- apply(vcf[c("Gene", "AAChange")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  }
  return(vcf)
}

# Gene-AA_pos
labelSite <- function(vcf){
  if (inherits(vcf, "VCF")){
    info(header(vcf))["site",] <- list("1" ,"String", "site")
    info(vcf)$site <- apply(cbind(as.character(info(vcf)$Gene), as.character(info(vcf)$AA_pos)), 1, function(x) paste0(x, collapse="-")) ;
  }else{
    if (!all(c("Gene", "AA_pos") %in% colnames(vcf))) stop("labelSite: vcf missing columns")
    vcf$site <- apply(vcf[c("Gene", "AA_pos")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  }
  return(vcf)
}

# use ensembl_gff to correct Gene names in SnpEff
fixSnpEffGene <- function(vcf){
  if(!exists("table_HGNC")){
    load("Results/table_HGNC.RData")
  }
  if(!exists("ensembl_gff")){
    load("Results/biomart_ensembl_gff.RData")
  }
  if (class(vcf) == "VCF"){
    fix_pos <- !info(vcf)$Gene %in% table_HGNC$Approved.Symbol
    if(sum(fix_pos)>0){
      ensembl <- plyr::rename( ensembl_gff[c("name", "symbol")], c("name"="Transcript"))
      rescued <- plyr::join(as.data.frame(info(vcf)[c("Gene", "Transcript")])[fix_pos,], ensembl, by="Transcript", type="left", match="first")
      rescued$Gene[!is.na(rescued$symbol)] <- rescued$symbol[!is.na(rescued$symbol)]
      info(vcf)$Gene[fix_pos] <- rescued$Gene
    }
  }else{
    fix_pos <- vcf$Gene %in% table_HGNC$Approved.Symbol
    if(sum(fix_pos)>0){
      ensembl <- plyr::rename( ensembl_gff[c("name", "symbol")], c("name"="Transcript"))
      rescued <- plyr::join(vcf[c("Gene", "Transcript")][fix_pos,], ensembl, by="Transcript", type="left", match="first")
      rescued$Gene[!is.na(rescued$symbol)] <- rescued$symbol[!is.na(rescued$symbol)]
      vcf$Gene[fix_pos] <- rescued$Gene
    }
  }
  return(vcf)
}

#extract Gene name from Clinvar record when GeneSymbol field is empty
fixClinGene <- function(vcf){
  Name<-vcf$Name[vcf$GeneSymbol=="-"]
  refGene <- sapply(Name, function(x) {ref<-strsplit(x, split=":", fixed=T)[[1]][1]; 
                                       front<- gregexpr("(", ref, fixed=T)[[1]]; 
                                       end <- gregexpr(")", ref, fixed=T)[[1]]; 
                                      return(ifelse(front>1 & end >1, substr(ref, front+1, end-1), "-")) })
  vcf$GeneSymbol[vcf$GeneSymbol=="-"] <- refGene
  return(vcf)
}

fixVEPCCDS <- function(vcf){
  if(!exists("CCDS_enst")){
    load("Results/CCDS_summary.RData")
  }
  #fix_pos <- is.na(info(vcf)$CCDS) | (! info(vcf)$CCDS %in% CCDS_enst$ccds_id) 
  #if(sum(fix_pos)>0){
  #  rescued <- merge(info(vcf)[c("CCDS", "Feature")][fix_pos,], CCDS_enst[c("ccds_id", "nucleotide_ID")], by.x="Feature", by.y="nucleotide_ID", all.x=T)
  #  rescued$CCDS[!is.na(rescued$ccds_id)] <- rescued$ccds_id[!is.na(rescued$ccds_id)]
  #  info(vcf)$CCDS[fix_pos] <- rescued$CCDS
  #}
  CCDS <- plyr::rename(CCDS_enst[c("ccds_id", "nucleotide_ID")], c("nucleotide_ID"="Feature"))
  print(class(CCDS))
  rescued <- plyr::join(as.data.frame(info(vcf)[c("Feature")]), CCDS, by="Feature", type="left", match="first")
  info(header(vcf))["ccds_id",] <- list("1" ,"String", "ccds_id")
  info(vcf)$ccds_id <- rescued$ccds_id
  return(vcf)
}

# split VEP SIFT score
convSIFT <- function(vcf){
  info(header(vcf))["SIFT_score",] <- list("1" ,"Foat", "SIFT_score")
  info(vcf)$SIFT_score <- sapply(info(vcf)$SIFT, function(x) as.numeric(gsub(")", "", strsplit(x, split="(", fixed=T)[[1]][2], fixed=T)))
  info(vcf)$SIFT <- sapply(info(vcf)$SIFT, function(x) gsub(")", "", strsplit(x, split="(", fixed=T)[[1]][1], fixed=T))
  return(vcf)
}

# based on strand, flip nucleotide sequence
flipCodon <- function(N, strand){
  if(strand==1){
    return(N)
  }else if(strand==-1){
    dict = list("A"="T", "T"="A", "G"="C", "C"="G")
    return(dict[[N]])
  }else{
    stop('invalid strand')
  }  
}  
labelDomiAllele  <- function(vcf){
  info(header(vcf))["domi_allele",] <- list("0" ,"Flag", "Dominant allele")
  info(vcf)$domi_allele <- FALSE
  info(vcf)$domi_allele[!is.na(info(vcf)$ALT_num)] <- with(subset(info(vcf), !is.na(ALT_num)), 
                               mapply(function(ACs, ALT_idx){AC <- as.numeric(strsplit(ACs, split="|", fixed=T)[[1]]); return(AC[as.numeric(ALT_idx)]/sum(AC)>=0.9) },  ACs_orig, ALT_idx))
  return(vcf)
}

getGTPat <- function(uid, vcf){
 Var_Pat <- samples(header(vcf))[grep("1", geno(vcf)$GT[uid,], fixed=T, useBytes=T)] 
 All_Pat <- samples(header(vcf))[!grepl(".", geno(vcf)$GT[uid,], fixed=T, useBytes=T)] 
 Var_pat <- unique(all_tcga$Patient[all_tcga$SM%in%Var_Pat])
 All_pat <- unique(all_tcga$Patient[all_tcga$SM%in%All_Pat])
 return(list(Var_pat=Var_pat, All_pat=All_pat))
}

  vcf<-nonsyn_GT#[1:50,]
  gt_df <- as.data.frame(geno(vcf)$GT)
  gt_df$id <- rownames(gt_df)
  gt_df <- melt(gt_df, id.vars="id", variable.name="SM", value.name="GT")
  ad_df <- as.data.frame(geno(vcf)$AD)
  ad_df$id <- rownames(ad_df)
  ad_df <- melt(ad_df, id.vars="id", variable.name="SM", value.name="GT")
  
sumVCFGT <- function(uid, vcf){
 idx <- info(vcf)[uid, "ALT_idx"][[1]]
 if(idx =="."){
   idx <- 1
 }else{
   idx <- as.numeric(idx)
 }
 SAMPLE_ADs <- geno(vcf)$AD[uid, ][grep("1", geno(vcf)$GT[uid,], fixed=T, useBytes=T)] 
 if(length(SAMPLE_ADs)==0){
   return(data.frame(uid=uid, lowq_AC=0, AD_med=0, AD_max=0, AD_iqr=0, AB_med=0, AB_max=0))
 }else{
   SAMPLE_DPs <- sapply(SAMPLE_ADs, sum)
   print(uid)
   ALT_ADs <- sapply(SAMPLE_ADs, function(x) x[idx+1])
   ALT_ABs <- ALT_ADs/SAMPLE_DPs
   ALT_ABs[is.nan(ALT_ABs)] <- 0
   lowq_AC <- sum( ALT_ADs < 3 | ALT_ABs < 0.15)
   AD_med <- median(ALT_ADs)
   AD_max <- max(ALT_ADs)
   AD_iqr <- IQR(ALT_ADs)
   AB_med <- signif(median(ALT_ABs), 2)
   AB_max <- signif(ALT_ABs[order(ALT_ADs)[length(ALT_ADs)]], 2)
   names(AB_max) <- NULL
   return(data.frame(uid=uid, lowq_AC=lowq_AC, AD_med=AD_med, AD_max=AD_max, AD_iqr=AD_iqr, AB_med=AB_med, AB_max=AB_max))
 }
}


mergeVCF <- function(vcf, df, by="uid"){
  if (inherits(vcf, "VCF")){
    tmp <- plyr::join(as.data.frame(info(vcf)[c(by)]), df, by=by)
    tmp[[by]] <- NULL
    for (coln in colnames(tmp)){
      info(header(vcf))[coln, ] <- list(1, class(tmp[[coln]]), coln)
      info(vcf)[, coln] <- tmp[[coln]]
    }
    return(vcf)
  }else{
    stop("vcf not VCF")
  }
}

varBurden <-function(NAC, AN, CNTR_AC, CNTR_AN, prefix="CNTR", test="t.test",alt="two.sided"){
  if(test=="t.test"){
    library(broom)
    data1 <- c(rep(1, NAC), rep(0, AN-NAC))
    data2 <- c(rep(1, CNTR_AC),rep(0, CNTR_AN-CNTR_AC))
    ttest <- tidy(t.test(data1, data2, alternative = alt))
    #returns<-c( signif(ttest$p.value, 2), sum(sign(ttest$conf.int)), signif(ttest$conf.int[1], 2), signif(ttest$conf.int[2], 2 ))
    returns<-c( signif(ttest$p.value, 2), signif(ttest$conf.low, 2), signif(ttest$conf.high, 2 ))
  }else{
    data <- matrix(c(NAC, AN-NAC, CNTR_AC, CNTR_AN) , nrow = 2)
    ftest <- fisher.test(data, alternative = alt)
    #print(names(ftest))
    returns<-c( signif(ftest$p.value, 2), sum(sign(log(ftest$conf.int[is.finite(ftest$conf.int)]))), signif(ftest$conf.int[1], 2), signif(ftest$conf.int[2], 2 ))
  }
  names(returns) <-paste(prefix, c(".pval", ".conf.lo", ".conf.hi"), sep="") 
  return(returns)
}

writeVCF <- function(vcf, filename){
  if(class(vcf)=="data.frame"){
    write("##fileformat=VCFv4.1", file=filename)
    write('##INFO=<ID=var_uid,Number=1,Type=String,Description="variant uid">', file=filename, append=T)
    write(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), ncolumns=8, file=filename, append=T, sep="\t")
    out_df <- vcf[colnames(vcf) %in% c("CHROM", "POS", "ID", "REF", "ALT", "FILTER", "QUAL")]
    if(!"ID" %in% colnames(vcf)){
      out_df$ID <- "."
    }
    if(!"FILTER" %in% colnames(vcf)){
      out_df$FILTER <- "."
    }
    if(!"QUAL" %in% colnames(vcf)){
      out_df$QUAL <- "."
    }
    out_df$INFO <- paste("var_uid=", vcf$var_uid, sep="")
    out_df <- out_df[c("CHROM", "POS", "ID", "REF", "ALT", "FILTER", "QUAL", "INFO")]
  write.table(arrange(out_df, CHROM, POS), file=filename, col.names=F, row.names=F, quote=F, sep="\t", append=T, na = ".")
  }
}  

annoESPX2kG <-function(vcf){
  if(!all(sapply(c("ESP_goi", "X2kG_goi"), exists))){
    load("DataSet_Results/ESP_X1kG_X2kG_goi.RData")
  }    
  #vcf <- mergeVCF(vcf, ESP_goi[c("var_uid", "ESP_AC", "ESP_AN", "ESP_AF", "ESP_EA_AC", "ESP_EA_AN", "ESP_EA_AF")], by="var_uid")
  vcf <- mergeVCF(vcf, ESP_goi[c("var_uid", "ESP_AC", "ESP_AN", "ESP_AF", "ESP_fAC", "ESP_fAN")], by="var_uid")
  vcf <- mergeVCF(vcf, X2kG_goi[c("var_uid", "X2kG_AC", "X2kG_AN", "X2kG_AF")], by="var_uid")
  info(vcf)$ESP_AC[is.na(info(vcf)$ESP_AC)] <- 0
  info(vcf)$ESP_AF[is.na(info(vcf)$ESP_AF)] <- 0
  info(vcf)$ESP_AN[is.na(info(vcf)$ESP_AN)] <- 13006
  info(vcf)$ESP_fAC[is.na(info(vcf)$ESP_fAC)] <- 0
  #info(vcf)$ESP_EA_AF[is.na(info(vcf)$ESP_EA_AF)] <- 0
  info(vcf)$ESP_fAN[is.na(info(vcf)$ESP_fAN)] <- 9555
  #info(vcf)$ESP_AN[is.na(info(vcf)$ESP_AN)] <- 8600
  info(vcf)$X2kG_AC[is.na(info(vcf)$X2kG_AC)] <- 0
  info(vcf)$X2kG_AF[is.na(info(vcf)$X2kG_AF)] <- 0
  info(vcf)$X2kG_AN[is.na(info(vcf)$X2kG_AN)] <- 5008
  return(vcf)
}

# 
# annoESPX2kG <-function(vcf){
#   if(!all(sapply(c("ESP_goi", "X1kG_goi", "X2kG_goi"), exists))){
#     load("DataSet_Results/ESP_X1kG_X2kG_goi.RData")
#   }    
#   call_set <- call_set[, !colnames(call_set)%in% c("ESP_AC", "ESP_AN", "ESP_AF", "ESP_EA_AC", "ESP_EA_AN", "ESP_EA_AF", "X1kG_AF", "X1kG_ERATE", "X1kG_LDAF", "X2kG_AC", "X2kG_AF", "X2kG_AN")]
#   call_set <- merge(call_set, ESP_goi[c("var_uid", "ESP_AC", "ESP_AN", "ESP_AF", "ESP_EA_AC", "ESP_EA_AN", "ESP_EA_AF")], by="var_uid", all.x=T)
#   call_set <- merge(call_set, X2kG_goi[c("var_uid", "X2kG_AC", "X2kG_AN", "X2kG_AF")], by="var_uid", all.x=T)
#   call_set <- merge(call_set, X1kG_goi[c("var_uid", "X1kG_AF", "X1kG_ERATE", "X1kG_LDAF")], by="var_uid", all.x=T)
#   call_set$ESP_AC[is.na(call_set$ESP_AC)] <- 0
#   call_set$ESP_AF[is.na(call_set$ESP_AF)] <- 0
#   call_set$ESP_AN[is.na(call_set$ESP_AN)] <- 13006
#   call_set$ESP_EA_AC[is.na(call_set$ESP_EA_AC)] <- 0
#   call_set$ESP_EA_AF[is.na(call_set$ESP_EA_AF)] <- 0
#   call_set$ESP_EA_AN[is.na(call_set$ESP_EA_AN)] <- 8600
#   call_set$X1kG_AF[is.na(call_set$X1kG_AF)] <- 0
#   call_set$X1kG_LDAF[is.na(call_set$X1kG_LDAF)] <- 0
#   call_set$X2kG_AC[is.na(call_set$X2kG_AC)] <- 0
#   call_set$X2kG_AF[is.na(call_set$X2kG_AF)] <- 0
#   call_set$X2kG_AN[is.na(call_set$X2kG_AN)] <- 5008
#   return(call_set)
# }

# hard filtering
hardFilter <- function(vcf){
  if (!inherits(vcf, "VCF")){
    stop("object not VCF")
  }else{
    var_set <- as.data.frame(info(vcf))
    #SNPS
    low_mq <- with(var_set, (is.na(ALT_num)| domi_allele ) & ((VARTYPE=="SNP" & MQ <36) | (VARTYPE!="SNP" & MQ <36)))
    high_mqrs <- with(var_set, (is.na(ALT_num)| domi_allele ) & VARTYPE=="SNP" & !is.na(MQRankSum) & MQRankSum < -12.5)
    high_rprs <- with(var_set, (is.na(ALT_num)| domi_allele ) &VARTYPE=="SNP" & !is.na(ReadPosRankSum) & ReadPosRankSum < -8)
    high_fs <- with(var_set, (is.na(ALT_num)| domi_allele ) & ((VARTYPE=="SNP" & FS > 100) | (VARTYPE!="SNP" & FS > 300 ) ))
    high_ic <- with(var_set, (is.na(ALT_num)| domi_allele ) & (abs(InbreedingCoeff) > 0.5 & QD < 8))
    #high_fs <- with(var_set, (VARTYPE=="SNP" & FS > 60 & is.na(ALT_idx)) | (VARTYPE!="SNP" & FS > 200 & is.na(ALT_idx)))
    #general
    #low_qd <- with(var_set, (is.na(ALT_idx) | ALT_idx ==1)& QD< 1.5)
    low_qd <- with(var_set, ((is.na(ALT_num)| domi_allele ) & QD< 2 )| QD < 1)
    low_ad <- with(var_set, AD_max < 3 | (AD_med <3 & AB_max < 0.3))
    low_ab <- with(var_set, (STR_match & STR_times >=8 & AC>20 & AB_med < 0.3) | (AC > 20 & AB_med < 0.25) | AB_med < 0.15  )
    fail_filter <- low_mq | high_mqrs | high_rprs | high_fs | high_ic | low_qd | low_ad | low_ab
    print(c("number of variants failing filter ", sum(fail_filter)), quote=F)
    info(header(vcf))["pass",] <- list("0" ,"Flag", "Whether pass hard filter")
    info(vcf)$pass <- !fail_filter
    return(vcf)
  }
}  

reduceCOL <- function(vcf){
  if (inherits(vcf, "VCF")){
    for (coln in c("CCC", "HWP","VariantType", "HaplotypeScore","BioType","HET","HOM","GTNum", "DB", "DS", "END", "RPA", "CANONICAL", "INTRON", "Amino_acids", "Codons", "DISTANCE",
                   "MNP", "INS", "DEL", "MIXED", "Coding","LoF_info",  "LoF_flags",  "LoF_filter",  "LoF")){
      info(vcf)[[coln]] <- NULL 
    }
    return(vcf)
  }
}


# take variant set and output bed file
writeFPfilter <- function(call_set, filename){
  write.table(call_set@GT[c("mut_uid", "var_uid","CHROM", "POS", "REF", "ALT", "SAMPLE")], file=filename, quote=F, row.names=F, col.names=T, sep="\t")
}

# Out put MAF file
writeMAF <- function(x, filename){
  if (class(x) != "VariantCallSet"){
    stop("object not VariantCallSet")
  }else{
    GT <- merge(x@GT, x@VAR[c("Entrez.Gene", "EFF", "VARTYPE", "ID", "CHROM", "POS","REF", "ALT", "var_uid")], by="var_uid", all=T) 
    # calculate variant length, SNP==0
    GT$length <- abs(nchar(GT$REF) - nchar(GT$ALT))
    
    # trim REF in INS
    if (any(GT$VARTYPE=="INS")) {
      GT$REF[GT$VARTYPE=="INS"]<- "-"
      GT$ALT[GT$VARTYPE=="INS"]<- sapply(GT$ALT[GT$VARTYPE=="INS"], function(x) substr(x, 2, nchar(x)))
      GT$POS[GT$VARTYPE=="INS"]<- GT$POS[GT$VARTYPE=="INS"] + 1
    }
    # trim REF in DEL
    if (any(GT$VARTYPE=="DEL")) {
      GT$REF[GT$VARTYPE=="DEL"]<- sapply(GT$REF[GT$VARTYPE=="DEL"], function(x) substr(x, 2, nchar(x)))
      GT$ALT[GT$VARTYPE=="DEL"]<-  "-"  
      GT$POS[GT$VARTYPE=="DEL"]<- GT$POS[GT$VARTYPE=="DEL"] + 1
    }
    out_df <- data.frame(Hugo_Symbol = GT$Gene, Entrez_Gene_Id = GT$Entrez.Gene, Center = GT$center, Ncbi_Build = 37, Chrom = GT$CHROM, 
                         Start_Position = GT$POS, End_Position = GT$POS+GT$length, Strand = "+", 
                         Variant_Classification = GT$EFF, Variant_Type = GT$VARTYPE,  Reference_Allele = GT$REF,  
                         Norm_Sample_Allele1 = mapply(function(x,y,z)ifelse(x=="1/1", y, z), GT$GT, GT$ALT, GT$REF),
                         Norm_Sample_Allele2 = GT$ALT, 
                         Dbsnp_Rs = GT$ID, Norm_Sample_Barcode = GT$SAMPLE,  Bam_File = GT$filename, Tumor_Sample_UUID = GT$analysis )
    if(missing(filename)){
      return(out_df)
    }else{
      write.table(out_df, filename, quote=F, row.names=F, col.names=F, sep="\t")
    }
  }
}  

# take variant set and output bed file
writeBED <- function(vcf, filename){
  vcf <- vcf[c("CHROM", "POS")]
  vcf$BEG <- vcf$POS -1
  write.table(vcf[c("CHROM", "BEG", "POS")], file=filename, quote=F, row.names=F, col.names=F, sep="\t")
}
# 
# Output VCF file
writeVCF <- function(x, filename){
  write("##fileformat=VCFv4.1", file=filename)
  write('##INFO=<ID=var_uid,Number=1,Type=String,Description="variant uid">', file=filename, append=T)
  write(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), ncolumns=8, file=filename, append=T, sep="\t")
  out_df <- x@VAR[c("CHROM", "POS", "ID", "REF", "ALT")]
  out_df$FILTER <- "."
  out_df$QUAL <- "."
  out_df$INFO <- paste("var_uid=", x@VAR$var_uid, sep="")
  write.table(arrange(out_df, CHROM, POS), file=filename, col.names=F, row.names=F, quote=F, sep="\t", append=T, na = ".")
}  

chisqTestCallSet <- function(df, factor_name, bg=all_tcga, output=T){
  if(is.character(df)){
    Patients <- df
  } else{
    if(!"Patient" %in% colnames(df)){
      stop("supplied dataframe doesn't contain Patient names")
    }else{
      Patients <- df$Patient      
    }
  }
  #check factor_name is in df
  if( !eval(factor_name) %in% colnames(bg)){
    stop("factor_name not in bg")
  }
  # tally the control distribution
  bg <- bg[c(eval(factor_name), "Patient")]
  # name the factor foi
  colnames(bg)[1] <- "fact"
  bg <- bg[!is.nan(bg$fact), ]
  bg_tbl <- dplyr::summarise(group_by(bg, fact), bg_patient = length(unique(Patient)))
  # tally the supplied subset distribution
  df2 <- subset(bg, Patient %in% Patients)
  colnames(df2)[1] <- "fact"
  df_tbl <- dplyr::summarise(group_by(df2, fact), df_patient = length(unique(Patient)))
  # merge
  conting_tbl <- merge(bg_tbl, df_tbl, all.x=T)
  conting_tbl$df_patient[is.na(conting_tbl$df_patient)] <- 0
  conting_tbl$res_patient <- conting_tbl$bg_patient - conting_tbl$df_patient
  if(output){
    print(c("bg total Patient", sum(bg_tbl$bg_patient)), quote=F)
    print(c("df total Patient", sum(df_tbl$df_patient)), quote=F)
    print(c("lose Patient", length(setdiff(Patients, bg$Patient))), quote=F)
    print(conting_tbl[c("df_patient", "res_patient")])
  }
  # chisq test
  #print(contig_tbl)
  test <- chisq.test(as.matrix(conting_tbl[c("df_patient", "res_patient")]), simulate.p.value=T, B=2000)
  return(test$p.value)
}    


oneWayTestCallSet <- function(df, value_name, bg, para=F, output=T){
  if(is.character(df)){
    Patients <- df
  } else{
    if(!"Patient" %in% colnames(df)){
      stop("supplied dataframe doesn't contain Patient names")
    }else{
      Patients <- df$Patient      
    }
  }
  #check value_name is in df
  if( !eval(value_name) %in% colnames(bg)){
    stop("value_name not in bg")
  }
  # check value_name is numeric
  if( !is.numeric(bg[[eval(value_name)]]) ){
    stop("value_name is not numeric")
  }
  # get the total distribution
  bg <- bg[c(eval(value_name), "Patient")]
  # name the value of interest
  colnames(bg)[1] <- "value"
  bg <- subset(bg, !is.na(value))
  if(output){
    print(c("bg total Patient", nrow(bg)), quote=F)
    print(c("df total Patient", nrow(bg[bg$Patient %in% Patients, ])), quote=F)
    print(c("lose Patient", length(setdiff(Patients, bg$Patient))), quote=F)
  }
  # tally the supplied subset distribution
  if(para){
    test <- t.test(subset(bg, Patient %in% Patients)$value, subset(bg, !Patient %in% Patients)$value)
    return(c(test$p.value, test$estimate))
  }else{
    test <- wilcox.test(subset(bg, Patient %in% Patients)$value, subset(bg, !Patient %in% Patients)$value)
    return(test$p.value)
  }
}    


queryCCDS <- function(Gene, POS, REF, ALT){
  require(IRanges)
  if(Gene %in% names(CCDS_IRanges)){
    POS <- as.integer(POS)
    query <- IRanges(POS,  POS+ abs(nchar(REF) - nchar(ALT)))
    subject <- CCDS_IRanges[[Gene]]
    hits <- as.list(findOverlaps(query, maxgap=2, subject))[[1]]
    returns <- ifelse(length(hits)==0, 0, length(unique(names(subject)[hits])))
    return(returns)
  }else{
    return(0)
  }
}



# Merge somatic mutation with germline GT
mergeNormSomaGT <- function(x, soma_mut){
  if (class(x) != "VariantCallSet"){
    stop("object not VariantCallSet")
  }else{
    require(plyr)
    gt_calls <- merge(x@GT, x@VAR[c("var_uid", "Gene", "EFF")], by="var_uid")
    soma_mut <- subset(soma_mut, Patient %in% gt_calls$Patient)
    gt_calls <- rbind.fill(gt_calls, soma_mut)
    return(gt_calls)
    #return(var_set[fail_filter, ])
  }
}  

#tmp <- clinvar_calls@GT
#tmp$log2CNA <- apply(tmp[c("Patient", "Gene")], 1, function(x) ifelse(x[1] %in% rownames(all_cnv), ifelse(x[2] %in% colnames(all_cnv), all_cnv[x[1], x[2]], NA), NA))
#tmp$mrnaz <- apply(tmp[c("Patient", "Gene")], 1, function(x) ifelse(x[1] %in% rownames(all_mrnaz), ifelse(x[2] %in% colnames(all_mrnaz), all_mrnaz[x[1], x[2]], NA), NA))

#var_set <- subset(var_set, (QD >= 1.5&AC<10) | (AC>=10 &QD>=2) )

#load("Results/TCGA_all_mutations.RData")
#tcga_mut
