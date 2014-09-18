
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

varBurden <-function(NAC, AN, CNTR_AC, CNTR_AN, prefix="CNTR", alt="two.sided"){
  #if(is.na(CNTR_AC) | is.na(CNTR_AN)){
  #  returns<- rep(NA, 3)
  #  names(returns) <-paste(prefix, c(".pval", ".conf.lo", ".conf.hi"), sep="" ) 
  #  return(returns)
  #}else{
  data <- matrix(c(NAC, AN-NAC, CNTR_AC, CNTR_AN) , nrow = 2)
  ftest <- fisher.test(data, alternative = alt)
  #print(names(ftest))
  returns<-c( signif(ftest$p.value, 2), sum(sign(log(ftest$conf.int[is.finite(ftest$conf.int)]))), signif(ftest$conf.int[1], 2), signif(ftest$conf.int[2], 2 ))
  names(returns) <-paste(prefix, c(".pval", ".sign", ".conf.lo", ".conf.hi"), sep="") 
  return(returns)
  #}
}

createVarUid <- function(call_set){
  # check call_set has the info
  if (length(setdiff(c("Gene", "POS", "REF", "ALT"), colnames(call_set))) > 0 ) stop("call_set missing columns")
  call_set$var_uid <- apply(call_set[c("Gene", "POS", "REF", "ALT")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  return(call_set)
}

annoESPX2kG <-function(call_set){
  if(!all(sapply(c("ESP_goi", "X1kG_goi", "X2kG_goi"), exists))){
    load("Results/ESP_X1kG_X2kG_goi.RData")
  }    
  call_set <- call_set[, !colnames(call_set)%in% c("ESP_AC", "ESP_AN", "ESP_AF", "X1kG_AF", "X1kG_ERATE", "X1kG_LDAF", "X2kG_AC", "X2kG_AF", "X2kG_AN")]
  call_set <- merge(call_set, ESP_goi[c("var_uid", "ESP_AC", "ESP_AN", "ESP_AF")], by="var_uid", all.x=T)
  call_set <- merge(call_set, X2kG_goi[c("var_uid", "X2kG_AC", "X2kG_AN", "X2kG_AF")], by="var_uid", all.x=T)
  call_set <- merge(call_set, X1kG_goi[c("var_uid", "X1kG_AF", "X1kG_ERATE", "X1kG_LDAF")], by="var_uid", all.x=T)
  call_set$ESP_AC[is.na(call_set$ESP_AC)] <- 0
  call_set$ESP_AF[is.na(call_set$ESP_AF)] <- 0
  call_set$ESP_AN[is.na(call_set$ESP_AN)] <- 13010
  call_set$X1kG_AF[is.na(call_set$X1kG_AF)] <- 0
  call_set$X1kG_LDAF[is.na(call_set$X1kG_LDAF)] <- 0
  call_set$X2kG_AC[is.na(call_set$X2kG_AC)] <- 0
  call_set$X2kG_AF[is.na(call_set$X2kG_AF)] <- 0
  call_set$X2kG_AN[is.na(call_set$X2kG_AN)] <- 5008
  return(call_set)
}

# hard filtering
hardFilter <- function(x, return_pass = TRUE){
  if (class(x) != "VariantCallSet"){
    stop("object not VariantCallSet")
  }else{
    var_set <- x@VAR
    #SNPS
    low_mq <- with(var_set, (is.na(ALT_idx)| domi_allele ) & ((VARTYPE=="SNP" & MQ <36) | (VARTYPE!="SNP" & MQ <36)))
    high_mqrs <- with(var_set, (is.na(ALT_idx)| domi_allele ) & VARTYPE=="SNP" & !is.na(MQRankSum) & MQRankSum < -12.5)
    high_rprs <- with(var_set, (is.na(ALT_idx)| domi_allele ) &VARTYPE=="SNP" & !is.na(ReadPosRankSum) & ReadPosRankSum < -8)
    high_fs <- with(var_set, (is.na(ALT_idx)| domi_allele ) & ((VARTYPE=="SNP" & FS > 100) | (VARTYPE!="SNP" & FS > 300 ) ))
    #high_fs <- with(var_set, (VARTYPE=="SNP" & FS > 60 & is.na(ALT_idx)) | (VARTYPE!="SNP" & FS > 200 & is.na(ALT_idx)))
    #general
    #low_qd <- with(var_set, (is.na(ALT_idx) | ALT_idx ==1)& QD< 1.5)
    low_qd <- with(var_set, (is.na(ALT_idx)| domi_allele ) & QD< 2)
    low_ad <- with(var_set, AD_max < 3 | (AD_med <3 & AB_max < 0.3))
    low_ab <- with(var_set, (STR_match & STR_times >=8 & AC>20 & AB_med < 0.3) | (AC > 20 & AB_med < 0.25) | AB_med < 0.15  )
    fail_filter <- low_mq | high_mqrs | high_rprs | high_fs | low_qd | low_ad | low_ab
    print(c("number of variants failing filter ", sum(fail_filter)), quote=F)
    if (return_pass){
      x@VAR <- var_set[!fail_filter, ]
    }else{
      x@VAR <- var_set[fail_filter, ]
    }
    x@GT <- subset(x@GT, var_uid %in% x@VAR$var_uid)
    return(x)
    #return(var_set[fail_filter, ])
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

chisqTestCallSet <- function(df, factor_name, bg=all_tcga){
  #check factor_name is in df
  #if( !eval(factor_name) %in% colnames(df)){
  #  stop("factor_name not in df")
  #}
  if( !eval(factor_name) %in% colnames(bg)){
    stop("factor_name not in bg")
  }
  # tally the control distribution
  bg <- bg[c(eval(factor_name), "Patient")]
  # name the factor foi
  colnames(bg)[1] <- "foi"
  bg <- bg[!is.na(bg$foi), ]
  bg_tbl <- dplyr::summarise(group_by(bg, foi), bg_patient = length(unique(Patient)))
  print(c("bg total Patient", sum(bg_tbl$bg_patient)), quote=F)
  # tally the supplied subset distribution
  df2 <- subset(bg, Patient %in% df$Patient)
  colnames(df)[1] <- "foi"
  df_tbl <- dplyr::summarise(group_by(df2, foi), df_patient = length(unique(Patient)))
  print(c("df total Patient", sum(df_tbl$df_patient)), quote=F)
  print(c("lose Patient", length(setdiff(df$Patient, bg$Patient))), quote=F)
  # merge
  conting_tbl <- merge(bg_tbl, df_tbl, all.x=T)
  conting_tbl$df_patient[is.na(conting_tbl$df_patient)] <- 0
  conting_tbl$res_patient <- conting_tbl$bg_patient - conting_tbl$df_patient
  # chisq test
  #print(contig_tbl)
  test <- chisq.test(as.matrix(conting_tbl[c("df_patient", "res_patient")]), simulate.p.value=T, B=2000)
  return(test$p.value)
}    


wilcoxTestCallSet <- function(df, factor_name, bg){
  #check factor_name is in df
  if( !eval(factor_name) %in% colnames(bg)){
    stop("factor_name not in bg")
  }
  # check factor_name is numeric
  if( !is.numeric(bg[[eval(factor_name)]]) ){
    stop("factor_name is not numeric")
  }
  # get the control distribution
  bg <- bg[c(eval(factor_name), "Patient")]
  # name the factor foi
  colnames(bg)[1] <- "foi"
  bg <- bg[!is.na(bg$foi),]
  print(c("bg total Patient", nrow(bg)), quote=F)
  print(c("df total Patient", nrow(bg[bg$Patient %in% df$Patient, ])), quote=F)
  print(c("lose Patient", length(setdiff(df$Patient, bg$Patient))), quote=F)
  # tally the supplied subset distribution
  test <- wilcox.test(subset(bg, Patient %in% df$Patient)$foi, subset(bg, !Patient %in% df$Patient)$foi)
  return(test$p.value)
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
