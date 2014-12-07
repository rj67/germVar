
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
library(VariantAnnotation)
setMethod("View", "VCF",
          function(vcf, ...) {
            View(info(vcf))
          } 
)
setMethod("subset", "VCF",
          function(x, ...) {
            print(c("before nrow ", nrow(info(x))), quote=F);
            x <- x[with(info(x), ...),];
            print(c("after nrow ", nrow(info(x))), quote=F);
            return(x)
          } 
)
# CHROM-POS-REF-ALT
labelUid <- function(vcf){
  if(inherits(vcf, "VCF")){
    if(inherits(vcf, "CollapsedVCF")){
      ALT <- sapply(rowData(vcf)$ALT, function(x) as.character(x[[1]]))
    }else{
      ALT <- as.character(rowData(vcf)$ALT)
    };
    uid <- apply( cbind(as.character(seqnames(rowData(vcf))), start(ranges(rowData(vcf))), as.character(rowData(vcf)$REF), ALT), 1, function(x) paste0(x, collapse="-"));
    names(ranges(rowData(vcf))) <- uid;
    info(header(vcf))["uid",] <- list("1" ,"String", "uid")
    info(vcf)$uid <- uid; 
  }else{
    vcf$uid <- apply(vcf[c("CHROM", "POS", "REF", "ALT")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  }
  return(vcf)
}

# GENE-POS-REF-ALT
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

# extract VEP annotation uid, depending on whether nonsyn or trunc, convert 3 AA letter
labelVEPUid <- function(vcf, dataset="nonsyn"){
  # parse HGVSp for nonsyn
  info(header(vcf))["vep_uid",] <- list("1" ,"String", "vep_uid")
  if(dataset=="nonsyn"){
    info(header(vcf))["vep_aa",] <- list("1" ,"String", "vep_uid")
    AA_change <- sapply(as.character(info(vcf)$HGVSp), function(x) gsub(":", "", strsplit(x, split=":")[[1]][2], fixed=T ))
    # record which rows are actually coding change
    is_p <- sapply(AA_change, function(x) substr(x, 1, 1)=="p") 
    is_p <- is_p &  grepl("missense", info(vcf)$Consequence)
    AA_change <- substr(AA_change, 3, nchar(AA_change))
    # translate VEP 3 letter amino acid into 1 letter
    library(seqinr)
    ref <- a(substr(AA_change[is_p], 1, 3) )
    alt <- a(substr(AA_change[is_p], nchar(AA_change[is_p])-2, nchar(AA_change[is_p])))
    AA_change[is_p] <- paste(paste(ref, substr(AA_change[is_p], 4, nchar(AA_change[is_p])-3), sep=""), alt , sep="")
    info(vcf)$vep_aa <- AA_change
    info(vcf)$vep_uid <- apply(cbind(as.character(info(vcf)$var_uid), AA_change), 1, function(x) paste0(x, collapse="-")) 
  }else if (dataset=="trunc"){
    info(vcf)$vep_uid <- apply(as.data.frame(info(vcf)[c("var_uid", "HGVSc")]), 1, function(x) paste0(x, collapse="-")) 
  }
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
    CCDS <- plyr::rename(CCDS_enst[c("ccds_id", "nucleotide_ID")], c("nucleotide_ID"="Feature"))
    print(class(CCDS))
    rescued <- plyr::join(as.data.frame(info(vcf)[c("Feature")]), CCDS, by="Feature", type="left", match="first")
    info(header(vcf))["ccds_id",] <- list("1" ,"String", "ccds_id")
    info(vcf)$ccds_id <- rescued$ccds_id
    # some CCDS related transcripts not in CCDS_enst, but CCDS contains them
    rescue_var <- as.data.frame(info(vcf)[c("var_uid", "CCDS", "ccds_id")]) %>% group_by(var_uid) %>% 
      dplyr::summarise(., num_ccds_id = length(unique(ccds_id[!is.na(ccds_id)])) ,num_ccds = length(unique(CCDS[!is.na(CCDS)]))) %>%
      subset(., num_ccds>0 & num_ccds_id ==0)
    print("analomous CCDS field")
    print(table(subset(info(vcf), var_uid %in% rescue_var$var_uid)$Gene))
    info(vcf)$ccds_id[info(vcf)$var_uid %in% rescue_var$var_uid] <- info(vcf)$CCDS[info(vcf)$var_uid %in% rescue_var$var_uid]
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

getVarPat <- function(uid, vcf){
  Var_Pat <- samples(header(vcf))[grep("1", geno(vcf)$GT[uid,], fixed=T, useBytes=T)] 
  Var_pat <- unique(all_tcga$Patient[all_tcga$SM%in%Var_Pat])
  return(Var_pat)
}

# summarise variant info and variant patient info from vcf_GT
summaryVCFGT <- function(vcf){
  library(reshape2)
  idx<- as.numeric(gsub(".","1" ,unlist(info(vcf)$ALT_idx), fixed=T))
  ## get AD matrix
  AD_mat <- as.matrix(geno(vcf)$AD)
  # remove colnames to save space
  colnames(AD_mat) <- NULL
  # strip matrix rowwise and stack each row into a matrix
  AD_rows <- apply(AD_mat, 1, function(x) do.call(rbind, c(x)))
  # extract relevant AD using idx
  AD <- mapply(function(AD, idx){ return(cbind(AD[,1], AD[,idx+1] ))}, AD_rows, idx , SIMPLIFY=F)
  AD <- as.data.frame(do.call(rbind, AD))
  colnames(AD) <- c("REF_AD", "ALT_AD")
  #AD$uid <- rep(rownames(vcf), each = dim(vcf)[2])
  #AD$SM <- rep(colnames(vcf), times = dim(vcf)[1])
  
  ## get the Genotype df
  GT <- as.data.frame(geno(vcf)$GT)
  # use row label uid to convert wide to long
  GT$uid <- rownames(GT)
  GT <- melt(GT, id.vars="uid", variable.name="SM", value.name="GT")  
  GT$uid <- factor(GT$uid, levels=rownames(vcf))
  GT$SM <- factor(GT$SM, levels=colnames(vcf))
  
  # arrange GT and AD in same order
  GT <- arrange(GT, uid, SM)
  #AD <- arrange(AD, uid, SM)
  GT <- cbind(GT, AD[c("REF_AD", "ALT_AD")])
  
  ## tally AC
  # remove no  call sample
  GT <- subset(GT, GT!=".")
  # refactor Patient so Variant equal 2/3 
  GT$GT <- as.numeric(factor(GT$GT, levels=c("0/0", "0/1", "1/1")))
  # merge with patient info
  GT <- plyr::join(GT, all_tcga[c("SM", "Patient", "cauca1")], by="SM")
  GT$Patient <- factor(GT$Patient)
  # calculate EAC
  tally1 <- GT %>% group_by(uid) %>% do(calcPat(.))
  
  ## tally Variant AD/AB
  GT <- subset(GT, GT!=1)
  GT <- GT %>% mutate( DP = REF_AD+ALT_AD, AB= ALT_AD/DP)
  GT$AB[is.nan(GT$AB)] <- 0
  tally2 <- GT %>% group_by(uid) %>% dplyr::summarise(
    lowq_AC = sum( ALT_AD < 3 | AB < 0.15),
    AD_med = median(ALT_AD),
    AD_max = max(ALT_AD),
    AD_iqr = IQR(ALT_AD),
    AB_med = signif(median(AB), 2),
    AB_max = signif(AB[order(ALT_AD)[length(ALT_AD)]], 2)
  )
  
  tally <- plyr::join(tally1, tally2, by="uid")
  #some variant dont have any reads, convert NA to 0
  tally[c("lowq_AC", "AD_med", "AD_max", "AD_iqr", "AB_med", "AB_max")] <- sapply(tally[c("lowq_AC", "AD_med", "AD_max", "AD_iqr", "AB_med", "AB_max")], replaZero)
  
  GT <-plyr::join(GT, subset(all_tcga, !duplicated(Patient))[c("Patient", "disease")], by="Patient")
  return(list(tally=tally, GT=GT))
}

# summarise variant info and variant patient info from vcf_GT
summaryJointGT <- function(GT){
  GT$is_var <- with(GT, GT %in% c("./1", "0/1", "1/1"))
  AD <- sapply(GT$AD, function(x) as.numeric(strsplit(x, split=",")[[1]]))
  GT$DP <- sapply(AD, sum)
  GT$DP[is.na(GT$DP)] <- 0
  GT$ALT_AD <- mapply(function(idx, AD) {if(is.na(idx)){return(AD[2])}else{return(AD[idx+1])}}, GT$ALT_idx, AD)
  GT$ALT_AD[is.na(GT$ALT_AD)] <- 0
  GT$AB <-  mapply(function(ALT_AD, DP) {if(DP==0){return(0)}else{return(ALT_AD/DP)}}, GT$ALT_AD, GT$DP)    
  GT <- merge(GT, all_tcga[c("SM", "Patient", "disease", "ToN")], by.x="SAMPLE", by.y="SM")
  GT$mut_uid <- apply(GT[c("Patient", "var_uid")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  #GT$event_uid <- apply(GT[c("Patient", "Gene")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  var_tally <- dplyr::summarise(group_by(GT, var_uid), NAC = length(unique( mut_uid [is_var & ToN =="N"])), TAC = length(unique( mut_uid [is_var & ToN =="T"])) ,
                                   NACq = length(unique( mut_uid [is_var & ToN =="N" & ALT_AD>=3 & AB >= 0.15 ])), 
                                   TACq = length(unique( mut_uid [is_var & ToN =="T" & ALT_AD>=3 & AB >= 0.15 ])) ,
                                   NALT_AD_med = median(ALT_AD[ToN=="N"]), TALT_AD_med = median(ALT_AD[ToN=="T"]),
                                   NAB_med = median(AB[ToN=="N"]), TAB_med = median(AB[ToN=="T"]), Nonly = all(ToN=="N"))
  return(list(tally=var_tally, GT=GT))
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

# output for Mutation Assessor
writeMA <- function(vcf, filename){
  if(inherits(vcf, "VCF")){
    out <- paste("hg19", info(vcf)$uid, sep=",")
    out <- gsub("-", ",", out)
    strands <- sapply(info(vcf)$STRAND, function(x)ifelse(x=="1", "+", "-"))
    #out <- paste(out, strands, sep=" ")
    #out <- paste(out, info(vcf)$var_uid, sep="\t")
    write(out, filename)
  }
}  

annoESPX2kG <-  function(vcf){
  if(!all(sapply(c("ESP_goi", "X2kG_goi"), exists))){
    load("DataSet_Results/ESP_X2kG_goi.RData")
  }    
  #vcf <- mergeVCF(vcf, ESP_goi[c("var_uid", "ESP_AC", "ESP_AN", "ESP_AF", "ESP_EA_AC", "ESP_EA_AN", "ESP_EA_AF")], by="var_uid")
  # merge with ESP
  anno_df <- plyr::join(as.data.frame(info(vcf)[c("uid", "AC")]), ESP_goi[c("uid", "ESP_AC", "ESP_AN", "ESP_AF", "ESP_EA_AN", "ESP_EA_AC", "ESP_AA_AN", "ESP_AA_AC")], by="uid")
  # merge with X2kG
  anno_df <- plyr::join(anno_df, X2kG_goi[c("uid", "X2kG_AC", "X2kG_AN", "X2kG_AF", "EAS_AF", "SAS_AF", "AMR_AF")], by="uid")
  # fix NAs
  anno_df[c("ESP_AC", "ESP_AF", "ESP_EA_AC", "ESP_AA_AC", "X2kG_AC", "X2kG_AF", "EAS_AF", "SAS_AF", "AMR_AF")] <- 
    sapply(anno_df[c("ESP_AC", "ESP_AF", "ESP_EA_AC", "ESP_AA_AC", "X2kG_AC", "X2kG_AF", "EAS_AF", "SAS_AF", "AMR_AF")], replaZero)
  anno_df$X2kG_AN[is.na(anno_df$X2kG_AN)] <- 5008
  anno_df$ESP_AN[is.na(anno_df$ESP_AN)] <- 13006
  anno_df$ESP_EA_AN[is.na(anno_df$ESP_EA_AN)] <- 8600
  anno_df$ESP_AA_AN[is.na(anno_df$ESP_AA_AN)] <- 4406
  # combine ESP_EA with scaled back ESP_AA and X2kG EAS/SAS
  #based on composition EA:4488, AA:455, AS:394, AM:171
  anno_df <- anno_df %>% mutate(
    AS_fAN = ESP_EA_AN/11.4,
    AS_fAC = AS_fAN*(EAS_AF+SAS_AF),
    AM_fAN = ESP_EA_AN/26.2,
    AM_fAC = AM_fAN*AMR_AF,
    AA_fAN = ESP_EA_AN/9.9,
    AA_fAC = AA_fAN*ESP_AA_AC/ESP_AA_AN,
    ESP_fAC = round(ESP_EA_AC + AS_fAC + AM_fAC + AA_fAC), 
    ESP_fAN = round(ESP_EA_AN + AS_fAN + AM_fAN + AA_fAN), 
  )
  vcf <- mergeVCF(vcf, anno_df[c("uid", "ESP_AC", "ESP_AN", "ESP_AF", "ESP_EA_AC","ESP_EA_AN", "X2kG_AC", "X2kG_AF", "ESP_fAC", "ESP_fAN")], by="uid")
  return(vcf)
}


replaZero <- function(x){
  x[is.na(x)] <- 0
  return(x)
}


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

reduceCOL <- function(vcf, dataset="nonsyn"){
  if (inherits(vcf, "VCF")){
    if(dataset=="nonsyn"){
      for (coln in c("CCC", "HWP","VariantType", "HaplotypeScore","BioType","HET","HOM","GTNum", "DB", "DS", "END", "RPA", "CANONICAL", "INTRON", "Amino_acids", "Codons", "DISTANCE",
                     "MNP", "INS", "DEL", "MIXED", "Coding","LoF_info",  "LoF_flags",  "LoF_filter",  "LoF")){
        info(vcf)[[coln]] <- NULL }
    }else if (dataset=="trunc"){
      for (coln in c("CCC", "HWP","VariantType", "HaplotypeScore","BioType","HET","HOM","GTNum", "DB", "DS", "END", "RPA", "CANONICAL", "SIFT", "PolyPhen", "Codons", "DISTANCE", "MIXED", "Coding")){
        info(vcf)[[coln]] <- NULL }
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
