# 
# setClass("VariantCallSet", slots=c(GT="data.frame", VAR="data.frame"))
# 
# setMethod("show", signature(object="VariantCallSet"),
#           function(object) { print(head(object@VAR)); print(head(object@GT))}
# )
# setMethod("subset", "VariantCallSet",
#           function(x, ...) {
#             print(c("before nrow ", nrow(x@VAR)), quote=F);
#             x@VAR <- subset(x@VAR, ...);
#             print(c("after nrow ", nrow(x@VAR)), quote=F);
#             x@GT <- subset(x@GT, var_uid %in% x@VAR$var_uid);
#             return(x)
#           } 
# )
# setMethod("merge", "VariantCallSet",
#           function(x, y, ...) {
#             print(c("before nrow ", nrow(x@VAR)), quote=F);
#             x@VAR <- merge(x@VAR, y, ...);
#             print(c("after nrow ", nrow(x@VAR)), quote=F);
#             x@GT <- subset(x@GT, var_uid %in% x@VAR$var_uid);
#             return(x)
#           } 
# )
# setMethod("View", "VariantCallSet",
#           function(x, ...) {
#             View(x@VAR)
#           } 
# )
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
      ALT <- as.character(unlist(rowData(vcf)$ALT))
    }else{
      ALT <- as.character(rowData(vcf)$ALT)
    }
    CHROM <- as.character(seqnames(rowData(vcf)))
    POS <- start(ranges(rowData(vcf))) 
    pos_uid <- paste(CHROM, POS, sep="-")
    uid <- apply( cbind(pos_uid, as.character(rowData(vcf)$REF), ALT), 1, function(x) paste0(x, collapse="-"))
    names(ranges(rowData(vcf))) <- uid
    info(header(vcf))["uid",] <- list("1" ,"String", "uid")
    info(vcf)$uid <- uid 
    info(header(vcf))["pos_uid",] <- list("1" ,"String", "pos_uid")
    info(vcf)$pos_uid <- pos_uid 
    info(header(vcf))["POS",] <- list("1" ,"Integer", "POS")
    info(vcf)$POS <- POS 
    info(header(vcf))["CHROM",] <- list("1" ,"String", "CHROM")
    info(vcf)$CHROM <- CHROM
  }else{
    vcf$uid <- apply(vcf[c("CHROM", "POS", "REF", "ALT")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  }
  return(vcf)
}

# GENE-POS-REF-ALT
labelVarUid <- function(vcf){
  if (inherits(vcf, "VCF")){
    info(header(vcf))["var_uid",] <- list("1" ,"String", "var_uid")
    if(inherits(vcf, "CollapsedVCF")){
      ALT <- as.character(unlist(rowData(vcf)$ALT))
    }else{
      ALT <- as.character(rowData(vcf)$ALT)
    };
    info(vcf)$var_uid <- apply(cbind(as.character(info(vcf)$Gene), start(ranges(rowData(vcf))), as.character(rowData(vcf)$REF), ALT), 1, function(x) paste0(x, collapse="-")) ;
  }else{
    if (!all(c("Gene", "POS", "REF", "ALT") %in% colnames(vcf))) stop("vcf missing columns")
    vcf$var_uid <- apply(vcf[c("Gene", "POS", "REF", "ALT")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  }
  return(vcf)
}

labelSnpEffUid <- function(vcf){
  if (inherits(vcf, "VCF")){
    info(header(vcf))["snpeff_uid",] <- list("1" ,"String", "snpeff_uid")
    info(header(vcf))["AAChange.p",] <- list("1" ,"String", "AAChange.p")
    info(header(vcf))["AAChange.c",] <- list("1" ,"String", "AAChange.c")
    info(vcf)$AAChange.p <- as.character(sapply(info(vcf)$AAChange, function(x) strsplit( strsplit(x, split = '/', fixed = T)[[1]][1], split="p.", fixed =T)[[1]][2]), USE.NAMES=F)
    #names(info(vcf)$AAChange.p ) <- NULL
    info(vcf)$AAChange.c <- as.character(sapply(info(vcf)$AAChange, function(x) strsplit(x, split = 'c.', fixed = T)[[1]][2]), USE.NAMES=F)
    #names(info(vcf)$AAChange.c ) <- NULL
    #info(vcf)$snpeff_uid <- apply(as.data.frame(info(vcf)[c("Transcript", "AAChange.c")]), 1, function(x) paste0(x, collapse="-"))
    info(vcf)$snpeff_uid <- apply(as.data.frame(info(vcf)[c("Transcript", "var_uid")]), 1, function(x) paste0(x, collapse="-"))
    # Amino Acid position
    AA_pos <- sapply(info(vcf)$AAChange.p, function(x) ifelse(grepl("_", x), strsplit(x, split="_")[[1]][1], x), USE.NAMES=F)
    AA_pos <- sapply(AA_pos, function(x) {query<-regexpr("[0-9]+", x); substr(x, query[1], query[1] + attributes(query)$match.length - 1)}, USE.NAMES=F)
    info(header(vcf))["AA_pos",] <- list("1" ,"Integer", "AA_pos")
    info(vcf)$AA_pos <- as.integer(AA_pos)
    info(header(vcf))["AA_frac",] <- list("1" ,"Float", "AA_frac")
    info(vcf)$AA_frac <- round(info(vcf)$AA_pos/info(vcf)$AALength, 3)
    # Nucleotide position
    CDS_pos <- sapply(info(vcf)$AAChange.c, function(x) ifelse(grepl("_", x), strsplit(x, split="_")[[1]][1], x), USE.NAMES=F)
    CDS_pos <- sapply(CDS_pos, function(x) {query<-regexpr("[0-9]+", x); substr(x, query[1], query[1] + attributes(query)$match.length - 1)}, USE.NAMES=F)
    info(header(vcf))["CDS_pos",] <- list("1" ,"Integer", "CDS_pos")
    info(vcf)$CDS_pos <- as.integer(CDS_pos)
    info(header(vcf))["CDS_frac",] <- list("1" ,"Float", "CDS_frac")
    info(vcf)$CDS_frac <- round(info(vcf)$CDS_pos/info(vcf)$AALength/3, 3)
  }
  return(vcf)
}
# 
# # extract VEP annotation uid, depending on whether nonsyn or trunc, convert 3 AA letter
# labelVEPUid <- function(vcf, dataset="nonsyn"){
#   # parse HGVSp for nonsyn
#   info(header(vcf))["vep_uid",] <- list("1" ,"String", "vep_uid")
#   if(dataset=="nonsyn"){
#     info(header(vcf))["vep_aa",] <- list("1" ,"String", "vep_uid")
#     AA_change <- sapply(as.character(info(vcf)$HGVSp), function(x) gsub(":", "", strsplit(x, split=":")[[1]][2], fixed=T ))
#     # record which rows are actually coding change
#     is_p <- sapply(AA_change, function(x) substr(x, 1, 1)=="p") 
#     is_p <- is_p &  grepl("missense", info(vcf)$Consequence)
#     AA_change <- substr(AA_change, 3, nchar(AA_change))
#     # translate VEP 3 letter amino acid into 1 letter
#     library(seqinr)
#     ref <- a(substr(AA_change[is_p], 1, 3) )
#     alt <- a(substr(AA_change[is_p], nchar(AA_change[is_p])-2, nchar(AA_change[is_p])))
#     AA_change[is_p] <- paste(paste(ref, substr(AA_change[is_p], 4, nchar(AA_change[is_p])-3), sep=""), alt , sep="")
#     info(vcf)$vep_aa <- AA_change
#     info(vcf)$vep_uid <- apply(cbind(as.character(info(vcf)$var_uid), AA_change), 1, function(x) paste0(x, collapse="-")) 
#   }else if (dataset=="trunc"){
#     info(vcf)$vep_uid <- apply(as.data.frame(info(vcf)[c("var_uid", "HGVSc")]), 1, function(x) paste0(x, collapse="-")) 
#   }
#   return(vcf)
# }

# Gene-AAChange(by SnpEff)
labelAAUid <- function(vcf){
  if (inherits(vcf, "VCF")){  
    info(header(vcf))["aa_uid",] <- list("1" ,"String", "aa_uid")
    info(vcf)$aa_uid <- apply(cbind(as.character(info(vcf)$Gene), as.character(info(vcf)$AAChange.p)), 1, function(x) paste0(x, collapse="-")) ;
  }else{
    if (!all(c("Gene", "AAChange.p") %in% colnames(vcf))) stop("labelAAUid: vcf missing columns")
    vcf$aa_uid <- apply(vcf[c("Gene", "AAChange.p")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
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

# 
# fixVEPCCDS <- function(vcf){
#     if(!exists("CCDS_enst")){
#       load("Results/CCDS_summary.RData")
#     }
#     CCDS <- plyr::rename(CCDS_enst[c("ccds_id", "nucleotide_ID")], c("nucleotide_ID"="Feature"))
#     print(class(CCDS))
#     rescued <- plyr::join(as.data.frame(info(vcf)[c("Feature")]), CCDS, by="Feature", type="left", match="first")
#     info(header(vcf))["ccds_id",] <- list("1" ,"String", "ccds_id")
#     info(vcf)$ccds_id <- rescued$ccds_id
#     # some CCDS related transcripts not in CCDS_enst, but CCDS contains them
#     rescue_var <- as.data.frame(info(vcf)[c("var_uid", "CCDS", "ccds_id")]) %>% group_by(var_uid) %>% 
#       dplyr::summarise(., num_ccds_id = length(unique(ccds_id[!is.na(ccds_id)])) ,num_ccds = length(unique(CCDS[!is.na(CCDS)]))) %>%
#       subset(., num_ccds>0 & num_ccds_id ==0)
#     print("analomous CCDS field")
#     print(table(subset(info(vcf), var_uid %in% rescue_var$var_uid)$Gene))
#     info(vcf)$ccds_id[info(vcf)$var_uid %in% rescue_var$var_uid] <- info(vcf)$CCDS[info(vcf)$var_uid %in% rescue_var$var_uid]
#     return(vcf)
# }

# VEP splice site FP
splice_loftee <- with(info(vcf), grepl("NON_CAN_SPLICE", LoF_filter) | grepl("EXON_INTRON_UNDEF", LoF_filter) |(grepl("ANC_ALLELE", LoF_filter) & grepl("splice", EFF)))

### annotate the effect of LoF variant
annoLoF <- function(vcf){
  # store var_uid
  ids <- unique(info(vcf)$var_uid)
   # tally HC LOF variants
  info(header(vcf))["dist_3edge",] <- list("1" ,"Integer", "dist_3edge")
  info(vcf)$dist_3edge <- with(info(vcf), ifelse(strand==1, exon_chrom_end-POS, POS-exon_chrom_start))
  info(header(vcf))["dist_5edge",] <- list("1" ,"Integer", "dist_5edge")
  info(vcf)$dist_5edge <- with(info(vcf), ifelse(strand==1, POS-exon_chrom_start, exon_chrom_end-POS))
  # SnpEff has bug, insertion more than 2bp away from the reverse strand donor/forward strand acceptor site is predicted LoF, even though they shouldn't affect
  splice_fp <- with(info(vcf), (dist_3edge < -2|dist_5edge < -2) & length>0)
  # remove donor sites that insert GT at the beginning
  insert_seq <- sapply(info(vcf)$AAChange.c, function(x) strsplit(x, split="ins")[[1]][2], USE.NAMES=F)
  # figure out the relevant part of sequence
  ins_fp <- mapply(splice_edge, info(vcf)$EFF, info(vcf)$length, info(vcf)$strand, info(vcf)$dist_5edge, info(vcf)$dist_3edge, insert_seq )
  # splice variant that has involve noncoding exons, these are splice junction between 5UTR and first coding exons or 3UTR with last exon 
  splice_utr <- with(info(vcf), grepl("splice_", EFF) & is.na(cds_start))
  # additionally, if splice acceptor and cds_start is 1, there might be some distance between splice junction and start codon, lots of room for alternative splice junction.
  splice_met <- with(info(vcf), grepl("splice_acceptor", EFF) & CDS_pos==1 )
  # combine the splice filter
  splice_fail <- splice_fp | splice_utr 
  # predict NMD
  coding_nmd <- with(info(vcf), (grepl("frameshift_variant", EFF)| grepl("stop_gained", EFF)) & ((TotalExon - ExonRank)>2 | ((TotalExon - ExonRank)==1 & edge_dist>50)))
  # splice nmd has to exclude splice_met
  splice_nmd <- with(info(vcf), grepl("splice", EFF) & !splice_met & !splice_fail & (TotalExon - ExonRank)>2 )
  # flag for high confidence LoF
  info(header(vcf))["LoF_HC",] <- list("1" ,"Logical", "LoF_HC")
  info(vcf)$LoF_HC <- with(info(vcf), ( CDS_frac < 0.95 | coding_nmd | splice_nmd ) & !splice_fail)
  info(header(vcf))["NMD_HC",] <- list("1" ,"Logical", "NMD_HC")
  info(vcf)$NMD_HC <- coding_nmd | splice_nmd
  info(header(vcf))["splice_fail",] <- list("1" ,"Logical", "splice_fail")
  info(vcf)$splice_fail <- splice_fail
  info(header(vcf))["ins_fp",] <- list("1" ,"Logical", "ins_fp")
  info(vcf)$ins_fp <- ins_fp
  return(vcf)
}

splice_edge <- function(EFF, length, strand, dist_5edge, dist_3dge, insert_seq){
  if(grepl("splice", EFF) & length>0){
    if(grepl("splice_donor", EFF)){
      if(strand==1){ # donor positive strand
        if(dist_3dge==0){
          return(substr(insert_seq, 1, min(2, nchar(insert_seq)))=="GT")
        }else if(dist_3dge==-1){
          return(substr(insert_seq, 1, 1)=="T")
        }
      }else{ # donor negative strand
        if(dist_3dge==-2){
          return(substr(insert_seq, 1, 1)=="T")
        }else if(dist_3edge < -2){
          return(T)
        }else{
          return(F)
        }
      }
    }else if (grepl("splice_acceptor", EFF)){
      if(strand==1){ # acceptor positive strand
        if(dist_5dge==-2){
          return(substr(insert_seq, nchar(insert_seq), nchar(insert_seq))=="A" & substr(insert_seq, 1, 1) !="G" & !grepl("AG", insert_seq))
        }else if(dist_5edge < -2){
          return(!grepl("AG", insert_seq))
        }else{
          return(F)
        }
      }else{ # acceptor negative strand
        if(dist_5edge==-1){
          return(substr(insert_seq, nchar(insert_seq), nchar(insert_seq))=="A" & substr(insert_seq, 1, 1)!="G" & !grepl("AG", insert_seq))
        }else{
          return(F)
        }
      }
    }
  }else{
    return(F)
  }
}


### Summarise LoF, retain the longest transcript.
summariseLoF <- function(vcf){
  # store var_uid
  ids <- unique(info(vcf)$var_uid)
  # record tier1 variants
  tier1 <- as.data.frame(info(vcf)[c("var_uid", "longest", "LoF_HC")]) %>% subset(., longest & LoF_HC)
  print(nrow(tier1))
  tier2 <- as.data.frame(info(vcf)[c("var_uid", "LoF_HC")]) %>% subset(., LoF_HC &  ( !var_uid %in% tier1$var_uid)) 
  # remove redundant 
  uniq_snpeff <- as.data.frame(info(vcf)[c("var_uid", "longest", "snpeff_uid")]) %>% arrange(., !longest) %>%  
    subset(., !duplicated(var_uid)) %>% dplyr::select(., matches("snpeff_uid"))
  info(vcf)$keep <- info(vcf)$snpeff_uid %in% uniq_snpeff[["snpeff_uid"]]
  vcf <- subset(vcf, keep)
  info(vcf)$keep <- NULL
  vcf <- subset(vcf, !duplicated(var_uid))
  # flag for high confidence LoF
  info(header(vcf))["tier1",] <- list("1" ,"Logical", "tier1")
  info(vcf)$tier1 <- info(vcf)$var_uid %in% tier1$var_uid
  info(header(vcf))["tier2",] <- list("1" ,"Logical", "tier2")
  info(vcf)$tier2 <- info(vcf)$var_uid %in% tier2$var_uid
  # check didn't lose LOF var_uid
  loss <- setdiff(ids, unique(info(vcf)$var_uid))
  print(c("Number of variants lost in summariseLoF:", length(loss)), quote=F)
  return(vcf)
}

### pick longest unique CCDS transcripts for SnpEFF annotation
pickSnpEffTranscript <- function(vcf, CCDS){
  # store var_uid
  ids <- unique(info(vcf)$var_uid)
  # merge CCDS Transcript with vcf
  vcf <- mergeVCF(vcf, list_gff_CCDS[c("Transcript", "cds_length", "longest")], by="Transcript")
  # tally HC LOF variants, AA_frac < 0.95 | ExonTotal - ExonRank
  #LOF_tally <- as.data.frame(info(vcf)[c("var_uid", "ccds_id", "CDS_frac", "ExonRank", "ExonTotal")]) %>%
  #  group_by(., var_uid) %>% dplyr::summarise(., LOF_N = length(unique(ccds_id[CDS_frac<0.95 | ExonTotal-ExonRank >=2])))
  # merge with vcf
  #vcf <- mergeVCF(vcf, LOF_tally, by="var_uid")
  # rank snpeff_uid according to AALength, only take one snpeff_uid per var_uid
  uniq_snpeff <- as.data.frame(info(vcf)[c("var_uid", "cds_length", "longest", "snpeff_uid")]) %>% arrange(., !longest, -cds_length) %>%  
    subset(., !duplicated(var_uid)) %>% dplyr::select(., matches("snpeff_uid"))
  info(vcf)$keep <- info(vcf)$snpeff_uid %in% uniq_snpeff[["snpeff_uid"]]
  vcf <- subset(vcf, keep)
  info(vcf)$keep <- NULL
  vcf <- subset(vcf, !duplicated(var_uid))
  # check didn't lose LOF var_uid
  loss <- setdiff(ids, unique(info(vcf)$var_uid))
  print(c("Number of variants lost in pickSnpEffTranscript:", length(loss)), quote=F)
  return(vcf)
}

# ### Simplify VEP Consequence field, Remove VEP annotations with no Consequence, excluding variants that don't have any Consequence annotation
# simplifyVEPConseq <- function(vcf){
#   # store var_uid
#   ids <- unique(info(vcf)$var_uid)
#   retains <- c("frameshift_variant", "stop_lost", "initiator_codon_variant", "splice_acceptor_variant", "splice_donor_variant", "stop_gained" )
#   info(vcf)$Consequence <- sapply(info(vcf)$Consequence, function(x) paste0(intersect(strsplit(x, split="&")[[1]], retains), collapse="&"))
#   # variants where all the Consequence are NULL
#   no_conseq<-as.data.frame(info(vcf)[c("var_uid", "Consequence")]) %>% group_by(var_uid) %>% dplyr::summarise(., num_conseq = sum(Consequence!="", na.rm=T))  %>% subset(., num_conseq==0)
#   # subset vcf
#   info(vcf)$keep <- info(vcf)$var_uid %in% no_conseq[["var_uid"]] | info(vcf)$Consequence != ""
#   vcf <- subset(vcf, keep)
#   info(vcf)$keep <- NULL
#   # check didn't lose LOF var_uid
#   loss <- setdiff(ids, unique(info(vcf)$var_uid))
#   print(c("Number of variants lost in simplifyVEPConseq:", length(loss)), quote=F)
#   return(vcf)
# }

### Remove VEP annotation fields that aren't in CCDS, excluding variants that don't have any CCDS annotation
# For each variant, tally the number of VEP annotated fields that are in CCDS
# filterVEPTranscript <- function(vcf, Transcripts){
#   # store var_uid
#   ids <- unique(info(vcf)$var_uid)
#   no_ccds<-as.data.frame(info(vcf)[c("var_uid", "Feature")]) %>% group_by(var_uid) %>% dplyr::summarise(., num_ccds = sum(Feature %in% Transcripts, na.rm=T))  %>% subset(., num_ccds==0)
#   print(c("Number of variants with no CCDS Transcripts:", nrow(no_ccds)), quote=F)
#   info(vcf)$keep <- info(vcf)$Feature %in% Transcripts | info(vcf)$var_uid %in% no_ccds[["var_uid"]]
#   vcf <- vcf %>% subset(.,  keep)
#   info(vcf)$keep <- NULL
#   loss <- setdiff(ids, unique(info(vcf)$var_uid))
#   print(c("Number of variants lost in filterVEPTranscript:", length(loss)), quote=F)
#   return(vcf)
# }
# 
# ### Rank vep_uid according to LoF, only take one vep_uid per var_uid
# pickVEPTranscript <- function(vcf){
#   # store var_uid
#   ids <- unique(info(vcf)$var_uid)
#   # for ANC_ALLELE filter, check ANC_Codon is the same as ALT_Codon
#   info(vcf)$LoF[with(info(vcf), LoF=="LC" & LoF_filter=="ANC_ALLELE" & !is.na(FunClass) & FunClass=="NONSENSE" & !is.na(ALT_Codon) & !is.na(ANC_Codon) & (ALT_Codon != ANC_Codon) )] <- "HC"
#   # set LoF NA equal unknown
#   info(vcf)$LoF[is.na(info(vcf)$LoF)] <- "unknown"
#   info(vcf)$LoF <- factor(info(vcf)$LoF, levels=c("HC", "LC", "unknown"))
#   #extract the position info and exon
#   info(header(vcf))["T_EXON",] <- list("1" ,"Integer", "T_EXON")
#   info(vcf)$T_EXON <- sapply(info(vcf)$EXON, function(x) as.numeric(strsplit(x, split="/", fixed=T)[[1]][2]))
#   info(header(vcf))["T_INTRON",] <- list("1" ,"Integer", "T_INTRON")
#   info(vcf)$T_INTRON <- sapply(info(vcf)$INTRON, function(x) as.numeric(strsplit(x, split="/", fixed=T)[[1]][2]))
#   # rank vep_uid by LoF and then by number of EXON and INTRON
#   uniq_vep <- as.data.frame(info(vcf)[c("var_uid", "LoF", "T_EXON", "T_INTRON", "vep_uid")]) %>% arrange(., LoF) %>% arrange(., -T_EXON) %>% arrange(., -T_INTRON) %>% subset(., !duplicated(var_uid)) %>% dplyr::select(., matches("vep_uid"))
#   # get the uniq var_uid set
#   info(vcf)$keep <- info(vcf)$vep_uid %in% uniq_vep[["vep_uid"]]
#   vcf <- vcf %>% subset(.,  keep) %>% subset(., !duplicated(var_uid))
#   info(vcf)$keep <- NULL
#   # check didn't lose LOF var_uid
#   loss <- setdiff(ids, unique(info(vcf)$var_uid))
#   print(c("Number of variants lost in pickVEPTranscript:", length(loss)), quote=F)
#   return(vcf)
# }
# 
# # split VEP SIFT score
# convSIFT <- function(vcf){
#   info(header(vcf))["SIFT_score",] <- list("1" ,"Foat", "SIFT_score")
#   info(vcf)$SIFT_score <- sapply(info(vcf)$SIFT, function(x) as.numeric(gsub(")", "", strsplit(x, split="(", fixed=T)[[1]][2], fixed=T)))
#   info(vcf)$SIFT <- sapply(info(vcf)$SIFT, function(x) gsub(")", "", strsplit(x, split="(", fixed=T)[[1]][1], fixed=T))
#   return(vcf)
# }

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
  info(vcf)$domi_allele[!is.na(info(vcf)$ALT_num)] <- with(as.data.frame(info(subset(vcf, !is.na(ALT_num)))), 
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

calcPat <- function(df){
  df <- arrange(df, -GT)
  df <- df[!duplicated(df$Patient), ]
  tally <- data.frame(AC2 = as.integer(sum(df$GT-1)), EAC = as.integer(sum(df$GT[df$EA]-1)) )
  #tally <- data.frame(AC2 = as.integer(sum(df$GT-1)),  AN2 = as.integer(nrow(df)*2),
  #                    EAC = as.integer(sum(df$GT[df$EA]-1)), EAN = as.integer(nrow(subset(df, EA))*2 ))
  #tally$AF2 <- tally$AC2/tally$AN2
  #tally$EAF <- tally$EAC/tally$EAN
  return(tally)
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
  #AD_rows <- apply(AD_mat, 1, function(x) do.call(rbind, c(x)))
  AD_rows <- lapply(split(AD_mat, seq(nrow(AD_mat))), function(x) do.call(rbind, c(x)))
  # extract relevant AD using idx
  AD <- mapply(function(AD, idx){ return(cbind(AD[,1], AD[,idx+1] ))}, AD_rows, idx , SIMPLIFY=F)
  AD <- as.data.frame(do.call(rbind, AD))
  colnames(AD) <- c("REF_AD", "ALT_AD")
  AD$REF_AD <- as.numeric(AD$REF_AD)
  AD$ALT_AD <- as.numeric(AD$ALT_AD)
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
  GT <- plyr::join(GT, all_tcga[c("SM", "Patient", "EA", "suspect")], by="SM")
  GT$Patient <- factor(GT$Patient)
  # calculate AC2 etc, remove suspect normal samples
  tally1 <- GT %>% subset(., !suspect) %>% group_by(uid) %>% do(calcPat(.))
  
  ## tally Variant AD/AB
  GT <- subset(GT, GT!=1)
  GT <- GT %>% mutate( DP = REF_AD+ALT_AD, AB= ALT_AD/DP)
  GT$AB[is.nan(GT$AB)] <- 0
  tally2 <- GT %>% subset(., !suspect) %>% group_by(uid) %>% dplyr::summarise(
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

summaryVCFAN <- function(vcf){
  pos_uid <- paste(as.character(seqnames(rowData(vcf))), start(ranges(rowData(vcf))), sep="-")
  names(ranges(rowData(vcf))) <- pos_uid
  ## get the Genotype df
  GT <- as.data.frame(geno(vcf)$GT)
  # use row label uid to convert wide to long
  GT$pos_uid <- rownames(GT)
  GT <- melt(GT, id.vars="pos_uid", variable.name="SM", value.name="GT")  
  GT$pos_uid <- factor(GT$pos_uid, levels=rownames(vcf))
  GT$SM <- factor(GT$SM, levels=colnames(vcf))
  # arrange GT in same order, remove no call, merge with all_tcga
  GT %<>% arrange(., pos_uid, SM) %>% subset(., GT!=".") %>% plyr::join(., all_tcga[c("SM", "Patient", "EA", "suspect")], by="SM")
  # calculate the AN per pos_uid, and also for EA sub population
  tally <- GT %>% subset(., !suspect) %>% group_by(pos_uid) %>% dplyr::summarise(., AN2 = length(unique(Patient))*2, EAN = length(unique(Patient[EA]))*2 )
  return(tally)
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
  GT <- merge(GT, all_tcga[c("SM", "Patient", "disease", "ToN", "sample_type", "analyte", "suspect","Specimen")], by.x="SAMPLE", by.y="SM")
  # remove the Patients with no Normal sample
  GT <- subset(GT, Patient %in% subset(all_tcga, ToN=="N" & !suspect)$Patient)
  GT$mut_uid <- apply(GT[c("Patient", "uid")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  #GT$event_uid <- apply(GT[c("Patient", "Gene")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
  var_tally <- dplyr::summarise(group_by(GT, uid), NAC = length(unique( mut_uid [is_var & ToN =="N"])), TAC = length(unique( mut_uid [is_var & ToN =="T"])) ,
                                   NACq = length(unique( mut_uid [is_var & ToN =="N" & ALT_AD>=3 & AB >= 0.15 ])), 
                                   TACq = length(unique( mut_uid [is_var & ToN =="T" & ALT_AD>=3 & AB >= 0.15 ])) ,
                                   NALT_AD_med = median(ALT_AD[ToN=="N"]), TALT_AD_med = median(ALT_AD[ToN=="T" & is_var]),
                                   NAB_med = median(AB[ToN=="N"]), TAB_med = median(AB[ToN=="T" & is_var]), Nonly = all(ToN=="N"))
  MUT <- dplyr::summarise(group_by(GT, uid, Patient, mut_uid), N_het = any(GT=="0/1"& (DP - ALT_AD) >= 3 & ToN=="N"), 
                   T_hom = any( AB>=0.7 & DP >= 16), LOH = N_het & T_hom, PIT = any(is_var & ToN=="T"), N_DIS = length(unique(is_var[ToN=="N"]))>1, T_DIS = length(unique(is_var[ToN=="T"]))>1)
  return(list(tally = var_tally, GT = GT, MUT = MUT))
}

mergeVCF <- function(vcf, df, by="uid"){
  if (inherits(vcf, "VCF")){
    tmp <- plyr::join(as.data.frame(info(vcf)[c(by)]), df, by=by)
    if (length(by)>1){
      for (coln in by){
        tmp[[coln]] <- NULL
      }
    }else{
      tmp[[by]] <- NULL
    }
    for (coln in colnames(tmp)){
      info(header(vcf))[coln, ] <- list(1, class(tmp[[coln]]), coln)
      info(vcf)[, coln] <- tmp[[coln]]
    }
    return(vcf)
  }else{
    stop("vcf not VCF")
  }
}

# varBurden <-function(NAC, AN, CNTR_AC, CNTR_AN, prefix="CNTR", test="t.test",alt="two.sided"){
#   if(test=="t.test"){
#     library(broom)
#     data1 <- c(rep(1, NAC), rep(0, AN-NAC))
#     data2 <- c(rep(1, CNTR_AC),rep(0, CNTR_AN-CNTR_AC))
#     ttest <- tidy(t.test(data1, data2, alternative = alt))
#     #returns<-c( signif(ttest$p.value, 2), sum(sign(ttest$conf.int)), signif(ttest$conf.int[1], 2), signif(ttest$conf.int[2], 2 ))
#     returns<-c( signif(ttest$p.value, 2), signif(ttest$conf.low, 2), signif(ttest$conf.high, 2 ))
#   }else{
#     data <- matrix(c(NAC, AN-NAC, CNTR_AC, CNTR_AN) , nrow = 2)
#     ftest <- fisher.test(data, alternative = alt)
#     #print(names(ftest))
#     returns<-c( signif(ftest$p.value, 2), sum(sign(log(ftest$conf.int[is.finite(ftest$conf.int)]))), signif(ftest$conf.int[1], 2), signif(ftest$conf.int[2], 2 ))
#   }
#   names(returns) <-paste(prefix, c(".pval", ".conf.lo", ".conf.hi"), sep="") 
#   return(returns)
# }

varBurden <-function(x, prefix="CNTR", alt="greater"){
    library(broom)
    NAC <- x[1]
    AN <- x[2]
    CNTR_AC <- x[3]
    CNTR_AN <- x[4]
    ### first do a ttest
    data1 <- c(rep(1, NAC), rep(0, AN-NAC))
    data2 <- c(rep(1, CNTR_AC),rep(0, CNTR_AN-CNTR_AC))
    ttest <- tidy(t.test(data1, data2, alternative = alt))
    ### also perform fisher exact test
    data <- matrix(c(NAC, AN-NAC, CNTR_AC, CNTR_AN) , nrow = 2)
    ftest <- fisher.test(data, alternative = alt)
    #returns<-c( signif(ttest$p.value, 2),  signif(ftest$p.value, 2), sum(sign(log(ftest$conf.int[is.finite(ftest$conf.int)]))), signif(ftest$conf.int[1], 2), signif(ftest$conf.int[2], 2 ))
    returns<- lapply(c( ttest$p.value, ftest$p.value, ftest$estimate, ftest$conf.int[1], ftest$conf.int[2]), function(x) signif(x, 2))
    returns <- as.data.frame(returns)
    colnames(returns) <-paste(prefix, c(".t.pval", ".f.pval", ".odds", ".conf.lo", ".conf.hi"), sep="") 
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
    strands <- sapply(info(vcf)$strand, function(x)ifelse(x==1, "+", "-"))
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
  anno_df <- plyr::join(subset(as.data.frame(info(vcf)[c("uid", "AC")]), !duplicated(uid)), subset(ESP_goi, !duplicated(uid))[c("uid", "ESP_AC", "ESP_AN", "ESP_AF", "ESP_EA_AN", "ESP_EA_AC", "ESP_AA_AN", "ESP_AA_AC")], by="uid")
  # merge with X2kG
  anno_df <- plyr::join(anno_df, subset(X2kG_goi, !duplicated(uid))[c("uid", "X2kG_AC", "X2kG_AN", "X2kG_AF", "EAS_AF", "SAS_AF", "AMR_AF")], by="uid")
  # fix NAs
  anno_df[c("ESP_AC", "ESP_AF", "ESP_EA_AC", "ESP_AA_AC", "X2kG_AC", "X2kG_AF", "EAS_AF", "SAS_AF", "AMR_AF")] <- 
    sapply(anno_df[c("ESP_AC", "ESP_AF", "ESP_EA_AC", "ESP_AA_AC", "X2kG_AC", "X2kG_AF", "EAS_AF", "SAS_AF", "AMR_AF")], replaZero)
  anno_df$X2kG_AN[is.na(anno_df$X2kG_AN)] <- 5008
  anno_df$ESP_AN[is.na(anno_df$ESP_AN)] <- 13006
  anno_df$ESP_EA_AN[is.na(anno_df$ESP_EA_AN)] <- 8600
  anno_df$ESP_AA_AN[is.na(anno_df$ESP_AA_AN)] <- 4406
  # combine ESP_EA with scaled back ESP_AA and X2kG EAS/SAS
  #based on composition EA:4488, AA:455, AS:394, AM:171
#   anno_df <- anno_df %>% mutate(
#     AS_fAN = ESP_EA_AN/11.4,
#     AS_fAC = AS_fAN*(EAS_AF+SAS_AF),
#     AM_fAN = ESP_EA_AN/26.2,
#     AM_fAC = AM_fAN*AMR_AF,
#     AA_fAN = ESP_EA_AN/9.9,
#     AA_fAC = AA_fAN*ESP_AA_AC/ESP_AA_AN,
#     ESP_fAC = round(ESP_EA_AC + AS_fAC + AM_fAC + AA_fAC), 
#     ESP_fAN = round(ESP_EA_AN + AS_fAN + AM_fAN + AA_fAN), 
#   )
  vcf <- mergeVCF(vcf, anno_df[c("uid", "ESP_AC", "ESP_AN", "ESP_AF", "ESP_EA_AC","ESP_EA_AN", "X2kG_AC", "X2kG_AF")], by="uid")
  return(vcf)
}


replaZero <- function(x){
  x[is.na(x)] <- 0
  return(x)
}

sigSymbol <- function(pval) {
  return(do.call(c, sapply(pval, function(x) ifelse(x>0.05, "", ifelse(x>0.01, "*", ifelse(x>0.001, "**", "***"))), simplify=F, USE.NAMES=F)))
}

# hard filtering
hardFilter <- function(vcf, vartype="SNP"){
  if (!inherits(vcf, "VCF")){
    stop("object not VCF")
  }else{
    var_set <- as.data.frame(info(vcf))
    #SNPS
    low_mq <- with(var_set, (is.na(ALT_num)| domi_allele ) & ((VARTYPE=="SNP" & MQ <36) | (VARTYPE!="SNP" & MQ <36)))
    high_mqrs <- with(var_set, (is.na(ALT_num)| domi_allele ) & VARTYPE=="SNP" & !is.na(MQRankSum) & MQRankSum < -12.5)
    high_rprs <- with(var_set, (is.na(ALT_num)| domi_allele ) &VARTYPE=="SNP" & !is.na(ReadPosRankSum) & ReadPosRankSum < -8)
    high_fs <- with(var_set, (is.na(ALT_num)| domi_allele ) & ((VARTYPE=="SNP" & SOR > 5) | (VARTYPE!="SNP" & SOR > 5 ) ))
    high_ic <- with(var_set, (is.na(ALT_num)| domi_allele ) & (abs(InbreedingCoeff) > 0.6 & AF2>=0.01 & CHROM!="X"))
    #high_fs <- with(var_set, (VARTYPE=="SNP" & FS > 60 & is.na(ALT_idx)) | (VARTYPE!="SNP" & FS > 200 & is.na(ALT_idx)))
    #general
    #low_qd <- with(var_set, (is.na(ALT_idx) | ALT_idx ==1)& QD< 1.5)
    low_qd <- with(var_set, (is.na(ALT_num)| domi_allele ) & QD< 2 )
    #low_qd <- with(var_set, QD < 2)
    low_ad <- with(var_set, AD_max < 3 )
    if (vartype=="LoF"){
      low_ab <- with(var_set, (STR_match & STR_times >=8 & AC2>20 & AB_med < 0.3) | (AC2 > 20 & AB_med < 0.25) | AB_med < 0.15  )
    }else if(vartype=="SNP"){
      low_ab <- with(var_set, (AF2 > 0.001 & AB_med < 0.2) | AB_max < 0.15  )
    }
    fail_filter <- cbind(low_mq, high_mqrs, high_rprs, high_fs, high_ic, low_qd, low_ad, low_ab)
    print(c("number of variants failing filter ", sum(apply(fail_filter, 1,  any))), quote=F)
    info(header(vcf))["pass",] <- list("0" ,"Flag", "Whether pass hard filter")
    info(vcf)$pass <-  !apply(fail_filter, 1,  any)
    info(header(vcf))["fail_flags",] <- list("0" ,"String", "what filter failed")
    info(vcf)$fail_flags <- apply(fail_filter, 1,  function(x) paste0(names(x)[x], collapse="|"))
    return(vcf)
  }
}  

reduceCOL <- function(vcf, dataset="nsSNP"){
  if (inherits(vcf, "VCF")){
    if(dataset=="nsSNP"){
      for (coln in c("CCC", "HWP","VariantType", "HaplotypeScore","BioType","HET","HOM","GTNum", "DB", "DS", "END", "RPA", "CANONICAL", "INTRON", "Amino_acids", "Codons", "DISTANCE",
                     "MNP", "INS", "DEL", "MIXED", "Coding","LOF_Gene", "LOF_NT", "LOF_PT", "NMD_Gene", "NMD_N_Transcripts", "NMD_P_Transcripts","LoF_info",  "LoF_flags",  "LoF_filter",  "LoF")){
        info(vcf)[[coln]] <- NULL }
    }else if (dataset=="LoF"){
      for (coln in c("CCC", "HWP","VariantType", "HaplotypeScore", "MQ0", "BioType", "Coding", "FunClass","HET","HOM","GTNum", "Warning","DB", "DS", "END", "RPA", 
          "CANONICAL", "SIFT", "PolyPhen", "Codons", "DISTANCE", "MNP","MIXED", "NMD_N_Transcripts", "NMD_P_Transcripts")){
        info(vcf)[[coln]] <- NULL }
    }
    return(vcf)
  }
}


# take variant set and output bed file
# writeFPfilter <- function(call_set, filename){
#   write.table(call_set@GT[c("mut_uid", "var_uid","CHROM", "POS", "REF", "ALT", "SAMPLE")], file=filename, quote=F, row.names=F, col.names=T, sep="\t")
# }

# Out put MAF file
# writeMAF <- function(x, filename){
#   if (class(x) != "VariantCallSet"){
#     stop("object not VariantCallSet")
#   }else{
#     GT <- merge(x@GT, x@VAR[c("Entrez.Gene", "EFF", "VARTYPE", "ID", "CHROM", "POS","REF", "ALT", "var_uid")], by="var_uid", all=T) 
#     # calculate variant length, SNP==0
#     GT$length <- abs(nchar(GT$REF) - nchar(GT$ALT))
#     
#     # trim REF in INS
#     if (any(GT$VARTYPE=="INS")) {
#       GT$REF[GT$VARTYPE=="INS"]<- "-"
#       GT$ALT[GT$VARTYPE=="INS"]<- sapply(GT$ALT[GT$VARTYPE=="INS"], function(x) substr(x, 2, nchar(x)))
#       GT$POS[GT$VARTYPE=="INS"]<- GT$POS[GT$VARTYPE=="INS"] + 1
#     }
#     # trim REF in DEL
#     if (any(GT$VARTYPE=="DEL")) {
#       GT$REF[GT$VARTYPE=="DEL"]<- sapply(GT$REF[GT$VARTYPE=="DEL"], function(x) substr(x, 2, nchar(x)))
#       GT$ALT[GT$VARTYPE=="DEL"]<-  "-"  
#       GT$POS[GT$VARTYPE=="DEL"]<- GT$POS[GT$VARTYPE=="DEL"] + 1
#     }
#     out_df <- data.frame(Hugo_Symbol = GT$Gene, Entrez_Gene_Id = GT$Entrez.Gene, Center = GT$center, Ncbi_Build = 37, Chrom = GT$CHROM, 
#                          Start_Position = GT$POS, End_Position = GT$POS+GT$length, Strand = "+", 
#                          Variant_Classification = GT$EFF, Variant_Type = GT$VARTYPE,  Reference_Allele = GT$REF,  
#                          Norm_Sample_Allele1 = mapply(function(x,y,z)ifelse(x=="1/1", y, z), GT$GT, GT$ALT, GT$REF),
#                          Norm_Sample_Allele2 = GT$ALT, 
#                          Dbsnp_Rs = GT$ID, Norm_Sample_Barcode = GT$SAMPLE,  Bam_File = GT$filename, Tumor_Sample_UUID = GT$analysis )
#     if(missing(filename)){
#       return(out_df)
#     }else{
#       write.table(out_df, filename, quote=F, row.names=F, col.names=F, sep="\t")
#     }
#   }
# }  


# 
# # Output VCF file
# writeVCF <- function(x, filename){
#   write("##fileformat=VCFv4.1", file=filename)
#   write('##INFO=<ID=var_uid,Number=1,Type=String,Description="variant uid">', file=filename, append=T)
#   write(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), ncolumns=8, file=filename, append=T, sep="\t")
#   out_df <- x@VAR[c("CHROM", "POS", "ID", "REF", "ALT")]
#   out_df$FILTER <- "."
#   out_df$QUAL <- "."
#   out_df$INFO <- paste("var_uid=", x@VAR$var_uid, sep="")
#   write.table(arrange(out_df, CHROM, POS), file=filename, col.names=F, row.names=F, quote=F, sep="\t", append=T, na = ".")
# }  

ViewDup <- function(df, value){
  dups <- df[[eval(value)]][duplicated(df[[eval(value)]])]
  View(df[df[[eval(value)]] %in% dups, ])
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

getRNASeq <- function(x){
  is.error <- function(x) inherits(x, "try-error")
  df <- try(subset(RNASeq, Gene==x))
  if(is.error(df)){
    return(NULL)
  }else{
    return(as.data.frame(df))
  }
}

getCBio <- function(x){
  is.error <- function(x) inherits(x, "try-error")
  df <- try(subset(CBio_query, Gene==x))
  if(is.error(df)){
    return(NULL)
  }else{
    return( subset(as.data.frame(df), Patient %in% all_tcga$Patient))
  }
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

# 
# queryCCDS <- function(Gene, start, end){
#   require(IRanges)
#   if(Gene %in% names(CCDS_IRanges)){
#     #print(c(Gene, start, end))
#     query <- IRanges(start,  end)
#     subject <- CCDS_IRanges[[Gene]]
#     hits <- as.list(findOverlaps(query, maxgap=2, subject))[[1]]
#     returns <- ifelse(length(hits)==0, 0, length(unique(names(subject)[hits])))
#     return(returns)
#   }else{
#     return(0)
#   }
# }


# 
# # Merge somatic mutation with germline GT
# mergeNormSomaGT <- function(x, soma_mut){
#   if (class(x) != "VariantCallSet"){
#     stop("object not VariantCallSet")
#   }else{
#     require(plyr)
#     gt_calls <- merge(x@GT, x@VAR[c("var_uid", "Gene", "EFF")], by="var_uid")
#     soma_mut <- subset(soma_mut, Patient %in% gt_calls$Patient)
#     gt_calls <- rbind.fill(gt_calls, soma_mut)
#     return(gt_calls)
#     #return(var_set[fail_filter, ])
#   }
# }  

########################################################################################################################
###           Analysis tools
########################################################################################################################
mrnaz_test <- function(df){
  print(unique(df$Gene))
  bg <- getRNASeq(unique(df$Gene))  
  if( is.null(bg)){
    return( data.frame(N_all = N_all, N_var = 0, effect = 0, std = 0, pval = NA))
  }else{
    bg <- subset(bg, !is.na(mrnaz) & med!=0)
    #bg <- getCBio(unique(df$Gene))  
    bg <- subset(bg, !is.na(mrnaz) & Patient %in% all_tcga$Patient)
    # remove NAs and med==0 tissues
    N_all <- nrow(bg)
    bg$VAR <- bg$Patient %in% df$Patient
    N_var <- sum(bg$VAR)
    if( N_var > 0){
      lm.mod <- lm( mrnaz ~  VAR, bg) 
      # get average age diff
      effect <- broom::tidy(lm.mod)$estimate[broom::tidy(lm.mod)$term=="VARTRUE"]
      std <- broom::tidy(lm.mod)$std.error[broom::tidy(lm.mod)$term=="VARTRUE"]
      test <- wilcox.test(subset(bg, VAR)$mrnaz, subset(bg, !VAR)$mrnaz)#, alternative = "less")
      return( data.frame(N_all = N_all, N_var = N_var, effect = effect, std = std, pval=test$p.value))
    }else{
      return( data.frame(N_all = N_all, N_var = 0, effect = 0, std = 0, pval = NA))
    }
  }
}

cna_test <- function(df){
  bg <- getCBio(unique(df$Gene))  
  bg <- subset(bg, !is.na(gistic) & Patient %in% all_tcga$Patient)
  bg$gistic2 <- cut(bg$gistic, breaks=c(-Inf, -1.5, 2.5, Inf), labels=c("-", "0", "+"))
  # remove NAs and med==0 tissues
  N_all <- nrow(bg)
  bg$VAR <- bg$Patient %in% df$Patient
  N_var <- sum(bg$VAR)
  if( N_var > 0){
    bg_tbl <- dplyr::summarise(group_by(bg, gistic2), bg_patient = length(unique(Patient)))
    df_tbl <- dplyr::summarise(group_by(subset(bg, VAR), gistic2), df_patient = length(unique(Patient)))    
    conting_tbl <- merge(bg_tbl, df_tbl, all.x=T)
    conting_tbl$df_patient <- replaZero(conting_tbl$df_patient)
    conting_tbl$res_patient <- conting_tbl$bg_patient - conting_tbl$df_patient
    test <- chisq.test(as.matrix(conting_tbl[c("df_patient", "res_patient")]), simulate.p.value=T, B=2000)
    #return( data.frame(N_all = N_all, N_var = N_var, res.loss= test$residuals[1,1], res.gain= test$residuals[3,1], pval = test$p.value))
    return( data.frame(N_all = N_all, N_var = N_var, pval = test$p.value))
  }else{
    return( data.frame(N_all = N_all, N_var = 0, pval = NA))
  }
}

plot_RNASeq <- function(df){
  library(ggthemes)
  library(RColorBrewer)
  bg <- getRNASeq(unique(df$Gene))
  Patients <- df$Patient
  # drop diseases with no variant
  bg <- droplevels(subset(bg, study %in% droplevels(subset(bg, Patient %in% Patients))$study))
  p <- ggplot(aes(reorder(study, normalized_count, median), sqrt(normalized_count)), data = bg) + geom_boxplot(outlier.shape = NA)
  p <- p + geom_jitter(aes(color=study), data = subset(bg, Patient %in% Patients))
  p <- p + theme_few()
  p <- p + scale_color_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(length(unique(bg$study))), guide="none")
  p <- p + xlab("study") + ylab("sqrt(Normalized count)")
  p 
}

age_test <- function(df){
  bg <- subset(all_tcga[c("Patient", "agez")], !is.na(agez) & !duplicated(Patient))
  # remove NAs and med==0 tissues
  N_all <- nrow(bg)
  bg$VAR <- bg$Patient %in% df$Patient
  N_var <- sum(bg$VAR)
  if( N_var > 0){
    lm.mod <- lm( agez ~  VAR, bg) 
    # get average age diff
    effect <- broom::tidy(lm.mod)$estimate[broom::tidy(lm.mod)$term=="VARTRUE"]
    std <- broom::tidy(lm.mod)$std.error[broom::tidy(lm.mod)$term=="VARTRUE"]
    test <- wilcox.test(subset(bg, VAR)$age, subset(bg, !VAR)$age)#, alternative = "less")
    return( data.frame(N_all = N_all, N_var = N_var, effect = effect, std = std, pval=test$p.value))
  }else{
    return( data.frame(N_all = N_all, N_var = 0, effect = 0, std = 0, pval = NA))
  }
}

plot_age <- function(df){
  library(ggthemes)
  library(RColorBrewer)
  bg <- subset(all_tcga, !duplicated(Patient))[c("disease", "Patient", "age")]
  Patients <- df$Patient
  # drop diseases with no variant
  bg <- droplevels(subset(bg, disease %in% droplevels(subset(bg, Patient %in% Patients))$disease))
  p <- ggplot(aes(reorder(disease, age, median), age), data = bg) + geom_boxplot(outlier.shape = NA)
  p <- p + geom_jitter(aes(color=disease), data = subset(bg, Patient %in% Patients))
  p <- p + theme_few()
  p <- p + scale_color_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(length(unique(bg$disease))), guide="none")
  p <- p + xlab("study") + ylab("age")
  p 
}

#var_set <- subset(var_set, (QD >= 1.5&AC<10) | (AC>=10 &QD>=2) )

#load("Results/TCGA_all_mutations.RData")
#tcga_mut
