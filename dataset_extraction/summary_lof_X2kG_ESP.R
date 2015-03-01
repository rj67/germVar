################################################################
###  X2kG
################################################################
# read vcf using VariantAnnotation interface
lof_X2kG <- expand(readVcf(file="./DataSet_Results/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.goi.lof.vcf.gz", "hg19"))
for (coln in c("CIEND", "CIPOS", "CS", "END",	"IMPRECISE", "MC", "MEINFO", "MEND", "MLEN", "MSTART", "SVLEN", "SVTYPE", "TSD")){
  info(lof_X2kG)[[coln]] <- NULL 
}

lof_X2kG <- lof_X2kG %>% labelUid  %>%   
  subset(., !is.na(Transcript) ) %>%
  subset(., Coding == "CODING" & BioType == "protein_coding" ) %>%
  subset(., is.na(BIOTYPE) | BIOTYPE == "protein_coding" ) %>%
  fixSnpEffGene %>%
  subset(., Gene %in% list_goi$Gene) %>%
  subset(., Impact=="HIGH" ) %>%
  subset(., is.na(SYMBOL) | Gene==SYMBOL) %>% 
  labelVarUid %>% 
  labelSnpEffUid  %>% 
  labelVEPUid(., "trunc") 

### pick longest unique CCDS transcripts for SnpEFF annotation
lof_X2kG <- pickSnpEffTranscript(lof_X2kG, CCDS_15)

# fix LOF_P_Transcripts
info(lof_X2kG)$LOF_P_Transcripts <- replaZero(info(lof_X2kG)$LOF_P_Transcripts)
info(lof_X2kG)$NMD_P_Transcripts <- replaZero(info(lof_X2kG)$NMD_P_Transcripts)

### Simplify VEP Consequence field, Remove VEP annotations with no Consequence, excluding variants that don't have any Consequence annotation
lof_X2kG <- lof_X2kG %>% simplifyVEPConseq

### Remove VEP annotation fields that aren't in CCDS, excluding variants that don't have any CCDS annotation
lof_X2kG <- lof_X2kG %>% filterVEPTranscript(., CCDS_15[["nucleotide_ID"]])

### Rank vep_uid according to LoF, only take one vep_uid per var_uid
lof_X2kG <- lof_X2kG %>% pickVEPTranscript

# check vairant effects, remove all variants that contain stop_lost since it's near the end
table(info(lof_X2kG)$EFF)
lof_X2kG <- lof_X2kG %>% subset(., sapply(EFF, function(x) !grepl("stop_lost", x)))

#remove genes that are not in trunc_VAR 
lof_X2kG <- subset(lof_X2kG, Gene %in% info(trunc_VARu)$Gene)


# consolidate SnpEff LoF prediction
info(header(lof_X2kG))["LoF_snpeff",] <- list("1" ,"Flag", "LoF_snpeff")
info(lof_X2kG)$LoF_snpeff <- with(as.data.frame(info(lof_X2kG)[c("NMD_P_Transcripts", "LOF_N")]), NMD_P_Transcripts>0 | LOF_N>0) 

with(as.data.frame(info(lof_X2kG)[c("LoF_snpeff", "LoF")]), table(LoF_snpeff, LoF))
info(header(lof_X2kG))["tier1",] <- list("1" ,"Flag", "tier1")
info(header(lof_X2kG))["tier2",] <- list("1" ,"Flag", "tier2")
info(lof_X2kG)$tier1 <- with(as.data.frame(info(lof_X2kG)[c("LoF_snpeff", "LoF")]), LoF_snpeff & LoF!= "LC")
info(lof_X2kG)$tier2 <- with(as.data.frame(info(lof_X2kG)[c("LoF_snpeff", "LoF")]), (LoF_snpeff & LoF== "LC" )| (!LoF_snpeff & LoF=="HC" ))
#remove LCs


################################################################
###  ESP
################################################################
# read vcf using VariantAnnotation interface
lof_ESP <- expand(readVcf(file="./DataSet_Results/ESP6500SI-V2-SSA137.updatedProteinHgvs.snps_indels.goi.lof.vcf.gz", "hg19"))
for (coln in c("MAF", "GTS", "EA_GTC", "AA_GTC",  "GTC", "FG", "HGVS_CDNA_VAR", "HGVS_PROTEIN_VAR", "CDS_SIZES", 
              "PH", "CP", "CG", "GL", "GS", "CA", "EXOME_CHIP", "GWAS_PUBMED", "EA_AGE", "AA_AGE", "EA_AC", "AA_AC", "TAC")){
  info(lof_ESP)[[coln]] <- NULL 
}

lof_ESP <- lof_ESP %>% labelUid  %>%   
  subset(., !is.na(Transcript) ) %>%
  subset(., Coding == "CODING" & BioType == "protein_coding" ) %>%
  subset(., is.na(BIOTYPE) | BIOTYPE == "protein_coding" ) %>%
  fixSnpEffGene %>%
  subset(., Gene %in% list_goi$Gene) %>%
  subset(., Impact=="HIGH" ) %>%
  subset(., is.na(SYMBOL) | Gene==SYMBOL) %>% 
  labelVarUid %>% 
  labelSnpEffUid  %>% 
  labelVEPUid(., "trunc") 

# fix LOF_P_Transcripts
info(lof_ESP)$LOF_P_Transcripts <- replaZero(info(lof_ESP)$LOF_P_Transcripts)
info(lof_ESP)$NMD_P_Transcripts <- replaZero(info(lof_ESP)$NMD_P_Transcripts)


lof_ESP <- lof_ESP %>% pickSnpEffTranscript(., CCDS_15)  %>% ### pick longest unique CCDS transcripts for SnpEFF annotation
           simplifyVEPConseq %>%                ### Simplify VEP Consequence field, Remove VEP annotations with no Consequence, excluding variants that don't have any Consequence annotation
           filterVEPTranscript(., CCDS_15[["nucleotide_ID"]]) %>%  ### Remove VEP annotation fields that aren't in CCDS, excluding variants that don't have any CCDS annotation
           pickVEPTranscript %>%                ### Rank vep_uid according to LoF, only take one vep_uid per var_uid
           subset(., sapply(EFF, function(x) !grepl("stop_lost", x))) ## remove all variants that contain stop_lost since it's near the end


#remove genes that are not in trunc_VAR 
lof_ESP <- subset(lof_ESP, Gene %in% info(trunc_VARu)$Gene)

# consolidate SnpEff LoF prediction
info(header(lof_ESP))["LoF_snpeff",] <- list("1" ,"Flag", "LoF_snpeff")
info(lof_ESP)$LoF_snpeff <- with(as.data.frame(info(lof_ESP)[c("NMD_P_Transcripts", "LOF_N")]), NMD_P_Transcripts>0 | LOF_N>0) 

with(as.data.frame(info(lof_ESP)[c("LoF_snpeff", "LoF")]), table(LoF_snpeff, LoF))
info(header(lof_ESP))["tier1",] <- list("1" ,"Flag", "tier1")
info(header(lof_ESP))["tier2",] <- list("1" ,"Flag", "tier2")
info(lof_ESP)$tier1 <- with(as.data.frame(info(lof_ESP)[c("LoF_snpeff", "LoF")]), LoF_snpeff & LoF!= "LC")
info(lof_ESP)$tier2 <- with(as.data.frame(info(lof_ESP)[c("LoF_snpeff", "LoF")]), (LoF_snpeff & LoF== "LC" )| (!LoF_snpeff & LoF=="HC" ))

novel_ids <- setdiff(subset(info(lof_ESP), tier1)$var_uid, subset(info(lof_X2kG), tier1)$var_uid)
View(subset(info(lof_ESP), var_uid %in% novel_ids & ESP_AC>= 10))

to_plot <- as.data.frame(info(lof_ESP)[c("var_uid", "ESP_AF", "VARTYPE")])
to_plot$NOVEL <- !to_plot$var_uid %in% info(lof_X2kG)$var_uid
to_plot$ESP_AF <- do.call(c, to_plot$ESP_AF)

p <- ggplot(aes(VARTYPE, ESP_AF), data=to_plot) + geom_jitter(alpha=0.5, position=position_jitter(height=0.1)) + theme_few()
p <- p + scale_color_manual(values=c("grey20", "red")) + scale_y_log10() + facet_grid(.~NOVEL)
p
