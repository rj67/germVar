
load("/Volumes/Orchestra/group/TCGA_somatic/July_29/TCGA_all_somatic.merged.anno.RData")

call_set <- call_set %>% subset(., !is.na(Gene))  %>%
  fixSnpEffGene %>%   # rescue unrecognized Gene names
  subset(., Gene %in% list_goi$Gene  & Impact %in% c("MODERATE", "HIGH")) %>%  # remove LOW Impact, synonymous, intron, etc
  labelVarUid

# figure out unique variant set
call_set$mut_uid <- apply(call_set[c("Patient", "var_uid")],1, function(x) paste0(x, collapse="-"))
# remove duplicates
call_set <- call_set[!duplicated(call_set$mut_uid),]

#kinome <- scan("./input/up2ec_pkinfam_merged.txt", what="character")
#write.table(subset(var_set, Gene %in% kinome)[c("Gene", "CHROM", "POS", "REF", "ALT", "X1kG_AF", "AAChange", "CodonChange", "EFF", "ExonRank", "Impact", "PfamDom","Transcript","VARTYPE", "count")],
#      file = "output/tcga_somatic_mutation_kinome.tsv", sep="\t", quote=F, row.names=F)

# all nonsynonymous mutations in gene list of interest
muts_nonsyn <- subset(call_set, Gene %in% list_goi$Gene & EFF == "NON_SYNONYMOUS_CODING")
muts_nonsyn$site <- apply(muts_nonsyn[c("Gene", "AA_pos")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
site_tally <- dplyr::summarise(group_by(muts_nonsyn, site), tcga_scount = length(mut_uid))
# variant representation
tcga_nonsyn <- muts_nonsyn[!duplicated(muts_nonsyn$var_uid), ]
tcga_nonsyn$Patient <- NULL
tcga_nonsyn$matched_norm_sample_barcode<- NULL
tcga_nonsyn$tumor_sample_barcode<- NULL
tcga_nonsyn$mut_uid <- NULL
tcga_nonsyn <- merge(tcga_nonsyn, dplyr::summarise(group_by(muts_nonsyn, var_uid), tcga_vcount = length(mut_uid)))
tcga_nonsyn <- merge(tcga_nonsyn, site_tally, by="site") 
tcga_nonsyn$aa_uid <- apply(tcga_nonsyn[c("Gene", "AAChange")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))


recur_tcga <- droplevels(subset(tcga_nonsyn, tcga_scount >=3))
#recur_soma <- subset(recur_soma, !(var_uid %in% recur_cosm$var_uid))
#recur_soma <- subset(recur_soma, is.na(X1kG_AF)| X1kG_AF<0.01)
#recur_soma <- recur_soma[!duplicated(recur_soma$var_uid), ] 
save(recur_tcga, file="Results/TCGA_recurrent_variants.RData")

