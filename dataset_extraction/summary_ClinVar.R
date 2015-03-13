require(IRanges)
#########################################################
# download the vcf
# #########################################################
# load("/Volumes/Orchestra/seq/ClinVar/clinvar_20140902.goi.RData")
# clinvar_vcf <- rescueSnpEffGenes(call_set)
# # only take coding change
# clinvar_vcf <- subset(call_set, Gene %in% list_goi$Gene & Impact %in% c("MODERATE", "HIGH"))
# # remove inframe INDEL
# clinvar_vcf <- subset(clinvar_vcf, !EFF %in% c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION"))
# # only take CLNSIG 4 or 5
# clinvar_vcf$CLNSIG <- sapply(clinvar_vcf$CLNSIG, function(x) unique(strsplit(x, split="|", fixed=T)[[1]]))
# clinvar_vcf$patho <- sapply(clinvar_vcf$CLNSIG, function(x) length(intersect(c("4", "5"), x))>0)
# clinvar_vcf <- subset(clinvar_vcf, patho)
# # remove duplicates
# clinvar_vcf <- createVarUid(clinvar_vcf)
# clinvar_vcf <- arrange(clinvar_vcf, var_uid, EFF)
# clinvar_vcf <- clinvar_vcf[!duplicated(clinvar_vcf$var_uid), ]


#########################################################
# txt
#########################################################

clinvar_txt <- read.delim(file="../dataDump/ClinVar/variant_summary.20150109.txt")

fixClinGene <- function(vcf){
  Name<-vcf$Name[vcf$GeneSymbol=="-"]
  refGene <- sapply(Name, function(x) {ref<-strsplit(x, split=":", fixed=T)[[1]][1]; 
                                       front<- gregexpr("(", ref, fixed=T)[[1]]; 
                                       end <- gregexpr(")", ref, fixed=T)[[1]]; 
                                       return(ifelse(front>1 & end >1, substr(ref, front+1, end-1), "-")) })
  vcf$GeneSymbol[vcf$GeneSymbol=="-"] <- refGene
  return(vcf)
}

### extract Gene name from Clinvar record when GeneSymbol field is empty
clinvar_txt <- clinvar_txt %>% subset(., Assembly =="GRCh37" ) %>% fixClinGene %>%
               subset(., GeneSymbol %in% list_goi$Gene) %>% # intersect with candidate gene list
               plyr::rename(., replace=c("GeneSymbol" = "Gene"))

### tally the Gene names that have pathogenic entries in not single nucleotide variant. 
clinvar_lof <- clinvar_txt %>% subset(., !Type %in% c("single nucleotide variant",  "undetermined variant", "NT expansion", "copy number gain", "duplication")) %>%
               subset(., grepl("Pathogenic", ClinicalSignificance, ignore.case=T)|grepl("risk factor", ClinicalSignificance)) %>%
               group_by(., Gene) %>% dplyr::summarise(., N_indel_entry = length(Name))

clinvar_txt <- clinvar_txt %>%  subset(., Type=="single nucleotide variant" ) %>% # only take SNP
               subset(., grepl("Pathogenic", ClinicalSignificance, ignore.case=T)|grepl("risk factor", ClinicalSignificance)) %>% # only take variants with pathogenic annotation or risk factor
               subset(., grepl(">", Name, fixed=T)) %>% # make sure variant name has > in it
               droplevels

### figure out the codon change from the name
clinvar_txt$anchor <- sapply(clinvar_txt$Name,  function(x) gregexpr(pattern ='>', x)[[1]][1])
clinvar_txt <- clinvar_txt %>% mutate(., REF = substr(Name, anchor-1, anchor-1),  ALT = substr(Name, anchor+1, anchor+1)) %>%
               plyr::rename(., replace=c("Chromosome"="CHROM", "Start"="POS")) 

### figure out the nucleotide change in genome space
clinvar_txt <- plyr::join(clinvar_txt, list_gff_long[c("Gene","strand")], by="Gene")
clinvar_txt$REF <- mapply(flipCodon, clinvar_txt$REF, clinvar_txt$strand)
clinvar_txt$ALT <- mapply(flipCodon, clinvar_txt$ALT, clinvar_txt$strand)
clinvar_txt <-  labelVarUid(clinvar_txt)
clinvar_txt <- arrange(clinvar_txt, var_uid, -nchar(clinvar_txt$OtherIDs))
clinvar_txt <- clinvar_txt[!duplicated(clinvar_txt$var_uid),]

### write out the variant, run SNPEFF on it
writeVCF(clinvar_txt, "Output/clinvar_txt.20150109.vcf")

###
clinvar_vcf <- readVcf(file="./DataSet_Results/clinvar_txt.20150109.snpEff.vcf.gz", "hg19")

#load("/Volumes/Orchestra/seq/ClinVar/clinvar_txt.snpEff.RData")
#call_set$AA_pos <- sapply(call_set$AAChange, function(x) {y=gregexpr("[0-9]+", x)[[1]]; ifelse(attributes(y)$useBytes, substr(x, y, y+ attributes(y)$match.length-1), NA) })
clinvar_vcf <-  clinvar_vcf %>% labelUid   %>%   
  subset(., !is.na(Coding) & Coding == "CODING" ) %>%
  subset(., !is.na(BioType) & BioType == "protein_coding" ) %>%
  fixSnpEffGene %>%
  subset(., Gene %in% list_goi$Gene) 

# tally the lof variants
clinvar_lof <- merge(clinvar_lof, info(clinvar_vcf) %>% as.data.frame %>% subset(., Impact =="HIGH") %>% group_by(., Gene) %>% dplyr::summarise(., N_snp_entry = length(unique(var_uid))), all =T)
clinvar_lof[c("N_indel_entry", "N_snp_entry")] <- sapply(clinvar_lof[c("N_indel_entry", "N_snp_entry")], replaZero)

clinvar_vcf <- clinvar_vcf %>% subset(., EFF=="missense_variant" ) %>%
  labelSnpEffUid %>%
  labelAAUid %>% 
  labelSite 
  
clinvar_txt <- merge(clinvar_txt, as.data.frame(info(clinvar_vcf)[c("var_uid", "AAChange", "AAChange.p", "AAChange.c", "AA_pos", "AALength", "aa_uid", "Transcript", "site", "EFF", "snpeff_uid")]), by="var_uid")



# uid for variant
# remove duplicates

save( clinvar_txt, clinvar_lof, file="Results/ClinVar.RData")
rm(clinvar_vcf, fixClinGene)

