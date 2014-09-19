require(IRanges)
#########################################################
# download the vcf
#########################################################
load("/Volumes/Orchestra/seq/ClinVar/clinvar_20140902.goi.RData")
# only take coding change
clinvar_vcf <- subset(call_set, Gene %in% list_goi$Gene & Impact %in% c("MODERATE", "HIGH"))
# remove inframe INDEL
clinvar_vcf <- subset(clinvar_vcf, !EFF %in% c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION"))
# only take CLNSIG 4 or 5
clinvar_vcf$CLNSIG <- sapply(clinvar_vcf$CLNSIG, function(x) unique(strsplit(x, split="|", fixed=T)[[1]]))
clinvar_vcf$patho <- sapply(clinvar_vcf$CLNSIG, function(x) length(intersect(c("4", "5"), x))>0)
clinvar_vcf <- subset(clinvar_vcf, patho)
# remove duplicates
clinvar_vcf <- createVarUid(clinvar_vcf)
clinvar_vcf <- arrange(clinvar_vcf, var_uid, EFF)
clinvar_vcf <- clinvar_vcf[!duplicated(clinvar_vcf$var_uid), ]


#########################################################
# txt
#########################################################

clinvar_txt <- read.delim(file="../dataDump/ClinVar/variant_summary.20140902.txt")

#setdiff(clinvar_txt$GeneSymbol, table_HGNC$Approved.Symbol)
#setdiff(clinvar_txt$GeneID, table_HGNC$Entrez.Gene)

# intersect with candidate gene list
clinvar_txt <- merge(clinvar_txt, list_goi[c("Gene", "Entrez.Gene")], by.x="GeneID", by.y="Entrez.Gene")

# only take SNP
#clinvar_txt <- droplevels(subset(clinvar_txt, Type=="single nucleotide variant" & GeneSymbol %in% list_goi$Gene))
clinvar_txt <- droplevels(subset(clinvar_txt, Type=="single nucleotide variant" ))
# only take variants with pathogenic annotation or risk factor
clinvar_txt <- clinvar_txt[grepl("Pathogenic", clinvar_txt$ClinicalSignificance, ignore.case=T) | grepl("risk factor", clinvar_txt$ClinicalSignificance), ]
clinvar_txt <- subset(clinvar_txt, Assembly =="GRCh37" )
# remove the splicesite variants which dont have () in their name
clinvar_txt <- clinvar_txt[grepl("(", clinvar_txt$Name, fixed=T ),]
# figure out the CDS change and AA change
clinvar_txt <- cbind(clinvar_txt, t(sapply(clinvar_txt$Name, function(x) { splits1 <- strsplit(x, split=" (", fixed=T)[[1]]; 
                                                               splits11 <- strsplit(splits1[1], split=":")[[1]];
                                                               CDS_Change <- splits11[2];
                                                               REF <- strsplit(substr(CDS_Change, nchar(CDS_Change)-2, nchar(CDS_Change)), split=">")[[1]][1];
                                                               ALT <- strsplit(substr(CDS_Change, nchar(CDS_Change)-2, nchar(CDS_Change)), split=">")[[1]][2];
                                                               AA_Change <- strsplit(splits1[2], split=")")[[1]][1] ;
                                                               return(c("REF"=REF, "ALT"=ALT, "CDS_Change"=CDS_Change, "AA_Change"=AA_Change))} )))
clinvar_txt <- subset(clinvar_txt, ! (is.na(ALT)|is.na(REF)))
# remove the truncating variants, with Ter in AA change
clinvar_txt <- clinvar_txt[!grepl("Ter", clinvar_txt$AA_Change), ]
# uid for variant
clinvar_txt$var_uid <- apply(clinvar_txt[c("Gene", "Start", "REF", "ALT")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
# remove duplicates
clinvar_txt <- arrange(clinvar_txt, var_uid, -nchar(clinvar_txt$OtherIDs))
clinvar_txt <- clinvar_txt[!duplicated(clinvar_txt$var_uid),]

clinvar_txt$CHROM <- clinvar_txt$Chromosome
clinvar_txt$POS <- clinvar_txt$Start
writeBED(clinvar_txt, "./output/clinvar_variants.bed")

