setwd("/Users/snafu/Documents/Project/germVar")
table_HGNC_file <- "../dataDump/HGNC/HGNC_ID_table_Nov11.txt"

table_cgc500_file <- "./Gene_List/cancer_gene_census_v71.csv"
table_tsg716_file <- "../dataDump/TSGene/Human_716_TSGs.txt"
new_tsg145_file <- "./Gene_List/New_145_TSGs.txt"
hc_tsg206_file <- "./Gene_List/High_confidence_206_TSGs.txt"
table_smg260_file <- "./Gene_List/SMG_genes_260.csv"
table_sfe_file <- "./Gene_List/pancan16_mini_maf_annotated.txt"
fam_files_path <- "./Gene_List/"
rep_files_path <- "./Gene_List/"
table_tag700_file <- "./Gene_List/Tumor_Associated_Genes.csv"
ref_TUSON_file <- "../dataDump/TUSON/mmc2.xls"
table_TUSON_A_file <- "../dataDump/TUSON/mmc3.A.csv"
table_TUSON_B_file <- "../dataDump/TUSON/mmc3.B.csv"
table_TUSON_D_file <- "../dataDump/TUSON/mmc3.D.csv"
table_cgl_file <- "./Gene_List/Vogelstein_Cancer_Genome_Landscapes_1235122TablesS1-4.xls"

# -------------------------------------------------------------------------------------
#  HGNC table, useful for human gene name, gene id, uniprot mapping
# -------------------------------------------------------------------------------------
table_HGNC <- read.delim(table_HGNC_file)
# shortern the column names
colnames(table_HGNC)[grep("supplied", colnames(table_HGNC), fixed=T)] <- sapply(colnames(table_HGNC)[grep("supplied", colnames(table_HGNC), fixed=T)],
                                                                          function(x) paste(strsplit(x, ".", fixed = T)[[1]][1:2], collapse = "."))
rm(table_HGNC_file)

# -------------------------------------------------------------------------------------
#  COSMIC 500 gene
# -------------------------------------------------------------------------------------
table_cgc500 <- read.csv(table_cgc500_file, strip.white=T)
# BCL5 a locus
# C12orf9 is retracted in HGNC 
# C15orf21 is pseudogene
# RUNDC2A have  outdated GeneID
# DUX4 reference, but probably meant DUX4L1 which is a pseudogene
# STL/TCL6/MALAT1 are RNA gene, not in the HGNC table
subset(table_cgc500, Entrez.GeneId %in% setdiff(table_cgc500$Entrez.GeneId, table_HGNC$Entrez.Gene))
# fix RUNDC2A
table_cgc500$Entrez.GeneId[table_cgc500$Gene.Symbol=="RUNDC2A"] <- 92017
table_cgc500$Mutation.Types <- gsub( " ", "", table_cgc500$Mutation.Types )
#merge the HGNC gene symbols with the cosmic table
table_cgc500 <- merge(table_cgc500, table_HGNC[c("Approved.Symbol", "Entrez.Gene", "UniProt.ID")], by.x = "Entrez.GeneId", by.y = "Entrez.Gene")
rm(table_cgc500_file)


# -------------------------------------------------------------------------------------
#  716 Tumor Suppressor gene
# -------------------------------------------------------------------------------------
table_tsg716 <- read.delim(table_tsg716_file, strip.white = T, header = T)
# only keep the first 8 columns
table_tsg716 <- table_tsg716[,1:8] 
# remove miscRNA, only include protein coding gene
table_tsg716 <- subset(table_tsg716, Gene_type=="protein-coding")
# make sure all genes are in HGNC
setdiff(table_tsg716$GeneID, table_HGNC$Entrez.Gene)
#GSTT1 not in HGNC, for some reason, something about alternate reference loci HGNC:4641
# merge the HGNC gene symbols with the tsgene table
table_tsg716 <- merge(table_tsg716, table_HGNC[c("Approved.Symbol", "Entrez.Gene", "UniProt.ID")], by.x = "GeneID", by.y = "Entrez.Gene")
# 5 genes with different names in 716 vs HGNC
table_tsg716[c("Gene_symbol", "Approved.Symbol")][!apply(table_tsg716[c("Gene_symbol", "Approved.Symbol")],1, function(x) {x[1]==x[2]}),]

# additional 154
new_tsg145 <- read.table(new_tsg145_file, strip.white = T, header =T)
subset(new_tsg145, GeneID %in% setdiff(new_tsg145$GeneID, table_HGNC$Entrez.Gene))
new_tsg145 <- merge(new_tsg145[c("GeneID")], table_HGNC[c("Approved.Symbol", "Entrez.Gene", "UniProt.ID")], by.x = "GeneID", by.y = "Entrez.Gene")
#remove the ones without UniProt
new_tsg145 <- subset(new_tsg145, UniProt.ID!="")
table_tsg716<-plyr::rbind.fill(table_tsg716, new_tsg145)

# high confidence 206
hc_tsg206 <- scan(hc_tsg206_file, what = "character")
table_tsg716$High_conf <- table_tsg716$GeneID %in% hc_tsg206
rm(hc_tsg206, hc_tsg206_file, new_tsg145, new_tsg145_file, table_tsg716_file)



# -------------------------------------------------------------------------------------
#  Significantly mutated genes additional 260 genes
# -------------------------------------------------------------------------------------
table_smg260 <- read.csv(table_smg260_file)
# a few gene names are outdated
setdiff(table_smg260$gene, table_HGNC$Approved.Symbol)
table_smg260$gene[table_smg260$gene=='SFRS2'] <- 'SRSF2'
table_smg260$gene[table_smg260$gene=='MLL'] <- 'KMT2A'
table_smg260$gene[table_smg260$gene=='MLL2'] <- 'KMT2D'
table_smg260$gene[table_smg260$gene=='MLL3'] <- 'KMT2C'
table_smg260$gene[table_smg260$gene=='MLL4'] <- 'KMT2B'
rm(table_smg260_file)



# # -------------------------------------------------------------------------------------
# #  Giovanni Pan-Can paper significant function event (SFE)
# # -------------------------------------------------------------------------------------

#table_sfe <- read.delim(file=table_sfe_file)
#colnames(table_sfe)[1] <- "Gene"
#table_sfe$Gene[table_sfe$Gene=='MLL2'] <- 'KMT2D'
#table_sfe$Gene[table_sfe$Gene=='SFRS2'] <- 'SRSF2'
#setdiff(table_sfe$Gene, c(table_cgc500$Approved.Symbol, table_tsg716$Approved.Symbol, table_smg260$gene))
#table_sfe <- subset(table_sfe, hotspot_count >=3 & Variant_Classification !="Silent")
#table_sfe$var_uid <- apply(table_sfe[c("Gene", "ONCOTATOR_PROTEIN_CHANGE")], 1, function(x) paste0(x, collapse="-"))
#table_sfe <- table_sfe[!duplicated(table_sfe$var_uid),]
#rm(table_sfe_file)

# -------------------------------------------------------------------------------------
# Various Gene familis from HGNC
# -------------------------------------------------------------------------------------
fam_names <- gsub("_family.txt", "", list.files(path = fam_files_path, pattern="_family.txt", full.names=F), fixed=T)
fam_genes <- do.call(rbind, lapply(fam_names, function(fam) {
  fam_table <- read.delim(file=paste(fam_files_path, fam, "_family.txt", sep=""), head=T, strip.white=T, stringsAsFactors=F);
  returns <- data.frame("Gene" = fam_table[["Approved.Symbol"]]);
  returns$fam_id <- fam;
  print(dim(returns));
  return(returns)
}))
rm(fam_names, fam_files_path)


# -------------------------------------------------------------------------------------
# Various DNA repair genes 
# -------------------------------------------------------------------------------------
rep_names <- gsub("_REP.txt", "", list.files(path = rep_files_path, pattern="_REP.txt", full.names=F), fixed=T)
rep_genes <- do.call(rbind, lapply(rep_names, function(path) {
  rep_table <- read.delim(file=paste(rep_files_path, path, "_REP.txt", sep=""), skip = 1, head=F, strip.white=T, stringsAsFactors=F);
  returns <- data.frame("Gene" = rep_table[["V1"]]);
  returns$rep_id <- path;
  return(returns)
}))
rm(rep_names, rep_files_path)


# -------------------------------------------------------------------------------------
# manually add groups of genes, complexes
# -------------------------------------------------------------------------------------
# Pro-apoptotic factors
BH3 <- c( "BAX", "BAK1", "BOK", "BID", "BAD", "BIK", "BBC3", "PMAIP1", "BMF", "HRK")
# Nu-A4 complex
NUA4 <- c( "TRRAP",  "EP400", "BRD8", "EPC1", "KAT5", "DMAP1", "ING3", "VPS72", "RUVBL1", "RUVBL2", "MORF4L1", "YEATS4", "MRGBP", "MEAF6")
# PRC2
PRC2 <- c( "SUZ12", "EED", "EZH1", "EZH2", "RBBP4", "RBBP7", "AEBP2", "JARID2", "PHF1")
# NuRD
NURD <- c( "CHD3", "CHD4", "HDAC1", "HDAC2", "MBD2", "MBD3", "MTA1", "MTA2", "MTA3", "GATAD2A", "GATAD2B")
# SWISNF
SWISNF <- c( "ARID1A", "ARID1B", "ARID2", "SMARCA2", "SMARCA4", "PBRM1", "SMARCC2", "SMARCC1", "SMARCD1", "SMARCD2", "SMARCD3", "SMARCE1", "ACTL6A", "ACTL6B", "SMARCB1", "BRD7")
# Pre-mRNA processing
mRNA <- c("PRPF40B", "SF1", "ZRSR2", "U2AF1", "U2AF2", "SRSF2", "SF3B1", "SF3A1")
# miRNA processing
miRNA <- c("DICER1", "AGO1", "AGO2", "TARBP2", "DROSHA") 
# manual literature curation
LIT <- c("BCORL1", "PRDM5", "DNMT1","DNMT3B", "HAND2", "ADAMTS15", "FGFR4",
              "HIST1H1B", "HIST1H1C", "ELP4", "BAG6", "CHD1", "CHD2","PIK3R2", "TRIP12", 
               "BUB1", "BUB3","LZTR1", "SYNE1", 
              "PREX2", "SOX9", "AMER2", "LDLRAP1", "STMN2", "AGTR2", "ZIM2", "NALCN", "SLC16A4", "MAGEA6", "RPS6KA3", "ROBO2", 
              "HNRNPK", "BRINP3", "KLHL6", "MEF2B", "ACTB", "PCLO", "LRRK2",
              "MTF1", "GALNT1", "TINF2", "WRAP53", "ELANE", "DKC1", "EPCAM" , "FBP1", "JMJD1C", "USP28", "PTPRF", "BRSK1", "AR", "CSF1R")

list_add <- list(BH3 = BH3, NUA4 = NUA4, PRC2 = PRC2, NURD = NURD, SWISNF = SWISNF, mRNA = mRNA, miRNA = miRNA, LIT=LIT)

grp_genes <-  do.call(rbind, lapply(c("BH3", "NUA4", "PRC2", "NURD", "SWISNF", "mRNA", "miRNA", "LIT"), function(x) data.frame(Gene=list_add[[x]], grp_id = x)))
rm(BH3, NUA4, PRC2, NURD, SWISNF, mRNA, miRNA, LIT, list_add)


# -------------------------------------------------------------------------------------

# genes in S500 not S200 but should be in S200
#list_s200_ext <-c("SUZ12", "CBFB", "CCND1", "CCND3", "ETV6", "ITK", "NTRK1", "NTRK3", "RUNX1", "TCF7L2", "SEPT9", "FGFR1", "PER1", "RAD51B" )
#the reverse
#list_s200_exc <-c("KIAA1549")

# -------------------------------------------------------------------------------------
# list of genes with activating mutation or overexpression that leads to tumorgenesis
#list_active <- c("EGFR", "FGFR1", "FGFR2", "FGFR3", "GLDC", "GNAS", "KRAS", "HRAS", "NRAS", "IDH1", "IDH2", "ALK", "ERBB2", "ERBB3", "BRAF", "PIK3CA", "PIK3R1", "PIK3R2", "CTNNB1", "AKT1", "MAP2K1", "MYC",
#                 "FLT3", "IL7R",  "PTPN11", "KDR", "FGFR4","JAK1", "JAK2", "JAK3", "MYCN", "CRCLF2", "ABL1", "NT5C2", "PDGFRA", "MAP2K2", "MPL", "MET", "SETBP1", "TRRAP", "MYD88", "KIT", "ITK",
#                 "ALK", "NOTCH1", "NOTCH2", "CARD11", "ECT2L", "FOXL2", "CBL", "KLF4", "SMO", "STAT3", "STAT5B","MTOR", "BCL6", "XPO1", "TSHR", "U2AF1", "SF3B1", "SRSF2", "GNAQ", "GNA11", "NTRK3" ,"NTRK1", 
#                 "RAC1", "SEPT9", "TERT", "TBL1XR1", "CDK4", "CDK12", "CCND1", "CCND3", "CD79B", "CALR", "CRLF2", "CACNA1D", "ATP1A1", "DNMT3A", "KCNJ5", "NFE2L2", "NPM1", "PHOX2B", "PPP2R1A",
#                 "H3F3A", "H3F3B", "GATA1", "CBLB", "RET", "TRAF7", "FOXA1", "DNM2", "RPL5")

# list of genes suspected TS
#list_ts <- c("ATRX", "CIC", "BCOR", "BIRC3","BMPR1A", "CBFB", "CNOT3", "DAXX", "DICER1", "EP300", "EZH2", "FBXO11", "FUBP1", "KDM5C", "KMT2A", "KMT2C", "KMT2D", "POT1", "RNF43", "RUNX1", "SDHC", "SETD2", "SH2B3", "STAG2", "SUZ12", "UBR5", "WAS", "WRN", "ZRSR2",
#             "FAM46C", "FAS", "GATA2", "GATA3", "TNFRSF14", "HNF1A", "RECQL4", "RAD21", "KDM6A", "ARHGAP26")

# -------------------------------------------------------------------------------------
list_goi <-c(table_tsg716$Approved.Symbol, table_cgc500$Approved.Symbol, table_smg260$gene,  rep_genes$Gene, fam_genes$Gene, grp_genes$Gene)
list_goi <- unique(list_goi)
list_goi <- data.frame("Gene" =list_goi)
list_goi$cgc500 <- list_goi$Gene %in% as.character(table_cgc500$Approved.Symbol)
list_goi$cgc200 <- with(list_goi, cgc500 & !(Gene %in% subset(table_cgc500, Mutation.Types %in% c("A", "T", "A,T", "T,A"))$Approved.Symbol))
#list_goi$Dom <- list_goi$Gene %in% as.character(with(table_cgc500, Approved.Symbol[grep('Dom', Cancer.Molecular.Genetics)]))
list_goi$tsg716 <- list_goi$Gene %in% as.character(table_tsg716$Approved.Symbol)
list_goi$tsg206 <- list_goi$Gene %in% as.character(subset(table_tsg716, High_conf)$Approved.Symbol)
list_goi$smg260 <- list_goi$Gene %in% table_smg260$gene
list_goi$smg33 <- list_goi$Gene %in% subset(table_smg260,Discussed.in.text..33.==1)$gene
#list_goi$sfe160 <- list_goi$Gene %in% table_sfe$Gene
list_goi$rep <- list_goi$Gene %in% rep_genes$Gene
list_goi <- merge(list_goi, rep_genes, by="Gene", all.x=T)
list_goi$fam <- list_goi$Gene %in% fam_genes$Gene
list_goi <- merge(list_goi, fam_genes, by="Gene", all.x=T)
list_goi$grp <- list_goi$Gene %in% grp_genes$Gene
list_goi <- merge(list_goi, grp_genes, by="Gene", all.x=T)


#get the full name
list_goi <- merge(list_goi, table_HGNC[c("Approved.Symbol", "Approved.Name", "Entrez.Gene", "UniProt.ID", "Ensembl.ID", "CCDS.IDs", "Locus.Type", "Locus.Group")], by.x = "Gene", by.y = "Approved.Symbol")

# remove RNA genes and pseudogene
print( subset(list_goi, Locus.Group!= "protein-coding gene" ))
list_goi <- subset(list_goi, Locus.Type == "gene with protein product")
list_goi$Locus.Type <- NULL
list_goi$Locus.Group <- NULL

# remove clusters
list_goi <- subset(list_goi, ! Gene %in% c("IGH", "IGL", "IGK", "TRA", "TRB", "TRD", "HLA-A", "HLA-B"))

rm(fam_genes, grp_genes, rep_genes, table_cgc500, table_smg260, table_tsg716)
#list_goi$s200[list_goi$Gene %in% list_s200_ext] <- TRUE 
#list_goi$s200[list_goi$Gene %in% list_s200_exc] <- FALSE



# -------------------------------------------------------------------------------------
#  Tumor Associated Gene 700 genes
# -------------------------------------------------------------------------------------
table_tag700 <- read.csv(table_tag700_file, stringsAsFactors =F)
# keep only genes with clear Onco/TS annotation
table_tag700 <- subset(table_tag700, !Category %in% c("--", "Other"))
#table_tag700$Symbol[table_tag700$Symbol=="CDC2L1"] <- "CDK11B"
#table_tag700$Symbol[table_tag700$Symbol=="MYCL1"] <- "MYCL"
#table_tag700$Symbol[table_tag700$Symbol=="KIAA1245"] <- ""

table_tag700_in <- subset(table_tag700, Symbol %in% table_HGNC$Approved.Symbol)
table_tag700_out <- subset(table_tag700, !Symbol %in% table_HGNC$Approved.Symbol)

# a few gene names are outdated
tmp<-subset(table_HGNC, Approved.Symbol %in% list_goi$Gene)[c("Approved.Symbol", "Previous.Symbols")]
tmp<-subset(tmp, Previous.Symbols!="")
tmp$Previous.Symbols <- levels(tmp$Previous.Symbols)[tmp$Previous.Symbols]
tmp <- do.call(rbind, apply(tmp, 1, function(x) data.frame( Approved.Symbol=x[1], Alias=strsplit(gsub(" ", "", x[2]), split=",")[[1]])))
# merge back
table_tag700_out <- merge(table_tag700_out, tmp, by.x="Symbol", by.y="Alias")
table_tag700_out$Symbol <- table_tag700_out$Approved.Symbol
table_tag700_out$Approved.Symbol <- NULL
table_tag700 <- rbind(table_tag700_in, table_tag700_out)
colnames(table_tag700) <- gsub("Symbol", "Gene", colnames(table_tag700))

# -------------------------------------------------------------------------------------
#  Cancer Genome Landscape review paper
# -------------------------------------------------------------------------------------
library(XLConnect)
# 
table_cgl <- readWorksheet(loadWorkbook(table_cgl_file), sheet=6, startRow=2)
table_cgl <- subset(table_cgl, !is.na(Classification.))
table_cgl$Gene.Symbol <- gsub("FAM123B", "AMER1", table_cgl$Gene.Symbol)
table_cgl$Gene.Symbol <- gsub("MLL2", "KMT2D", table_cgl$Gene.Symbol)
table_cgl$Gene.Symbol <- gsub("MLL3", "KMT2C", table_cgl$Gene.Symbol)
table_cgl2 <- readWorksheet(loadWorkbook(table_cgl_file), sheet=7, startRow=2)
table_cgl2 <- subset(table_cgl2, !is.na(Classification.))
table_cgl2$Gene.Symbol <- gsub("MYCL1", "MYCL", table_cgl2$Gene.Symbol)
rm(table_cgl_file)
colnames(table_cgl) <- gsub("Gene.Symbol", "Gene", colnames(table_cgl))
colnames(table_cgl2) <- gsub("Gene.Symbol", "Gene", colnames(table_cgl2))

# -------------------------------------------------------------------------------------
#  TUSON TSG and OG
# -------------------------------------------------------------------------------------
library(XLConnect)
# TSG/OG training 
ref_TUSON <- readWorksheet(loadWorkbook(ref_TUSON_file), sheet=1, startRow=3)
#MYCL1 -> MYCL
ref_TUSON$OG <- gsub("MYCL1", "MYCL", ref_TUSON$OG)
# TSG/OG results
#table_TUSON_tsg <- readWorksheet(loadWorkbook(table_TUSON_file, create=F), sheet=1, startRow=3)
table_TUSON_tsg <- read.csv(table_TUSON_A_file, skip=2)
table_TUSON_tsg <- subset(table_TUSON_tsg, TUSON_q_value_TSG< 0.01)
#table_TUSON_og <- readWorksheet(loadWorkbook(table_TUSON_file, create=F), sheet=2, startRow=3)
table_TUSON_og <- read.csv(table_TUSON_B_file, skip=2)
table_TUSON_og <- subset(table_TUSON_og, TUSON_q_value_OG< 0.01)
#table_TUSON_fam <- readWorksheet(loadWorkbook(table_TUSON_file, create=F), sheet=4, startRow=3)
table_TUSON_fam <- read.csv(table_TUSON_D_file, skip=2)

rm(ref_TUSON_file, table_TUSON_A_file, table_TUSON_B_file, table_TUSON_D_file)

# tag700 designation/ Cancer Genome landscape supplimentary / TUSON designation
list_goi$TSG <- list_goi$Gene %in% c( subset(table_tag700, Category=="Tumor suppressor gene")$Gene, 
                                      subset(table_cgl, Classification.=="TSG")$Gene, subset(table_cgl2, Classification.=="TSG")$Gene,
                                      ref_TUSON$TSG, table_TUSON_tsg$Gene, table_TUSON_fam$Gene )
list_goi$OG <- list_goi$Gene %in% c(subset(table_tag700, Category=="Oncogene")$Gene,
                                    subset(table_cgl, Classification.!="TSG")$Gene, subset(table_cgl2, Classification.!="TSG")$Gene,
                                    ref_TUSON$OG, table_TUSON_og$Gene )
# name contains proto-oncogene
#list_goi$OG[grep("oncogene", list_goi$Approved.Name)] <- T
# hereditary cancer syndrome
list_goi$HER <- list_goi$Gene %in% table_TUSON_fam$Gene
# manually add a few genes
list_goi$HER[list_goi$Gene %in% c("RAD51B", "RAD51C", "RAD51D", "MRE11A", "FAM175A")] <- TRUE

#conflicting
#SMO, PAX5, ELF3, NOTCH1, NOTCH2, DNMT3A, KLF4, EZH2, TGFB1

list_goi$Cat <- apply(list_goi[c("TSG", "OG")], 1, function(x) ifelse(x[1], "TSG", ifelse(x[2], "OG", "OTHER")))
list_goi$Cat[list_goi$Gene %in% c("SMO", "PAX5", "ELF3", "NOTCH1", "NOTCH2", "DNMT3A", "KLF4", "EZH2", "TGFB1", "ERBB4", "RHOB")] <- "DUAL"
list_goi$TSG <- NULL
list_goi$OG <- NULL
list_goi$ONCO <- NULL
# get gene membership in the different sources
list_goi$mem <- mapply(function(cgc, smg, tsg, rep){ifelse(cgc, "CGC", ifelse(smg, "SMG", ifelse(tsg, "TSG", ifelse(rep, "REP","OTHER"))))}, list_goi$cgc500, list_goi$smg260, list_goi$tsg716, list_goi$rep)
list_goi$mem <- factor(list_goi$mem, levels=c("CGC", "SMG", "TSG", "REP", "OTHER"))
rm(table_tag700, table_tag700_file, table_tag700_in, table_tag700_out, tmp, table_TUSON_tsg, table_TUSON_og, ref_TUSON, table_TUSON_fam, table_cgl, table_cgl2)

# output
write.table(list_goi, file="Output/candidate_gene_list.tsv", quote=F, row.names=F, sep="\t")
save(table_HGNC, file="Results/table_HGNC.RData")
save(list_goi, file="Results/candidate_gene_list.RData")

#goi_exons <- read.delim(file="input/candidate_gene_hg19_exons.bed", stringsAsFactors=F, header=F)
#goi_hap <- rbind(goi_exons[grep("hap", goi_exons$V1),], goi_exons[grep("chrUn", goi_exons$V1),])
#goi_hap$refseq <- sapply(goi_hap$V4, function(x) strsplit(x, split="_cds_")[[1]][1])
#goi_hap <- goi_hap[!duplicated(goi_hap$refseq),]

if(F){
library(GenomicFeatures)
supportedUCSCtables()
hg19.refgene.tx <- makeTranscriptDbFromUCSC(genome = "hg19", tablename = "refGene")
save(hg19.refgene.tx, file="Results/hg19_refgene_tx.RData")


list_goi_exons<-exons(tmp, list(gene_id = list_goi$Ensembl.ID),columns=c("gene_id","exon_id"))

#write exons
list_goi_exons<-exons(hg19.refgene.tx, list(gene_id = list_goi$Entrez.Gene),columns=c("gene_id","exon_id"))
# show Genes that are in alternative haplotype region
to_write <- Grange2bed(list_goi_exons)
tmp2 <- to_write[grep("hap", to_write$seqnames),]
subset(list_goi, Entrez.Gene %in% unique(unlist(tmp2$gene_id)))
subset(list_goi, Entrez.Gene %in% unique(unlist(subset(tmp, seqnames=="Y")$gene_id)))

#convert Grange object to 
Grange2bed<-function(GrangeObject){
  returns<-droplevels(as.data.frame(GrangeObject))
  returns$seqnames<-sapply(returns$seqnames,function(x) gsub("chr","",x))
  return(returns)
}

# write the exons, remove alternative haplotype, extend exons by 10bp in both direction
to_write <- to_write[-grep("hap", to_write$seqnames),]
to_write$start <- to_write$start - 10
to_write$end <- to_write$end + 10

writeBed <- function(to_write, outfile){
  # to_write is a dataframe starting with chrom, start, end columns.
  # sort the bed and merge interval, write to file
  infile <- tempfile()
  write.table(to_write, infile, quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  command=paste( "sort -k1,1 -k2,2n", infile, "|", "/opt/bedtools/bin/bedtools merge -i -",">",outfile,sep=" ")
  cat(command,"\n")
  try(system(command))
}
writeBed(to_write, "Output/candidate_gene_hg19_exons.bed")

}
