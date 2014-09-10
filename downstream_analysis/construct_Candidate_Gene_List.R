setwd("/Users/snafu/Documents/Project/germVar")
table_HGNC_file <- "../dataDump/HGNC/HGNC_ID_table_July_17.txt"
table_cgc500_file <- "./Gene_List/cancer_gene_census_v70.csv"
table_tsg716_file <- "../dataDump/TSGene/Human_716_TSGs.txt"
new_tsg145_file <- "./Gene_List/New_145_TSGs.txt"
hc_tsg206_file <- "./Gene_List/High_confidence_206_TSGs.txt"
table_smg260_file <- "./Gene_List/SMG_genes_260.csv"

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
# C12orf9 is retracted in HGNC 
# RUNDC2A have  outdated GeneID
# STL six twelve leukemia is a RNA gene, not in the HGNC table
subset(table_cgc500, Entrez.GeneId %in% setdiff(table_cgc500$Entrez.GeneId, table_HGNC$Entrez.Gene))
# fix RUNDC2A/DUX4
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
# merge the HGNC gene symbols with the tsgene table
table_tsg716 <- merge(table_tsg716, table_HGNC[c("Approved.Symbol", "Entrez.Gene", "UniProt.ID")], by.x = "GeneID", by.y = "Entrez.Gene")
# 5 genes with different names in 716 vs HGNC
table_tsg716[c("Gene_symbol", "Approved.Symbol")][!apply(table_tsg716[c("Gene_symbol", "Approved.Symbol")],1, function(x) {x[1]==x[2]}),]

# additional 154
new_tsg145 <- read.table(new_tsg145_file, strip.white = T, header =T)
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
table_smg260 <- read.csv(, stringsAsFactors =F)
# a few gene names are outdated
setdiff(table_smg260$gene, table_HGNC$Approved.Symbol)
table_smg260$gene[table_smg260$gene=='SFRS2'] <- 'SRSF2'
table_smg260$gene[table_smg260$gene=='MLL'] <- 'KMT2A'
table_smg260$gene[table_smg260$gene=='MLL2'] <- 'KMT2D'
table_smg260$gene[table_smg260$gene=='MLL3'] <- 'KMT2C'
table_smg260$gene[table_smg260$gene=='MLL4'] <- 'KMT2B'




# # -------------------------------------------------------------------------------------
# #  Giovanni Pan-Can paper significant function event (SFE)
# # -------------------------------------------------------------------------------------
#table_sfe <- read.csv("../dataDump/GiovanniSFE/ng.2762-S2.csv", stringsAsFactors =F)
#table_sfe <- droplevels(subset(table_sfe, Type == "MUTATION"))
# table_sfe <- droplevels(subset(table_sfe, Type %in% c("MUTATION", "METHYLATION")))
#colnames(table_sfe)[1] <- "Gene"
#table_sfe$Gene <- levels(table_sfe$Gene)[table_sfe$Gene]
#table_sfe$Gene[table_sfe$Gene=='C3orf63'] <- 'FAM208A'
#table_sfe$Gene[table_sfe$Gene=='FAM38B'] <- 'PIEZO2'
#table_sfe$Gene[table_sfe$Gene=='FAM5C'] <- 'BRINP3'
#table_sfe$Gene[table_sfe$Gene=='HEATR7A'] <- 'MROH1'
#table_sfe$Gene[table_sfe$Gene=='JUB'] <- 'AJUBA'
#table_sfe$Gene[table_sfe$Gene=='KIAA1804'] <- '??'
#table_sfe$Gene[table_sfe$Gene=='LPPR3'] <- '??'
#table_sfe$Gene[table_sfe$Gene=='MLL'] <- 'KMT2A'
#table_sfe$Gene[table_sfe$Gene=='MLL2'] <- 'KMT2D'
#table_sfe$Gene[table_sfe$Gene=='MYCL1'] <- 'MYCL'
#table_sfe$Gene[table_sfe$Gene=='ZNF498'] <- 'ZSCAN25'
table_sfe <- read.delim(file="input/pancan16_mini_maf_annotated.txt", stringsAsFactors=F)
colnames(table_sfe)[1] <- "Gene"
table_sfe$Gene[table_sfe$Gene=='MLL2'] <- 'KMT2D'
table_sfe$Gene[table_sfe$Gene=='SFRS2'] <- 'SRSF2'

# -------------------------------------------------------------------------------------
# Various Gene familis from HGNC
# -------------------------------------------------------------------------------------
fam_files <- list.files(path = "./Gene_List/", pattern="_family.txt", full.names=F)
fam_names <- gsub("_family.txt", "", fam_files, fixed=T)
fam_genes <- do.call(rbind, lapply(fam_names, function(fam) {
  fam_table <- read.delim(file=paste("./Gene_List/", fam, "_family.txt", sep=""), head=T, strip.white=T, stringsAsFactors=F);
  returns <- data.frame("Gene" = fam_table[["Approved.Symbol"]]);
  returns$fam_id <- fam;
  print(dim(returns));
  return(returns)
}))
# Pro-apoptotic factors
BH3 <- c( "BAX", "BAK1", "BOK", "BID", "BAD", "BIK", "BBC3", "PMAIP1", "BMF", "HRK")
# Nu-A4 complex
NUA4 <- c( "TRRAP",  "EP400", "BRD8", "EPC1", "KAT5", "DMAP1", "ING3", "VPS72", "RUVBL1", "RUVBL2", "MORF4L1", "YEATS4", "MRGBP", "MEAF6")
# PRC2
PRC2 <- c( "SUZ12", "EED", "EZH1", "EZH2", "RBBP4", "RBBP7", "AEBP2", "JARID2", "PHF1")
# NuRD
NURD <- c( "CHD3", "CHD4", "HDAC1", "HDAC2", "MBD2", "MBD3", "MTA1", "MTA2", "MTA3", "GATAD2A", "GATAD2B", "RBBP4", "RBBP7")
# SWISNF
SWISNF <- c( "ARID1A", "ARID1B", "ARID2", "SMARCA2", "SMARCA4", "PBRM1", "SMARCC2", "SMARCC1", "SMARCD1", "SMARCD2", "SMARCD3", "SMARCE1", "ACTL6A", "ACTL6B", "SMARCB1", "BRD7")
# Pre-mRNA processing
mRNA <- c("PRPF40B", "SF1", "ZRSR2", "U2AF1", "U2AF2", "SRSF2", "SF3B1", "SF3A1")
# miRNA processing
miRNA <- c("DICER1", "AGO1", "AGO2", "TARBP2", "DROSHA") 
# manual literature curation
LIT <- c("BCORL1", "PRDM5", "DNMT1","DNMT3B", "HAND2", "ADAMTS15", "FGFR4",
              "HIST1H1B", "HIST1H1C", "ELP4", "BAG6", "CHD1", "PIK3R2", "TRIP12", 
               "BUB1", "BUB3","LZTR1", "PLCG1", "SYNE1", 
              "PREX2", "SOX9", "AMER2", "LDLRAP1", "STMN2", "AGTR2", "ZIM2", "NALCN", "SLC16A4", "MAGEA6", "RPS6KA3", "IRF2", "ROBO2", 
              "HNRNPK", "BRINP3", "KLHL6", "MEF2B", "ACTB", "PCLO", "LRRK2",
              "MTF1", "GALNT1", "TINF2", "WRAP53", "ELANE", "DKC1", "EPCAM" )# to add "FBP1")

list_add <- list(BH3 = BH3, NUA4 = NUA4, PRC2 = PRC2, NURD = NURD, SWISNF = SWISNF, mRNA = mRNA, miRNA = miRNA, LIT=LIT)

fam_genes <- rbind(fam_genes, do.call(rbind, lapply(c("BH3", "NUA4", "PRC2", "NURD", "SWISNF", "mRNA", "miRNA", "LIT"), function(x) data.frame(Gene=list_add[[x]], fam_id = x))))


# -------------------------------------------------------------------------------------
# Various DNA repair genes 
# -------------------------------------------------------------------------------------
rep_files <- list.files(path = "./Gene_List/", pattern="_REP.txt", full.names=F)
rep_names <- gsub("_REP.txt", "", rep_files, fixed=T)
rep_genes <- do.call(rbind, lapply(rep_names, function(path) {
  rep_table <- read.delim(file=paste("./Gene_List/", path, "_REP.txt", sep=""), skip = 1, head=F, strip.white=T, stringsAsFactors=F);
  returns <- data.frame("Gene" = rep_table[["V1"]]);
  returns$rep_id <- path;
  return(returns)
}))

# -------------------------------------------------------------------------------------
#list_hr <- c("BRCA1", "BRCA2", "FAM175A", "TP53", "TP63", "CHEK1", "CHEK2", "ATM","ATR","NBN", "RAD50", "BRIP1", "PALB2", "BARD1", 
#             "MRE11A", "MSH6", "RAD51C", "RAD51B","RAD51D","PMS1","PTEN","BAX","MLH1","MSH2","PMS2","MUTYH","CDH1","STK11",
#             "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF","FANCG" ,"FANCI" ,"FANCL", "FANCM","FANCP","XPF"   )
# Homologous Recombination and Fanconi Anemia group
# -------------------------------------------------------------------------------------
#list_hr <- c( "FAM175A", "NBN", "RAD50",  "BAP1",  "BARD1",  "MRE11A", "RAD51B", "RAD51D", 
#              "RAD51C","SLX4", "BRCA1", "BRCA2", "BRIP1", "FANCA", "PALB2", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF","FANCG" ,"FANCI" ,"FANCL", "FANCM")

#list_hr_ext <- c( "XRCC2", "XRCC3", "MDC1", "DMC1", "RAD51", "RAD51AP1", "BRAP", "BRCC3", "FAN1", "FOXM1", "MUS81" )



# -------------------------------------------------------------------------------------
#  Gene mania gene list of various pathways 
# -------------------------------------------------------------------------------------
#table_sting <- read.delim("./input/GeneMania_sting_genes_list.txt", stringsAsFactors=F, strip.white=T)
#colnames(table_sting)[1] <- "Gene"


# -------------------------------------------------------------------------------------

# genes in S500 not S200 but should be in S200
#list_s200_ext <-c("SUZ12", "CBFB", "CCND1", "CCND3", "ETV6", "ITK", "NTRK1", "NTRK3", "RUNX1", "TCF7L2", "SEPT9", "FGFR1", "PER1", "RAD51B" )
#the reverse
#list_s200_exc <-c("KIAA1549")

# -------------------------------------------------------------------------------------
# list of genes with activating mutation or overexpression that leads to tumorgenesis
list_active <- c("EGFR", "FGFR1", "FGFR2", "FGFR3", "GLDC", "GNAS", "KRAS", "HRAS", "NRAS", "IDH1", "IDH2", "ALK", "ERBB2", "ERBB3", "BRAF", "PIK3CA", "PIK3R1", "PIK3R2", "CTNNB1", "AKT1", "MAP2K1", "MYC",
                 "FLT3", "IL7R",  "PTPN11", "KDR", "FGFR4","JAK1", "JAK2", "JAK3", "MYCN", "CRCLF2", "ABL1", "NT5C2", "PDGFRA", "MAP2K2", "MPL", "MET", "SETBP1", "TRRAP", "MYD88", "KIT", "ITK",
                 "ALK", "NOTCH1", "NOTCH2", "CARD11", "ECT2L", "FOXL2", "CBL", "KLF4", "SMO", "STAT3", "STAT5B","MTOR", "BCL6", "XPO1", "TSHR", "U2AF1", "SF3B1", "SRSF2", "GNAQ", "GNA11", "NTRK3" ,"NTRK1", 
                 "RAC1", "SEPT9", "TERT", "TBL1XR1", "CDK4", "CDK12", "CCND1", "CCND3", "CD79B", "CALR", "CRLF2", "CACNA1D", "ATP1A1", "DNMT3A", "KCNJ5", "NFE2L2", "NPM1", "PHOX2B", "PPP2R1A",
                 "H3F3A", "H3F3B", "GATA1", "CBLB", "RET", "TRAF7", "FOXA1", "DNM2", "RPL5")

# list of genes suspected TS
#list_ts <- c("ATRX", "CIC", "BCOR", "BIRC3","BMPR1A", "CBFB", "CNOT3", "DAXX", "DICER1", "EP300", "EZH2", "FBXO11", "FUBP1", "KDM5C", "KMT2A", "KMT2C", "KMT2D", "POT1", "RNF43", "RUNX1", "SDHC", "SETD2", "SH2B3", "STAG2", "SUZ12", "UBR5", "WAS", "WRN", "ZRSR2",
#             "FAM46C", "FAS", "GATA2", "GATA3", "TNFRSF14", "HNF1A", "RECQL4", "RAD21", "KDM6A", "ARHGAP26")

# -------------------------------------------------------------------------------------
list_goi <-c(as.character(table_tsg716$Approved.Symbol), as.character(table_cgc500$Approved.Symbol), table_smg260$gene,  levels(rep_genes$Gene)[rep_genes$Gene],
             levels(fam_genes$Gene)[fam_genes$Gene], table_sfe$Gene)
list_goi <- unique(list_goi)
list_goi <- data.frame("Gene" =list_goi)
list_goi$cgc500 <- list_goi$Gene %in% as.character(table_cgc500$Approved.Symbol)
list_goi$cgc200 <- with(list_goi, cgc500 & !(Gene %in% subset(table_cgc500, Mutation.Type %in% c("A", "T", "A,T", "T,A"))$Approved.Symbol))
#list_goi$Dom <- list_goi$Gene %in% as.character(with(table_cgc500, Approved.Symbol[grep('Dom', Cancer.Molecular.Genetics)]))
list_goi$tsg716 <- list_goi$Gene %in% as.character(table_tsg716$Approved.Symbol)
list_goi$tsg206 <- list_goi$Gene %in% as.character(subset(table_tsg716, High_conf)$Approved.Symbol)
list_goi$smg260 <- list_goi$Gene %in% table_smg260$gene
list_goi$smg33 <- list_goi$Gene %in% subset(table_smg260,Discussed.in.text..33.==1)$gene
list_goi$sfe160 <- list_goi$Gene %in% table_sfe$Gene
list_goi$rep <- list_goi$Gene %in% rep_genes$Gene
list_goi <- merge(list_goi, rep_genes, by="Gene", all.x=T)
list_goi$fam <- list_goi$Gene %in% fam_genes$Gene
list_goi <- merge(list_goi, fam_genes, by="Gene", all.x=T)

list_goi <- subset(list_goi, ! Gene %in% c("BCL5", "TCL6", "IGH", "IGL", "IGK", "TRA", "TRB", "TRD", "HLA-A", "HLA-B"))
#list_goi$sfe <- list_goi$Gene %in% table_sfe$Gene
#list_goi$active<- list_goi$Gene %in% list_active
#list_goi$bh3<- list_goi$Gene %in% list_bh3
#list_goi$nua4<- list_goi$Gene %in% list_nua4
#list_goi$mirna<- list_goi$Gene %in% list_mirna
#list_goi$ts<- list_goi$Gene %in% list_ts
#require(reshape2)

#list_goi$exp <- with(list_goi, !rep &(bh3 | nua4 | mirna | CASP | EPH | EPHRIN | KAT | KDM | KMT | MAP2K | MAP3K | MAP4K | MAPK |SRSF | TNFRSF |TNFSF))
#list_goi$ts2 <- with(list_goi, !s200 & (ts206 | ts716 & !ts206 & (s500 | smg260)))

#with(list_goi, sum(bh3|nua4| mirna| CASP|DNAREF|DNARIF|DNARTF|EPH | EPHRIN | KAT|KDM|KMT|MAP2K|MAP3K|MAP4K|MAPK|PARP|SRSF|TNFRSF|TNFSF))

#get the full name
list_goi <- merge(list_goi, table_HGNC[c("Approved.Symbol", "Approved.Name", "Entrez.Gene", "UniProt.ID", "Ensembl.ID", "Locus.Type", "Locus.Group")], by.x = "Gene", by.y = "Approved.Symbol")

# remove RNA genes and pseudogene
print( subset(list_goi, Locus.Group %in%c("non-coding RNA", "pseudogene")))
list_goi <- subset(list_goi, ! Locus.Group %in%c("non-coding RNA", "pseudogene"))
list_goi$Locus.Type <- NULL
list_goi$Locus.Group <- NULL

#list_goi$s200[list_goi$Gene %in% list_s200_ext] <- TRUE 
#list_goi$s200[list_goi$Gene %in% list_s200_exc] <- FALSE
# sfe only, mostly crap
#list_goi$sfeonly <- with(list_goi, sfe&!(s500|hr|ts716|smg260))
#remove frequently deleted and unmappable HLA-A/B

# -------------------------------------------------------------------------------------
#  Tumor Associated Gene 700 genes
# -------------------------------------------------------------------------------------
table_tag700 <- read.csv("./Gene_List/Tumor_Associated_Genes.csv", stringsAsFactors =F)
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

list_goi$TS <- list_goi$Gene %in% subset(table_tag700, Category=="Tumor suppressor gene")$Symbol
list_goi$ONCO <- list_goi$Gene %in% subset(table_tag700, Category=="Oncogene")$Symbol
list_goi$ONCO[list_goi$Gene %in% c("CBLB", "CBLC","GNA11", "GNAQ", "IDH1", "IDH2", "IL7R", "JAK1", "JAK2", "JAK3", "MYD88", "NOTCH1", "NOTCH2", "TSHR", "PIK3R1", "ALK")] <- T
list_goi$TS[with(list_goi, cgc200 & !ONCO & !TS & rep)] <- T
list_goi$TS[with(list_goi, cgc200 & !ONCO & !TS & tsg206)] <- T
list_goi$TS[list_goi$Gene %in% c("ATRX", "BUB1B", "DAXX", "EP300", "EZH2", "KMT2C", "KMT2D", "BMPR1A", "POT1", "RAD51B")] <- T

#setdiff(subset(table_tag700, Category=="Tumor suppressor gene")$Symbol, list_goi$Gene)
#setdiff(subset(table_tag700, Category=="Tumor suppressor gene")$Symbol, table_tsg716$Approved.Symbol)

#list_goi <- subset(list_goi, !Gene %in% c("HLA-B", "HLA-A", "PDE4DIP", "IGH", "IGK", "IGL", "TRA", "TRB", "TRD"))
write.table(list_goi, file="output/candidate_gene_list.tsv", quote=F, row.names=F, sep="\t")
rm(fam_genes, rep_genes, table_cgc500, table_smg260, table_sfe, table_tsg716, table_tag700, tsgene_206, new_145, LIT, miRNA, mRNA, NUA4, NURD, PRC2, SWISNF, fam_files, 
   fam_names, rep_files, rep_names, list_add, list_active, BH3, table_tag700_in, table_tag700_out)

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

#write stratton 500 exons
list_goi_exons<-exons(hg19.refgene.tx, list(gene_id = list_goi$Entrez.Gene),columns=c("gene_id","exon_id"))
#convert Grange object to 
Grange2bed<-function(GrangeObject){
  #require(dplyr)
  returns<-as.data.frame(GrangeObject)
  returns$seqnames<-sapply(returns$seqnames,function(x) gsub("chr","",x))
  #returns<- arrange(returns, seqnames, start)
  return(returns)
  #returns<-returns[order(returns$seqnames,returns$start),]
}
write.table(Grange2bed(list_goi_exons),file="output/candidate_gene_hg19R_exons.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

}

