
setwd("/Users/snafu/Documents/Project/germVar")

library(VariantAnnotation)
# variant info
cosm_goi <- expand(readVcf(file="./DataSet_Results/CosmicCodingMuts.v71.goi.vcf.gz", "hg19"))

# preprocess
cosm_goi <- cosm_goi %>% labelUid   %>%   
  subset(., !is.na(Coding) & Coding == "CODING" ) %>%
  subset(., !is.na(BioType) & BioType == "protein_coding" ) %>%
  subset(., !Impact %in% c("LOW", "MODIFIER") ) %>%
  fixSnpEffGene %>%
  subset(., Gene %in% list_goi$Gene) %>%
  labelVarUid %>%
  labelSnpEffUid %>%
  labelAAUid %>% 
  labelSite 

cosm_goi <- as.data.frame(info(cosm_goi))
# too  big for arrange
cosm_goi <- cosm_goi %>% arrange(., var_uid, -CNT) %>% subset(., !duplicated(var_uid)) %>% plyr::rename(., replace=c("CNT"="cosm_vcount"))

# put away  lof variants
cosm_LoF <- subset(cosm_goi, Impact =="HIGH") 

# tally nonsyn coding variants position
cosm_nsSNP <- subset(cosm_goi, VARTYPE %in% c("SNP", "MNP") & EFF =="missense_variant") 
site_tally <- dplyr::summarise(group_by(cosm_nsSNP, Gene, site), cosm_scount=sum(cosm_vcount)) #%>% subset(., cosm_scount>=3)
aa_tally <- dplyr::summarise(group_by(cosm_nsSNP, Gene, aa_uid), cosm_acount=sum(cosm_vcount)) %>% mutate(., cosm_arank = rank(-cosm_acount, ties.method="average"))
cosm_nsSNP <- cosm_nsSNP %>% plyr::join(., site_tally, by="site") %>% plyr::join(., aa_tally, by="aa_uid")


#clin_site <- subset(site_tally, site %in% clinvar_txt$site)
# 
# if(F){
# # take subset of sites, count>=3
# #site_tally <- dplyr::summarise(group_by(cosm_nsSNP, Gene, site), cosm_scount=sum(cosm_vcount)) #%>% subset(., cosm_scount>=3)
# gene_select <- group_by(site_tally, Gene) %>% dplyr::summarise(., N_site = length(site), N2_site = length(site[cosm_scount>=2]), N3_site = length(site[cosm_scount>=3]), N4_site = length(site[cosm_scount>=4]), N5_site = length(site[cosm_scount>=5])) %>%  
#               subset(., N3_site>0)  %>% plyr::join(., subset(cosm_nsSNP, !duplicated(Gene))[c("Gene", "AALength")]) %>%
#               mutate(., N_ratio = N_site/AALength, N3_ratio = N3_site/N_site, N_z = (N_ratio-median(N_ratio))/IQR(N_ratio)*1.349, N3_z = (N3_ratio-median(N3_ratio))/IQR(N3_ratio)*1.349) %>% arrange(., N3_ratio)
# # top mutator
# gene_select <- subset(gene_select , !Gene %in% c("TP53", "CDKN2A", "VHL", "PTEN", "KRAS", "PIK3CA"))
# subset(gene_select, Gene %in% c("FAT4", "PEG3", "FLG", "FLT3", "ATM", "XIRP2", "PARP1", "POLE", "RET", "DICER1", "PTPRJ", "BRIP1", "MPO", "SDHB"))
# 
# select_thresh <- function(N_z, N3_z){
#    if( N_z>=2 & N3_z>=2){
#      return(5)
#    } else if( N_z <0 | N3_z < 0 | (N_z <2 & N3_z <2)){
#      return(3)
#    } else {
#      return(4)
#    }     
# }
# 
# gene_select$thresh <- mapply(select_thresh, gene_select$N_z, gene_select$N3_z)              
# cosm_nsSNP <- merge(cosm_nsSNP, gene_select[c("Gene", "thresh")])
# cosm_nsSNP <- subset(cosm_nsSNP, cosm_scount>= thresh)
# }

cosm_nsSNP <- subset(cosm_nsSNP, cosm_scount>= 3)

### plot the gene tally, distribution of MHS per gene
gene_tally <- as.data.frame(table(site_tally$Gene)) %>% plyr::rename(., replace=c("Var1"="Gene"))
gene_tally$scount_bin <- cut(gene_tally$Freq, c(0:20,500)) 
levels(gene_tally$scount_bin) <- c(as.character(1:20), ">20")
p <- ggplot(gene_tally, aes(scount_bin)) + geom_bar(stat = "bin", width = .7,  position="dodge", fill="darkgrey") 
p <- p + theme_few()
p <- p + theme( axis.text.x = element_text(size= rel(.7), angle=0, hjust =.5), axis.title = element_text(size= rel(.7))) 
p <- p + ylab("Number of genes") + xlab("Number of MHS per gene")
p
ggsave(filename="COSMIC_mhs_gene_hist.png",height=4, width=4)

to_plot <- plyr::join(site_tally, list_goi[c("Gene", "mem")], by="Gene")
to_plot$mem <- factor(to_plot$mem, levels=c("OTHER", "REP", "TSG", "SMG", "CGC"))
p <- ggplot(to_plot, aes(mem)) + geom_bar(stat = "bin", width= 0.7, fill="darkgrey") + theme_few()
p <- p  + theme( axis.title= element_text(size=rel(0.7)), axis.text.y= element_text(size=rel(0.7)))
p <- p + scale_x_discrete(labels = c( "CGC"="Cancer Gene\nCensus", "SMG"="Significantly\nMutated Genes", "TSG" ="Tumor Suppressor\nGene Database", "REP"="Human DNA\nRepair Genes", "OTHER"="Others"))
p <- p + coord_flip() + xlab("") + ylab("Number of MHS") 
p
ggsave(filename="COSMIC_mhs_mem_hist.png",height=4, width=4)


snp2BED <- function(df){
  to_write <- data.frame(CHROM = df$CHROM )
  PHASES <- sapply(df$CDS_pos, function(x) x %% 3)
  to_write$START <- mapply(function(POS, PHASE, strand) ifelse(strand==1, POS - ((PHASE+2) %%3) - 1, POS - ((4-PHASE) %%3)), df$POS, PHASES, df$strand)
  to_write$END <- mapply(function(POS, PHASE, strand) ifelse(strand==1, POS + ((3-PHASE) %%3), POS + ((PHASE+2) %%3)), df$POS, PHASES, df$strand)
  return(to_write)
}

tmp <- rbind(cosm_nsSNP[c("CHROM", "POS", "Transcript", "CDS_pos", "site")],
      as.data.frame(info(clinvar_vcf)[c("CHROM", "POS", "Transcript", "CDS_pos", "site")]))
tmp <- plyr::join(tmp, list_gff_long[c("Transcript", "strand")], by="Transcript")

writeBed(snp2BED(tmp), "Output/nsSNP.MHS_Clinvar.bed")
#writeBed(snp2BED(subset(cosm_nsSNP, cosm_scount>=3 & Gene=="RB1")), "RB1.bed")

save(cosm_goi, file="./DataSet_Results/COSMIC_all_variants.RData")
rm(aa_tally, site_tally, cosm_goi)

```

calcEntro <- function(df){
  df$cosm_vcount[df$cosm_vcount>=50] <- 50
  site_tally <- dplyr::summarise(group_by(df, site), cosm_scount=sum(cosm_vcount))
  total <- sum(site_tally$cosm_scount)
  prob <- site_tally$cosm_scount/total
  entropy <- -(sum(prob*log(prob)))
  #hist(df$cosm_scount)
  #df$cuts <- cut(df$cosm_scount, breaks=c(0, 1.5, 2.5, Inf))
  #table(df$cuts)
  return(data.frame(Entropy =entropy ))#, Entropy0 = log(nrow(site_tally)), EntropyD = log(nrow(site_tally))-entropy))
}





