summary stat on COSMIC
========================================================

Read in the 1000G snps_indel file. 

```{r}
setwd("/Users/snafu/Documents/Project/germVar")
load("/Volumes/Orchestra/seq/COSMIC/CosmicCodingMuts_v69.anno.RData")
source("./loadCommData.R")
```


Total number of variants `r dim(call_set)[1]`.   

Show the columns of the COSMIC
```{r}
summary(call_set)
```


Intersect the call_set with genes of interest, remove duplicates
```{r}
call_set <- call_set[!is.na(call_set$Gene), ]
# rescue unrecognized Gene names
call_set<-rescueSnpEffGenes(call_set)
#cosm_goi <- merge(call_set, list_goi, by='Gene')
cosm_goi <- subset(call_set, Gene %in% list_goi$Gene  & EFF == "NON_SYNONYMOUS_CODING")
rm(call_set)
# creat unique identifier for variant, remove duplicates
cosm_goi$var_uid <- apply(cosm_goi[c("Gene", "POS", "REF", "ALT")], 1, function(x) gsub(" ","",paste0(x, collapse="-")))
cosm_goi <- arrange(cosm_goi, var_uid, -CNT)
cosm_goi <- cosm_goi[!duplicated(cosm_goi$var_uid), ]

# remove LOW Impact, synonymous, intron, etc
#cosm_goi <- droplevels(subset(cosm_goi, !Impact %in% c("LOW", "MODIFIER")))
# remove complex variants
#cosm_goi <- droplevels(subset(cosm_goi, !(VARTYPE=="MNP" & MNP!=1)))
print(nrow(cosm_goi))

# define the AA_pos for splice site mutations 
#cosm_goi$AA_pos[with(cosm_goi, EFF %in% c("SPLICE_SITE_DONOR", "SPLICE_SITE_ACCEPTOR"))] <-  
#  apply(cosm_goi[c("EFF", "ExonRank")][with(cosm_goi, EFF %in% c("SPLICE_SITE_DONOR", "SPLICE_SITE_ACCEPTOR")),], 1, function(x) ifelse(x[1]=="SPLICE_SITE_DONOR", paste("D", x[2], sep=""), paste0("A", x[2], collapse="-")))
# tally coding variants position
cosm_goi$site <- interaction(cosm_goi$Gene, cosm_goi$AA_pos, drop=T)
#cosm_goi$VariantType <- factor(sapply(cosm_goi$EFF, function(x) ifelse(x=="NON_SYNONYMOUS_CODING", "Nonsyn", "Other")))
#cosm_goi$site_id <- interaction(cosm_goi$site, cosm_goi$VariantType, drop=T)
#soma_goi <- merge(soma_goi, dplyr::summarise(group_by(soma_goi, Gene, site_id), soma_CNT = sum(count))) 

#cosm_goi<-subset(cosm_goi, CNT>1)
site_tally <- dplyr::summarise(group_by(cosm_goi, site), cosm_count = sum(CNT)) 
cosm_goi <- merge(cosm_goi, site_tally, by="site")

site_tally$count <- cut(site_tally$cosm_count, c(0:30,50000)) 
levels(site_tally$count) <- c(as.character(1:30), ">30")
p2 <- ggplot(site_tally, aes(count)) + geom_bar(stat = "bin", width = .7,  position="dodge") + scale_y_sqrt()
p2 <- p2 + theme(panel.background = element_rect(fill='white', colour='black'))
p2 <- p2 + theme( axis.text.x = element_text(size= rel(1.), angle=90, hjust =1., vjust=.5), axis.title = element_text(size= rel(1.))) 
p2 <- p2 + ylab("Count") + xlab("Number of times reported")
p2

#recur_cosm <- subset(cosm_goi, (VariantType=="Nonsyn" & cosm_CNT>=8) |(VariantType!="Nonsyn" & cosm_CNT>=5) )
#recur_cosm <- subset(cosm_goi, site %in% recur_sites$site)
recur_cosm <- droplevels(subset(cosm_goi, cosm_count>=5) )

save(recur_cosm, file="./Results/COSMIC_recurrent_variants.RData")
recur_cosm_bed <- recur_cosm[c("CHROM", "POS")]
recur_cosm_bed$BEG <- recur_cosm$POS -1
write.table(recur_cosm_bed[c("CHROM", "BEG", "POS")], "./output/COSMIC_recurrent_mutations.bed", quote=F, row.names=F, col.names=F, sep="\t")
cosm_goi <- cosm_goi[c("Gene", "var_uid", "CNT", "AAChange")]
save(cosm_goi, file="./output/COSMIC_all_variants.RData")
rm(site_tally, cosm_goi)

```



