
setwd("/Users/snafu/Documents/Project/germVar")
load("/Volumes/Orchestra/seq/COSMIC/CosmicCodingMuts_v70.goi.RData")
#source("./loadCommData.R")

call_set <- call_set[!is.na(call_set$Gene), ]
# rescue unrecognized Gene names
call_set<-rescueSnpEffGenes(call_set)
# remove LOW Impact, synonymous, intron, etc
cosm_goi <- subset(call_set, Gene %in% list_goi$Gene  & Impact %in% c("MODERATE", "HIGH"))
rm(call_set)
# creat unique identifier for variant, remove duplicates
cosm_goi <- createVarUid(cosm_goi)
cosm_goi <- arrange(cosm_goi, var_uid, -CNT)
cosm_goi <- cosm_goi[!duplicated(cosm_goi$var_uid), ]

# tally nonsyn coding variants position
cosm_nonsyn <- subset(cosm_goi, EFF=="NON_SYNONYMOUS_CODING")
cosm_nonsyn$site <- interaction(cosm_nonsyn$Gene, cosm_nonsyn$AA_pos, drop=T)
site_tally <- dplyr::summarise(group_by(cosm_nonsyn, site), cosm_count = sum(CNT)) 
cosm_nonsyn <- merge(cosm_nonsyn, site_tally, by="site")

site_tally$count <- cut(site_tally$cosm_count, c(0:30,50000)) 
levels(site_tally$count) <- c(as.character(1:30), ">30")
require(ggplot2)
p2 <- ggplot(site_tally, aes(count)) + geom_bar(stat = "bin", width = .7,  position="dodge") + scale_y_sqrt()
p2 <- p2 + theme(panel.background = element_rect(fill='white', colour='black'))
p2 <- p2 + theme( axis.text.x = element_text(size= rel(1.), angle=90, hjust =1., vjust=.5), axis.title = element_text(size= rel(1.))) 
p2 <- p2 + ylab("Count") + xlab("Number of times reported")
p2

recur_cosm <- droplevels(subset(cosm_nonsyn, cosm_count>=5) )

save(recur_cosm, file="./Results/COSMIC_recurrent_variants.RData")

save(cosm_goi, file="./Output/COSMIC_all_variants.RData")
rm(site_tally, cosm_goi)

```



