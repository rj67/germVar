
setwd("/Users/snafu/Documents/Project/germVar")
load("/Volumes/Orchestra/seq/COSMIC/CosmicCodingMuts_v69.goi.RData")
#source("./loadCommData.R")

cosm_goi <- call_set %>% subset(., !is.na(Gene))  %>%
            fixSnpEffGene %>%   # rescue unrecognized Gene names
            subset(., Gene %in% list_goi$Gene  & Impact %in% c("MODERATE", "HIGH"))  %>%  # remove LOW Impact, synonymous, intron, etc
            labelVarUid %>% arrange(., var_uid, -CNT) %>% subset(., !duplicated(var_uid)) %>%# creat unique identifier for variant, remove duplicates
            rename(., replace=c("CNT"="cosm_vcount"))

# tally nonsyn coding variants position
cosm_nonsyn <- subset(cosm_goi, EFF=="NON_SYNONYMOUS_CODING")
cosm_nonsyn$site <- apply(cosm_nonsyn[c("Gene", "AA_pos")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))
cosm_nonsyn$aa_uid <- apply(cosm_nonsyn[c("Gene", "AAChange")],1, function(x) gsub(" ", "", paste0(x, collapse="-")))

site_tally <- dplyr::summarise(group_by(cosm_nonsyn, site), cosm_scount = sum(cosm_vcount)) 
cosm_nonsyn <- merge(cosm_nonsyn, site_tally, by="site")

site_tally$count <- cut(site_tally$cosm_count, c(0:30,50000)) 
levels(site_tally$count) <- c(as.character(1:30), ">30")
require(ggplot2)
p2 <- ggplot(site_tally, aes(count)) + geom_bar(stat = "bin", width = .7,  position="dodge") + scale_y_sqrt()
p2 <- p2 + theme(panel.background = element_rect(fill='white', colour='black'))
p2 <- p2 + theme( axis.text.x = element_text(size= rel(1.), angle=90, hjust =1., vjust=.5), axis.title = element_text(size= rel(1.))) 
p2 <- p2 + ylab("Count") + xlab("Number of times reported")
p2
rm(p2)
recur_cosm <- droplevels(subset(cosm_nonsyn, cosm_scount>=5) )

save(recur_cosm, file="./Results/COSMIC_recurrent_variants.RData")

save(cosm_goi, file="./Output/COSMIC_all_variants.RData")
rm(site_tally, cosm_goi)

```



