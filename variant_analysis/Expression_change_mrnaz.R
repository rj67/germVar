
#####################################################
## LoFl low frequency variants, per variant mrnaz wilcox test
#####################################################

LoF_lf_mrnaz <- group_by(droplevels(subset(LoF_GTj, uid %in% LoF_lf$uid & ToN=="N" & is_varq & !Gene %in% c("CEP89", "CNBD1"))), uid) %>% do(mrnaz_test(.))
# check the minimum possible N_var to result in pval < 0.05
LoF_lf_mrnaz<-arrange(plyr::join(LoF_lf_mrnaz, LoF_lf[c("Gene", "EFF", "uid", "AAChange", "tier1", "tier2", "NMD_HC", "LoF_filter", "LOF_PT")], by="uid"), pval)
LoF_lf_mrnaz %<>% plyr::join(., list_goi[c("Gene", "driver", "Cat", "tsg716", "tsg206", "rep")], by="Gene") 

to_plot %<>% subset(., N_var>=8) %>% mutate(., padj = p.adjust(pval, method="BH") )
to_plot$pred <- with(to_plot, mapply(function(LoF, NMD) ifelse(LoF, ifelse(NMD, "LoF/NMD", "LoF"), "None" ), tier1, NMD_HC))
to_plot$pred <- factor(to_plot$pred, levels=c("LoF/NMD", "LoF", "None"))
p <- ggplot(to_plot, aes(effect, pval, label=Gene)) + geom_point(alpha=0.6, aes(size= sqrt(N_var), color=pred)) +  scale_size_continuous(range = c(2, 8), guide="none")
p <- p + scale_x_continuous(limits=c(-1.7, .4)) + scale_color_manual(values=c("#f46d43", "#74add1", "grey")) 
p <- p + geom_text(data = subset(to_plot, pval <=0.001), position=position_jitter(height=0.35, width=0.15), size=3, color="black", alpha=0.9)
p <- p + theme_few() + scale_y_continuous(trans=reverselog_trans(10))
p <- p + geom_hline(yintercept=0.015, color="grey", linetype=3) + geom_text(x=.4, y = 2., label="FDR 10%", hjust=1, color="grey", size=3)
p <- p + geom_hline(yintercept=0.001, color="grey", linetype=3) + geom_text(x=.4, y = 3.3, label="FDR 1%", hjust=1, color="grey", size=3)
p <- p + geom_vline(xintercept=0, color="lightgrey", linetype=3) 
p <- p + theme(axis.text.y=element_text(angle=90, hjust=0.5), axis.title=element_text(size=rel(0.7)))
p <- p + ylab("pval for altered mRNA level") + xlab("Average mRNA level change")
p <- p + theme(legend.position=c(0.87, 0.85), legend.title = element_blank())
p
ggsave(filename="LoF_lf_mrnaz_test.png",width=5, height=5)

#####################################################
## mrnaz wilcox test for rare LoF variants, per genes
#####################################################
LoF_rv_muts <- subset(LoF_GTj, ToN=="T" & is_var2 & uid %in% subset(LoF_rv, tier1)$uid)
LoF_rv_genes <- LoF_rv_muts %>% mutate(., event_uid = paste(Patient, Gene, sep="-")) %>% group_by(., Gene) %>% 
  dplyr::summarise(., N_pat = length(unique(event_uid))) %>% subset(., N_pat>=8)
LoF_rv_genes <- subset(LoF_rv_genes, !Gene %in% c("CEP89", "CNBD1","SEPT12", "CLIP1", "PPFIBP1"))
#LoF_rv_muts <- subset(LoF_rv_muts, Gene %in% LoF_rv_genes$Gene)
LoF_rv_muts %<>% arrange(., -DP) %>% subset(., !duplicated(mut_uid))

LoF_rv_mrnaz <- group_by(droplevels(subset(LoF_GTj, uid %in% subset(LoF_rv_muts, Gene %in% LoF_rv_genes$Gene)$uid & ToN=="T" & is_var2 )), Gene) %>% do(mrnaz_test(.))
LoF_rv_mrnaz %<>% plyr::join(., list_goi[c("Gene", "driver", "Cat", "tsg716", "tsg206", "rep")], by="Gene") 
# flag entries where not enough samples
LoF_rv_mrnaz$suff <- with(LoF_rv_mrnaz, N_var>=8 | (!is.na(pval) & pval <0.01))
LoF_rv_mrnaz$padj<-1
LoF_rv_mrnaz$padj[LoF_rv_mrnaz$suff ]<- p.adjust(LoF_rv_mrnaz$pval[LoF_rv_mrnaz$suff ], method="BH")

#trunc_mrnaz_test$sig <- cut(trunc_mrnaz_test$padj, breaks=c(0, 0.001, 0.01, 0.05, 1), labels=c("***", "**", "*", ""))
p <- ggplot(subset(LoF_rv_mrnaz, suff), aes(effect, pval, label=Gene)) + geom_point(alpha=0.6, aes(size= sqrt(N_var), color=Cat)) +  scale_size_continuous(range = c(2, 8), guide="none")
p <- p + scale_x_continuous(limits=c(-2.5, 2.)) + scale_color_manual(values=c("#f46d43", "grey", "#74add1")) 
p <- p + geom_text(data = subset(LoF_rv_mrnaz, padj <=0.002|driver), position=position_jitter(height=0.4, width=0.3), size=2, color="black", alpha=0.7)
p <- p + theme_few() + scale_y_continuous(trans=reverselog_trans(10))
p <- p + geom_hline(yintercept=0.03, color="grey", linetype=3) + geom_text(x=2., y = 1.5, label="FDR 10%", hjust=1, color="grey", size=3)
p <- p + geom_hline(yintercept=0.0015, color="grey", linetype=3) + geom_text(x=2., y = 3., label="FDR 1%", hjust=1, color="grey", size=3)
p <- p + geom_vline(xintercept=0, color="lightgrey", linetype=3) 
p <- p + theme(axis.text.y=element_text(angle=90, hjust=0.5), axis.title=element_text(size=rel(0.7)))
p <- p + ylab("pval for altered mRNA level") + xlab("Average mRNA level change")
p <- p + theme(legend.position=c(0.87, 0.85), legend.title = element_blank())
p
ggsave(filename="LoF_rv_mrnaz_test.png",width=5, height=5)


#####################################################
## mrnaz wilcox test for all trunc genes
#####################################################
gene_tally <- subset(gene_tally, ! Gene %in% c("SPO11", "SEPT12"))
#if(is.na(tryCatch(subset(RNASeq, Gene==x), error=function(cond) {message(cond);return(NA)}) )){
#test <- as.data.frame(t(sapply(intersect(gene_tally$Gene, subset(list_goi, cgc500)$Gene), mrnaz_test)))
trunc_mrnaz_test <- sapply(gene_tally$Gene, mrnaz_test)
trunc_mrnaz_test <- do.call(rbind, trunc_mrnaz_test)
trunc_mrnaz_test$padj <- p.adjust(trunc_mrnaz_test$pval, method="BH") 
trunc_mrnaz_test$sig <- cut(trunc_mrnaz_test$padj, breaks=c(0, 0.001, 0.01, 0.05, 1), labels=c("***", "**", "*", ""))

trunc_mrnaz <- group_by(droplevels(subset(trunc_MUTs, Gene %in% subset(trunc_gene, N_pat>=5)$Gene & !Gene %in% c("SPO11", "SEPT12"))), Gene) %>% do(mrnaz_test(.))
# check the minimum possible N_var to result in pval < 0.05
head(arrange(subset(trunc_mrnaz, pval<0.05), N_var))
trunc_mrnaz <- subset(trunc_mrnaz, N_var>=8)
trunc_mrnaz$padj <- p.adjust(trunc_mrnaz$pval, method="BH") 

View(arrange(trunc_mrnaz, pval))
#cgc_test <- subset(trunc_mrnaz_test, Gene %in% subset(list_goi, cgc500)$Gene & N_var >=5)
trunc_mrnaz <- plyr::join(trunc_mrnaz, list_goi[c("Gene", "driver", "Cat")], by="Gene")

p <- ggplot(trunc_mrnaz, aes(effect, pval, label=Gene)) + geom_point(alpha=0.6, aes(size= sqrt(N_var), color=driver)) +  scale_size_continuous(range = c(2, 8), guide="none")
p <- p + scale_x_continuous(limits=c(-2.3, 2.3)) + scale_color_manual(values=c("grey", "#d73027"), guide=F)
p <- p + geom_text(data = subset(trunc_mrnaz, pval <=0.00001 | driver), position=position_jitter(height=0.3, width=0.2), size=3, color="black", alpha=0.9)
p <- p + theme_few() + scale_y_continuous(trans=reverselog_trans(10))
p <- p + geom_hline(yintercept=0.027, color="grey", linetype=3) + geom_text(x=2.4, y = 1.8, label="FDR 10%", hjust=1, color="grey", size=3)
p <- p + geom_hline(yintercept=0.0012, color="grey", linetype=3) + geom_text(x=2.4, y = 3.1, label="FDR 1%", hjust=1, color="grey", size=3)
p <- p + geom_vline(xintercept=0, color="lightgrey", linetype=3) 
p <- p + theme(axis.text.y=element_text(angle=90, hjust=0.5), axis.title=element_text(size=rel(0.7)))
p <- p + ylab("pval for altered mRNA level") + xlab("Average mRNA level change")
p
ggsave(filename="trunc_mrnaz_test.png",width=5, height=7)


tmp <- subset(mut_stat, Gene %in% subset(list_goi, cgc500)$Gene)
tmp <- merge(tmp, list_goi[c("Gene", "TS", "ONCO")], by="Gene")
tmp$cat <- apply(tmp[c("TS", "ONCO")], 1, function(x) ifelse(x[1], "TS", ifelse(x[2], "ONCO", "OTHER")))
tmp2 <- arrange(dplyr::summarise(group_by(tmp, Gene), mut_count = length(unique(event_uid[PIT]))), -mut_count)
tmp <- subset(merge(tmp, tmp2, by="Gene"), mut_count >=10)
tmp <- tmp[!duplicated(tmp$mut_uid),]
#tmp <- merge(tmp, CBio_query[c("event_uid", "log2CNA",  "gistic", "mrnaz",  "study", "gistic2")], by="event_uid", all.x=T)
tmp <- merge(tmp, RNASeq[c("event_uid", "normalized_count", "mrnaz")], by="event_uid", all.x=T)
tmp <- merge(tmp, CBio_query[c("event_uid", "log2CNA",  "gistic",  "gistic2")], by="event_uid", all.x=T)
pvals <- sapply( unique(tmp$Gene),  function(x){
  print(x);
  to_test <- subset(tmp, Gene==x&PIT);
  bg <- subset(RNASeq, Gene==x)
  return(wilcoxTestCallSet(to_test, "mrnaz", bg))
})
pvals <- as.data.frame(pvals)
pvals$Gene <- rownames(pvals)
tmp <- merge(tmp, pvals, by="Gene")

tmp$sig <- sigSymbol(tmp$pval)
tmp$Gene2 <- apply(tmp[c("Gene", "sig")], 1, function(x) paste0(x, collapse=""))
p <- ggplot(aes(mrnaz, reorder(Gene2, mrnaz, mean)), data = subset(tmp, PIT&!is.na(mrnaz)))
p <- p  + geom_jitter( aes(color = Gene2, alpha=1, size=2)) + facet_grid(cat~., drop = T, scales="free_y", space="free_y")
p <- p + xlab("Tumor mRNA expression") + ylab("") + geom_vline(xintercept = 0, color="red")
p <- p + scale_color_manual(values=colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tmp$Gene))))
p <- p + theme( panel.background = element_rect(fill='white',colour='black'), axis.text.y= element_text(angle=0), legend.position="none")
p <- p + theme(panel.grid.major.y = element_line(linetype=3, color="darkgray"))
p <- p + theme( axis.text = element_text(size= rel(1.5))) 
p <- p + xlim(c(-4.1, 4))
p

