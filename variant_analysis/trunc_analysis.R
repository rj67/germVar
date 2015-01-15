########## 
## mosaic plot of CGC500 variants in TS vs ONCO
tmp<-list_goi
tmp$cat <- apply(tmp[c("TS", "ONCO")], 1, function(x) ifelse(x[1], "TS", ifelse(x[2], "ONCO", "OTHER")))
tmp2 <- merge(VAR, tmp[c("Gene", "cgc500", "cat")], by="Gene")
tmp3 <- rbind(table(subset(tmp, cgc500)$cat), table(subset(tmp2, cgc500)$cat))
rownames(tmp3) <- c("Gene", "Variant")
mosaicplot(tmp3, main="", color=c("pink", "lightgreen", "lightblue"), cex.axis=1.5)

########## 
## boxplot of top genes in CGC500 with truncating variants
tmp <- subset(mut_stat, Gene %in% subset(list_goi, cgc500)$Gene)
tmp <- merge(tmp, list_goi[c("Gene", "TS", "ONCO")], by="Gene")
tmp$cat <- apply(tmp[c("TS", "ONCO")], 1, function(x) ifelse(x[1], "TS", ifelse(x[2], "ONCO", "OTHER")))
tmp2 <- arrange(dplyr::summarise(group_by(tmp, Gene), mut_count = length(unique(event_uid[PIT]))), -mut_count)
tmp <- subset(merge(tmp, tmp2, by="Gene"), mut_count >=10)
tmp <- tmp[!duplicated(tmp$mut_uid),]
#tmp <- merge(tmp, CBio_query[c("event_uid", "log2CNA",  "gistic", "mrnaz",	"study", "gistic2")], by="event_uid", all.x=T)
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
tmp$sig <- sapply(tmp$pvals, function(x) ifelse(x>0.05, "", ifelse(x>0.01, "*", ifelse(x>0.001, "**", "***"))))
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



###############################
## boxplot of Patient's age harboring CGC500 TS truncating variants
tmp <- subset(mut_stat, Gene %in% subset(list_goi, cgc500)$Gene)
tmp <- merge(tmp, list_goi[c("Gene", "TS", "ONCO")], by="Gene")
tmp$cat <- apply(tmp[c("TS", "ONCO")], 1, function(x) ifelse(x[1], "TS", ifelse(x[2], "ONCO", "OTHER")))
tmp2 <- arrange(dplyr::summarise(group_by(tmp, Gene), mut_count = length(unique(event_uid[PIT]))), -mut_count)
tmp <- subset(merge(tmp, tmp2, by="Gene"), PIT&TS)
tmp <- tmp[!duplicated(tmp$mut_uid), ]
tmp <- merge(tmp, all_tcga[c("Patient", "disease")])
tmp$study <- tmp$disease
pvals <- sapply( unique(tmp$study),  function(x){
  print(x);
  to_test <- subset(tmp, study==x);
  bg <- subset(all_clin[c("Patient", "study", "age")], study==x)
  return(wilcoxTestCallSet(to_test, "age", bg))
})
pvals <- as.data.frame(pvals)
pvals$study <- rownames(pvals)
tmp <- merge(tmp, pvals, by="study")
tmp$sig <- sapply(tmp$pvals, function(x) ifelse(x>0.05, "", ifelse(x>0.01, "*", ifelse(x>0.001, "**", "***"))))
tmp$study2 <- apply(tmp[c("study", "sig")], 1, function(x) paste0(x, collapse=""))
p <- ggplot(aes(reorder(study, age, median), age), data = all_clin) + geom_boxplot(outlier.shape = NA)
p <- p + geom_jitter( aes(color = study, alpha=1, size=2), data= subset(all_clin, Patient %in% tmp$Patient))
p <- p + theme( panel.background = element_rect(fill='white',colour='black'), axis.text.y= element_text(angle=0), legend.position="none")
p <- p + theme(panel.grid.major.y = element_line(linetype=3, color="darkgray"))
p <- p + theme( axis.text.x = element_text(size= rel(1.5), angle=45, hjust=1, vjust=1)) 
p <- p + scale_color_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(length(unique(tmp$study))))
p <- p + xlab("") + ylab("")
p


############################
## bubble plot of exp.pval vs X2kG ratio
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
library(scales)
tmp <- merge(mut_stat, list_goi, by="Gene")
Gene_stat <- arrange(dplyr::summarise(group_by(tmp, Gene), mut_count = length(unique(event_uid[PIT])), loh_count = length(unique(event_uid[LOH]))), -mut_count)
Gene_stat <- merge(Gene_stat, dplyr::summarise(group_by(X2kG_trunc, Gene), CNTR_AC =sum(X2kG_AC)), by="Gene")
Gene_stat <- subset(Gene_stat, mut_count>=10)
Gene_stat$AN <- 15086
Gene_stat$X2kG_AN <- 5008
Gene_stat <- merge(Gene_stat, list_goi)
tmp <- subset(tmp, Gene %in% Gene_stat$Gene & PIT)
tmp <- tmp[!duplicated(tmp$mut_uid),]
tmp <- subset(tmp, Gene %in% RNASeq$Gene)

exp.pvals <- sapply( unique(tmp$Gene),  function(x){
  print(x);
  to_test <- subset(tmp, Gene==x);
  bg <- subset(RNASeq, Gene==x)
  return(wilcoxTestCallSet(to_test, "mrnaz", bg))
})
exp.pvals <- as.data.frame(exp.pvals)
exp.pvals$Gene <- rownames(exp.pvals)
Gene_stat <- merge(Gene_stat, exp.pvals, by="Gene")
Gene_stat <- cbind(Gene_stat, t(apply(Gene_stat[c("mut_count", "AN", "CNTR_AC", "X2kG_AN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "X2kG", "greater"))})))
Gene_stat$X2kG.pval[Gene_stat$X2kG.pval<10^(-6)] <- 10^(-6)
Gene_stat$exp.pvals[Gene_stat$exp.pvals<10^(-6)] <- 10^(-6)
Gene_stat$ratio <- log(Gene_stat$mut_count/Gene_stat$CNTR_AC*5008/15086, 2)
#p <- ggplot(Gene_stat, aes(exp.pvals, X2kG.pval, label=Gene)) + geom_jitter(aes(color="green", alpha=0.8, size=sqrt(mut_count))) + scale_size(range = c(3, 10))
p <- ggplot(Gene_stat, aes(exp.pvals, ratio, label=Gene)) + geom_jitter(aes( color="red", alpha=0.9, size=sqrt(mut_count))) + scale_size( range = c(3, 13))
#p <- p + scale_x_continuous(trans=reverselog_trans(10)) + scale_y_continuous(trans=reverselog_trans(10))
p <- p + scale_x_continuous(trans=reverselog_trans(10)) + scale_y_continuous(limits = c(-2.75, 4.8))
p <- p + geom_vline(xintercept = 0.05, color="red")
p <- p + theme( panel.background = element_rect(fill='white',colour='black'), axis.text.y= element_text(angle=0), legend.position="none")
p <- p + theme(panel.grid.major.y = element_line(linetype=3, color="darkgray"))
p <- p + theme( axis.text = element_text(size= rel(1.5))) 
p <- p + geom_text(data=subset(Gene_stat, X2kG.pval<0.01 ), color="blue", position=position_jitter(height=0.4, width=0.1))
p <- p + geom_text(data=subset(Gene_stat, X2kG.pval>0.01 & exp.pvals < 0.01 ), color="black", position=position_jitter(height=0.4, width=0.2))
p <- p + xlab("") + ylab("") 
p
#library(wordcloud)
#Gene_stat$x <- log10(Gene_stat$exp.pvals)
#Gene_stat$y <- log10(Gene_stat$X2kG.pval)
#mx <- apply(Gene_stat[c("x", "y")],2,max)
#mn <- apply(Gene_stat[c("x", "y")],2,min)
#Gene_stat2 <- subset(Gene_stat, x< -1.5 | y< -1.5)
#textplot(Gene_stat2$x, Gene_stat2$y, Gene_stat2$Gene, xlim=c(mn[1],mx[1]),ylim=c(mn[2],mx[2]))



############################
## beanplot mut_count per gene in CGC500
tmp <- subset(merge(mut_stat, list_goi, by="Gene"), cgc500)
tmp <- tmp[!duplicated(tmp$mut_uid), ] 
Gene_stat <- arrange(dplyr::summarise(group_by(tmp, Gene), mut_count = length(unique(event_uid[PIT]))), -mut_count)
Gene_stat <- merge(Gene_stat, list_goi[c("Gene", "TS", "ONCO")], by="Gene")
Gene_stat$cat <- apply(Gene_stat[c("TS", "ONCO")], 1, function(x) ifelse(x[1], "TS", ifelse(x[2], "ONCO", "OTHER")))
require(beanplot)
beanplot(mut_count ~ cat, data=Gene_stat, bw="SJ", log="", ll=0.02, col=c("purple", "lightblue", "black"), what=c(0,1,0,1), xlab="Gene Type",
         main="Allele frequence of CGC truncating variants")



############################
## 
tmp <- subset(merge(mut_stat, list_goi, by="Gene"), cgc500 & TS & PIT)
tmp <- tmp[!duplicated(tmp$mut_uid), ] 
tmp <- merge(tmp, all_tcga[c("Patient", "disease")])
tally_gene_disease <- dplyr::summarise(group_by(tmp, Gene, disease), num_sm = length(unique(Patient))  )
tally_disease <- dplyr::summarise(group_by(tmp, disease), num_sm = length(unique(Patient))  )
tally_disease <- merge(tally_disease, as.data.frame(table(all_tcga[!duplicated(all_tcga$Patient),]$disease)), by.x="disease", by.y="Var1")
tally_disease$ratio <- with(tally_disease, num_sm/Freq)
tally_disease <- arrange(tally_disease, -ratio)
tally_rep_disease <- dplyr::summarise(group_by(subset(tmp, rep), rep_id, disease), num_sm = length(unique(Patient))  )

tally_rep_disease$disease <- factor(tally_rep_disease$disease, levels=tally_disease$disease)
# plot Gene X disease
#p1<-ggplot(data=tally_gene_disease, aes(reorder(disease, -num_sm, sum), reorder(Gene, -num_sm, sum), size= sqrt(num_sm),  label=num_sm)) 
p1<-ggplot(data=tally_rep_disease, aes(disease, reorder(rep_id, -num_sm, sum), size= sqrt(num_sm),  label=num_sm)) 
p1 <- p1 + geom_point(shape=1, color="grey50", alpha=0.5) + geom_point(alpha=0.9, color="pink") 
p1 <- p1 + geom_text(color="black", size=5, alpha=0.9)
p1 <- p1 + theme( axis.text.x = element_text(size= rel(1.5), angle = 45, vjust = 0.9, hjust = 0.9)) 
p1 <- p1 + theme( panel.background = element_rect(fill='white',colour='black'), axis.text.y= element_text(angle=0), legend.position="none")
p1 <- p1 + theme(panel.grid.major.y = element_line(linetype=3, color="darkgray"))
p1 <- p1 + theme(legend.position="none")
p1 <- p1 + theme(panel.grid.major.y = element_line(linetype=3, color="darkgray"))
p1 <- p1 + scale_size_continuous(range = c(2,12), guide = "none") 
#p1 <- p1 + scale_colour_manual(values = rev(brewer.pal(9,"PuRd")))
#p1 <- p1 + scale_colour_manual(values = c("#49006a", "#ae017e", "#f768a1", "#fcc5c0", "#fff7f3"))
p1 <- p1 + theme( axis.title = element_blank())
#p1 <- p1 + ylab("Gene") + xlab("Disease")
p1


p2 <- ggplot(tally_disease, aes(reorder(disease, -ratio), ratio)) + geom_bar(stat="identity", fill="red", alpha=0.7)
p2 <- p2 + theme( axis.text.y = element_text(size= rel(1.5)))
p2 <- p2 + theme( axis.text.x = element_blank())
p2 <- p2 + theme( panel.background = element_rect(fill='white',colour='black'), axis.text.y= element_text(angle=0), legend.position="none")
p2 <- p2 + theme(panel.grid.major.y = element_line(linetype=3, color="darkgray"))
p2 <- p2 + xlab("") + ylab("")
p2
