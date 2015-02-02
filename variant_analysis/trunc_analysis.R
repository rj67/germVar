## remove complex variants, which results in no net LoF, because insertion and deletion cancel out, 
## "XRCC3-104174968-TA-T, XRCC3-104174895-C-CCT, XRCC3-104174897-G-GGGAC suspect
#"TNFSF12-7460169-ATCGG-A", "TNFSF12-7460164-GC-G", OK
# "ARID3B-74836316-A-AG, ARID3B-74836315-A-AGC, suspect
# BAD-64037722-C-CT, BAD-64037724-C-T, OK
# BCLAF1 OK
# BRCA2 EI-7004 OK
# DSC3 G3-A7M8 OK
# EP400 E9-A3Q9 EP400-132547138-A-AACAG suspect, compound with EP400-132547137-A-AGC?
# OGG1 DW-7839 OK
# RECQL5 EW-A6SD OK
# SERPINI2, SERPINI2-167183312-TTG-T, SERPINI2-167183316-ACC-A, SERPINI2-167183321-TTTTC-T OK
# SPO11 DJ-A13T OK

complex_list <- c("PEG3-57335976-C-T", "PEG3-57335978-AGGAC-A", "PEG3-57335984-T-A", "PEG3-57335985-C-CAGAT",
                  "TIMELESS-56817451-T-TCC", "TIMELESS-56817449-T-TC",
                  "TSC2-2132470-GCCAAGGA-G", "TSC2-2132467-GCT-G",
                  "CCND3-41904348-CA-C", "CCND3-41904346-GA-G", "CCND3-41904350-TGGAGCAG-T",
                  "CREBBP-3778439-T-TGC", "CREBBP-3778438-G-GC",
                  "KLK2-51379788-G-GAA", "KLK2-51379791-C-CACAG",
                  "NCOR1-16089934-T-TAACA", "NCOR1-16089930-CTGTT-C",
                  "NEURL1-105331516-ACGGCCGACCCGCTC-A", "NEURL1-105331531-T-TA", 
                  "NIT2-100067704-G-A", "NIT2-100067705-G-A",
                  "POU6F2-39379279-C-CA", "POU6F2-39379281-C-CAG")
trunc_VARu <- subset(trunc_VARu, !var_uid %in% complex_list)

#####################################################
# Mutational data
#####################################################
trunc_MUTs <-  joint_tally$MUT %>% merge(., as.data.frame(info(trunc_VARu))[c("var_uid", "Gene", "tier1")], by="var_uid") %>%
  merge(., all_tcga[!duplicated(all_tcga$Patient),][c("Patient", "age", "disease", "disease2")], by="Patient") %>%
  merge(., list_goi[c("Gene", "cgc500", "Cat", "HER", "rep", "rep_id")]) %>% subset(., tier1 & PIT)

trunc_gene <- group_by(trunc_MUTs, Gene) %>% dplyr::summarise( N_pat = length(unique(Patient)))


#####################################################
## lf frequency comparision
#####################################################
trunc_lf <- as.data.frame(subset(info(trunc_VARu), EAC>=10 & tier1))
# examine the variants
View(subset(trunc_lf, tier1 & EAC>=10 & Gene %in% subset(list_goi, cgc500)$Gene)[c("Gene", "EFF", "AALength", "AAChange", "tier1", "tier2", "EAC", "CAC1","NAC", "TAC", "X2kG_AC", "ESP_AC", "var_uid")])

# EA1 with ESP_EA
trunc_lf <- cbind(trunc_lf,  do.call(rbind, apply(trunc_lf[c("CAC1", "CAN1", "ESP_EA_AC", "ESP_EA_AN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "EA1" ))})))
trunc_lf$EA1.t.padj <- signif(p.adjust(trunc_lf$EA1.t.pval, method="BH"), 2)

# compare with X2kG
trunc_lf$X2kG_AN <- 5008
trunc_lf <- cbind(trunc_lf, do.call(rbind, apply(trunc_lf[c("EAC", "EAN", "X2kG_AC", "X2kG_AN")], 1, function(x) {names(x) <- NULL;do.call(varBurden, c(as.list(x), "X2kG"))})))
trunc_lf$X2kG.t.padj <- signif(p.adjust(trunc_lf$X2kG.t.pval, method="BH"), 2)

View(subset(trunc_lf, EA1.t.padj <0.2 & X2kG.t.padj<0.15  )[c("STR_match","EAC", "EAN", "EAF","NAC","TAC", "CAC1", "CAN1","ESP_AC", "ESP_AN", "ESP_EA_AC", "ESP_EA_AN", "X2kG_AF", "Gene", "uid", "AAChange.p", 
                                                                "EA1.f.pval", "EA1.t.padj","X2kG.f.pval","X2kG.t.padj" )])

lf_table <- subset(trunc_VARu_lf, (EA1.padj <0.2 & X2kG.padj<0.15)| (Gene %in% subset(list_goi, cgc500 & Cat=="TSG")$Gene))[c( "EAF", "CAC1", "CAN1", "ESP_EA_AC", "ESP_EA_AN", "X2kG_AF", "Gene", "var_uid", "AAChange.p", "EA1.padj", "X2kG.padj", "NAC", "TAC")]
lf_table <- lf_table %>% mutate( EAF = signif(EAF, 2), CAF1 = signif(CAC1/CAN1, 2), ESP_EA_AF = signif(ESP_EA_AC/ESP_EA_AN, 2), X2kG_AF = signif(X2kG_AF, 2), EA1.padj = signif(EA1.padj, 2), X2kG.padj = signif(X2kG.padj, 2))
lf_table <- arrange(lf_table[c("Gene", "AAChange.p", "EAF", "CAF1", "ESP_EA_AF", "X2kG_AF", "EA1.padj", "X2kG.padj", "var_uid")], EA1.padj)


#####################################################
## test age of onset for each gene
#####################################################
# test low frequency variant and age
trunc_lf_age <- group_by(subset(trunc_MUTs, var_uid %in% trunc_lf$var_uid), var_uid) %>% do(age_test(.))

# summ all variants in each Gene, test age of onset
trunc_age <- group_by(droplevels(subset(trunc_MUTs, Gene %in% subset(trunc_gene, N_pat>=5)$Gene)), Gene) %>% do(age_test(.))
trunc_age <- subset(trunc_age, N_var>=8)
trunc_age$padj <- p.adjust(trunc_age$pval, method="BH") 
View(arrange(trunc_age, pval))
library(scales)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
p <- ggplot(trunc_age, aes(effect, pval, label=Gene)) + geom_point(alpha=0.5, aes(size= sqrt(N_var))) +  scale_size_continuous(range = c(2, 8), guide="none")
p <- p + geom_text(data = subset(trunc_age, pval <=0.05), position=position_jitter(height=0.2, width=0.05), size=2)
p <- p + theme_few() + scale_y_continuous(trans=reverselog_trans(10))
p <- p + geom_hline(yintercept=0.0035, color="red") + geom_text(x=0.7, y = 2.6, label="FDR 10%", hjust=0.5, color="red", size=3)
p
ggsave(filename="trunc_age_test.png",width=5, height=4)

age_test <- function(df, bg=all_tcga){
  # get all patient age
  age_df <- subset(bg, !duplicated(Patient))[c("Patient", "disease", "agez")] %>% subset(., !is.na(agez) & is.finite(agez))
  age_df$VAR <- age_df$Patient %in% df$Patient
  N_all <- nrow(age_df)
  N_var <- sum(age_df$VAR)
  if (N_var==0) {
    warning("Not enough VAR patient")
    return(data.frame(N_all = N_all, N_var = 0, effect = 0, std = 0, pval = NA))  
  } else{
    # fit linear model  
    lm.mod <- lm(agez ~  VAR, age_df) 
    # get average age diff
    effect <- broom::tidy(lm.mod)$estimate[broom::tidy(lm.mod)$term=="VARTRUE"]
    std <- broom::tidy(lm.mod)$std.error[broom::tidy(lm.mod)$term=="VARTRUE"]
    # calc anova pval
    #pval <- summary.aov(lm.mod)[[1]][2,5]
    pval <- summary.aov(lm.mod)[[1]][1,5]
    return(data.frame(N_all = N_all, N_var = N_var, effect = effect, std = std, pval = pval))  
  }
}

#####################################################
## Test the number of variants per gene or per group of genes
#####################################################

trunc_bygene <- merge(as.data.frame(table(subset(info(trunc_VARu), tier1)$Gene)), as.data.frame(table(subset(info(lof_X2kG), tier1)$Gene)), by="Var1", all=T) %>% 
  arrange(., -Freq.x) %>%  mutate(., Freq.x=replaZero(Freq.x), Freq.y=replaZero(Freq.y)) %>% plyr::rename(., c("Var1"="Gene", "Freq.x"="TCGA", "Freq.y"="X2kG")) %>% 
  plyr::join(., list_goi[c("Gene", "rep", "rep_id", "fam", "fam_id", "grp", "grp_id", "HER", "Cat", "mem")], by="Gene")

trunc_bygene %>% subset(., !Gene %in% c("BRCA1", "BRCA2")) %>% subset(., mem=="CGC"|HER) %>% group_by(., Cat) %>% do(var_ftest(.))
trunc_bygene %>% subset(., !Gene %in% c("BRCA1", "BRCA2")) %>% subset(., mem=="CGC"|HER) %>% subset(., HER & !rep) %>% var_ftest
trunc_bygene %>% subset(., !Gene %in% c("BRCA1", "BRCA2"))  %>% subset(., rep) %>% var_ftest


rep_test <- trunc_bygene %>% subset(., !Gene %in% c("BRCA1", "BRCA2")) %>% subset(., rep &  (mem=="CGC"| HER )) %>% droplevels(.) %>% group_by(., rep_id)  %>% do(var_ftest(.))
rep_test$padj <- sigSymbol(rep_test$pval)
rep_test <- arrange(rep_test, -NV_TCGA)
rep_test$rep_id <- factor(rep_test$rep_id, levels=rep_test$rep_id)
rep_test$sig <- sigSymbol(rep_test$pval)

fam_test <- trunc_bygene %>% subset(., !Gene %in% c("BRCA1", "BRCA2")) %>% subset(., fam) %>% droplevels(.) %>% group_by(., fam_id)  %>% do(var_ftest(.))
grp_test <- trunc_bygene %>% subset(., !Gene %in% c("BRCA1", "BRCA2")) %>% subset(., grp) %>% droplevels(.) %>% group_by(., grp_id)  %>% do(var_ftest(.))

# test the number of variants in arbitrary group of genes, compared to X2kG
var_ftest <-function(df){
  NV_TCGA <- sum(df$TCGA)  
  NV_X2kG <- sum(df$X2kG)  
  ### also perform fisher exact test
  data <- matrix(c(NV_TCGA, NV_X2kG, 2524-NV_TCGA, 941-NV_X2kG), nrow=2)
  ftest <- fisher.test(data, alternative="greater")
  returns<- data.frame(NV_TCGA =NV_TCGA, NV_X2kG = NV_X2kG, pval = ftest$p.value, odds = ftest$estimate, conf.lo = ftest$conf.int[1], conf.hi = ftest$conf.int[2])
  returns[,3:6] <- sapply(returns[,3:6], function(x) signif(x, 2))
  return(returns)
}
to_plot <- plyr::join(reshape2::melt(rep_test[c("rep_id", "NV_TCGA", "NV_X2kG")], id.vars = "rep_id", value.name="Freq", variable.name="Source"), rep_test[c("rep_id", "sig")], by="rep_id")
p1 <- ggplot(aes(rep_id, Freq, fill=Source, label=sig), data= to_plot)
p1 <- p1 + theme_few() + geom_bar(stat="identity", position="dodge")
p1 <- p1 +  theme(axis.text.x= element_text(vjust=1, hjust=1, size=rel(0.7), angle=45), legend.position = c(.85, .8)) + scale_fill_manual(values=c("grey10", "grey80"), guide = guide_legend(title = ""), labels=c("TCGA", "1000G"))
p1 <- p1 + xlab("") + ylab("Variant count")  + theme(axis.title.y=element_text(size=rel(0.7)))
p1 <- p1 + geom_text(data=subset(to_plot, Source=="NV_TCGA"), vjust=0.5)
p1 <- p1 + scale_x_discrete(labels = c( "FA"="Fanconi\nanemia", "CP"="Damage\nresponse", "MISC"="Miscelaneous", "NER"="Nucleotide\nexcision repair", "MMR"="Mismatch\nexcision repair", "HR"="Homologous\nrecombination","BER"="Base\nexcision repair"))
p1
ggsave(filename="trunc_var_rep_hist.png",height=4, width=4)

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
trunc_mrnaz <- subset(trunc_mrnaz, N_var>=3)
trunc_mrnaz$padj <- p.adjust(trunc_mrnaz$pval, method="BH") 

View(arrange(trunc_mrnaz, pval))
#cgc_test <- subset(trunc_mrnaz_test, Gene %in% subset(list_goi, cgc500)$Gene & N_var >=5)

p <- ggplot(trunc_mrnaz, aes(effect, pval, label=Gene)) + geom_point(alpha=0.5, aes(size= sqrt(N_var)), fill="lightgrey") +  scale_size_continuous(range = c(2, 8), guide="none")
p <- p + scale_x_continuous(limits=c(-2.3, 2.3))
p <- p + geom_text(data = subset(trunc_mrnaz, pval <=0.001), position=position_jitter(height=0.25, width=0.15), size=2, color="black", alpha=0.9)
p <- p + theme_few() + scale_y_continuous(trans=reverselog_trans(10))
p <- p + geom_hline(yintercept=0.027, color="red", linetype=3) + geom_text(x=2.4, y = 1.8, label="FDR 10%", hjust=1, color="red", size=3)
p <- p + geom_hline(yintercept=0.0012, color="red", linetype=3) + geom_text(x=2.4, y = 3.1, label="FDR 1%", hjust=1, color="red", size=3)
p <- p + theme(axis.text.y=element_text(angle=90, hjust=0.5), axis.title=element_text(size=rel(0.7)))
p <- p + ylab("P-value") + xlab("Average mRNA level change")
p
ggsave(filename="trunc_mrnaz_test.png",width=4, height=7)


mrnaz_test <- function(df){
  bg <- getRNASeq(unique(df$Gene))  
  # remove NAs and med==0 tissues
  bg <- subset(bg, !is.na(mrnaz) & med!=0)
  N_all <- nrow(bg)
  bg$VAR <- bg$Patient %in% df$Patient
  N_var <- sum(bg$VAR)
  if( N_var > 0){
    lm.mod <- lm( mrnaz ~  VAR, bg) 
    # get average age diff
    effect <- broom::tidy(lm.mod)$estimate[broom::tidy(lm.mod)$term=="VARTRUE"]
    std <- broom::tidy(lm.mod)$std.error[broom::tidy(lm.mod)$term=="VARTRUE"]
    test <- wilcox.test(subset(bg, VAR)$mrnaz, subset(bg, !VAR)$mrnaz)#, alternative = "less")
    return( data.frame(N_all = N_all, N_var = N_var, effect = effect, std = std, pval=test$p.value))
  }else{
    return( data.frame(N_all = N_all, N_var = 0, effect = 0, std = 0, pval = NA))
  }
}


plot_RNASeq <- function(df){
  library(ggthemes)
  library(RColorBrewer)
  bg <- getRNASeq(unique(df$Gene))
  Patients <- df$Patient
  # drop diseases with no variant
  bg <- droplevels(subset(bg, study %in% droplevels(subset(bg, Patient %in% Patients))$study))
  p <- ggplot(aes(reorder(study, normalized_count, median), sqrt(normalized_count)), data = bg) + geom_boxplot(outlier.shape = NA)
  p <- p + geom_jitter(aes(color=study), data = subset(bg, Patient %in% Patients))
  p <- p + theme_few()
  p <- p + scale_color_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(length(unique(bg$study))), guide="none")
  p <- p + xlab("study") + ylab("sqrt(Normalized count)")
  p 
}

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


#####################################################

TS_muts$rep_id[is.na(TS_muts$rep_id)] <- "non_rep" 

gene_disease_tally <- TS_muts %>% group_by(., rep_id, Gene, disease) %>% dplyr::summarise(., N_pat = length(unique(Patient[PIT])))
#gene_disease_tally$Gene <- factor(gene_disease_tally$Gene, levels= arrange(droplevels(subset(gene_disease_tally, !duplicated(Gene))), rep_id)$Gene)
gene_disease_tally <- gene_disease_tally %>% merge(., subset(all_tcga, !duplicated(Patient)) %>% group_by(disease) %>% dplyr::summarise(study_size = length(unique(Patient))), by="disease") %>%
                      mutate(., P_pat = N_pat/study_size)
rep_disease_tally <- TS_muts  %>% group_by(., rep_id, disease) %>% dplyr::summarise(., N_pat = length(unique(Patient[PIT])))
gene_disease_tally$disease <- factor(gene_disease_tally$disease, hclust_factor(rep_disease_tally, "disease", "rep_id", "N_pat", T))
gene_disease_tally$Gene <- factor(gene_disease_tally$Gene, hclust_factor(gene_disease_tally, "Gene", "disease", "P_pat", T))

#gene_disease_tally$disease <- factor(gene_disease_tally$disease, levels= hclust_factor(gene_disease_tally, "disease", "Gene", "N_pat", T))

hclust_factor <- function(df, Var1, Var2, value, output=F){
  library(reshape2)
  library(gclus)
  df <- df[c(eval(Var1), eval(Var2), eval(value))]
  colnames(df) <- c("Var1", "Var2", "value")
  mat <- apply(acast(df, Var1 ~ Var2, value.var= "value" ), 2, replaZero)
  mat.clust <- hclust(dist(mat), method="ward.D")
  mat.reorder <- reorder(mat.clust, dist(mat))
  if(output){ plot(mat.reorder)}
  return(rownames(mat)[mat.reorder$order])
}

#gene_disease_tally$rep_id <- factor(gene_disease_tally$rep_id)
p2 <- ggplot(aes(disease, Gene), data=gene_disease_tally) + geom_point(aes(size=sqrt(N_pat)))
p2 <- p2 + theme( panel.background = element_rect(fill='white',colour='black'))
p2 <- p2 + facet_grid(rep_id~., drop=TRUE, scales="free_y", space = "free_y")
p2 <- p2 + theme( panel.grid.major.y = element_line(linetype=3, color="darkgray"),  panel.grid.major.x = element_line(linetype=1, color="lightgray"))
p2 <- p2 + theme(axis.text.x= element_text(hjust=1., vjust = 0.5, angle=90)) + xlab("") +ylab("")
p2

grid.arrange(p1, p2, ncol=1)


#####################################################
# test patient age, see if 
MUTS <- joint_tally$MUT %>% merge(., as.data.frame(info(trunc_VARu))[c("var_uid", "Gene", "tier1")], by="var_uid") %>%
  merge(., all_tcga[!duplicated(all_tcga$Patient),][c("Patient", "age", "disease")], by="Patient") %>%
  merge(., list_goi[c("Gene", "cgc500", "Cat", "HER", "rep", "rep_id")]) 
  
age_df <- subset(all_tcga, !duplicated(Patient))[c("Patient", "disease", "age")]
age_df$TS <- age_df$Patient %in% subset(MUTS, PIT & Gene=="LZTR1" )$Patient
summary.aov(lm(age~disease+TS, age_df))

#summary.aov(lm(age~disease+TS, subset(age_df, !Patient %in% subset(TS_muts, Gene %in% c("BRCA1", "BRCA2", "CHEK2", "ERCC3"))$Patient)))


########## 
## percentage of LOH by gene
MUTs <- merge(joint_tally$MUT, as.data.frame(info(trunc_VARu))[c("var_uid", "Gene")], by="var_uid")

LOH_tally <- MUTs %>% group_by(., Gene) %>% dplyr::summarise(., N_mut = length(unique(mut_uid)), per_LOH = length(unique(mut_uid[LOH]))/length(unique(mut_uid)))
View(arrange(LOH_tall, -N_mut))


########## 
## variant per gene by Cat
tmp<-merge(list_goi[c("Gene", "cgc500", "Cat")], as.data.frame(table(subset(info(trunc_VARu), tier1)$Gene)), by.x="Gene", by.y="Var1", all.x=T)
tmp$Freq[is.na(tmp$Freq)] <- 0


########## 
## boxplot of top genes in CGC500 with truncating variants


gene_tally$pvals <- pvals




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
