if (F){
  load("./DataSet_Results/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.goi.RData")
  call_set <- droplevels(subset(call_set, Gene %in% list_goi$Gene & Impact %in% c("MODERATE", "HIGH")))
  call_set <- call_set %>% labelVarUid %>% labelUid
  # remove duplicate entries
  call_set <- arrange(call_set, var_uid, Impact)
  call_set <- call_set[!duplicated(call_set$var_uid), ]
  call_set <- call_set %>% plyr::rename(., replace=c( "AC"="X2kG_AC", "AF"="X2kG_AF", "AN"="X2kG_AN"))
  X2kG_goi <- call_set
  #save(X2kG_goi, file="Results/X2kG_goi.RData")
  
  #load("./DataSet_Results/ALL.wgs.phase1_release_v3.20101123.snps_indels.sites.goi.RData")
  #call_set <- droplevels(subset(call_set, Gene %in% list_goi$Gene & Impact %in% c("MODERATE", "HIGH")))
  #call_set <- createVarUid(call_set)
  # remove duplicate entries
  #call_set <- arrange(call_set, var_uid, Impact)
  #call_set <- call_set[!duplicated(call_set$var_uid), ]
  #colnames(call_set)[colnames(call_set)=="X1kG_AF"] <- "X1kG_LDAF"
  #colnames(call_set)[colnames(call_set)=="AF"] <- "X1kG_AF"
  #colnames(call_set)[colnames(call_set)=="ERATE"] <- "X1kG_ERATE"
  #X1kG_goi <- call_set
  #save(X1kG_goi, file="Results/X1kG_goi.RData")
  
  load("./DataSet_Results/ESP6500SI-V2-SSA137.updatedProteinHgvs.snps_indels.goi.RData")
  call_set <- droplevels(subset(call_set, Gene %in% list_goi$Gene & Impact %in% c("MODERATE", "HIGH")))
  call_set <- call_set %>% labelVarUid %>% labelUid
  # remove duplicate entries
  call_set <- arrange(call_set, var_uid, Impact)
  call_set <- call_set[!duplicated(call_set$var_uid), ]
  ESP_goi <- call_set
  # hack ESP to mimic TCGA ethnic composition
  #tmp<-apply(ESP_goi[c("ESP_AA_AC", "ESP_AA_AN", "ESP_EA_AC", "ESP_EA_AN")], 1, function(x) return(c(ESP_fAC=floor(x[1]/9*x[4]/x[2]+x[3]), ESP_fAN=floor(x[4]*(1+1/9)))))
  #tmp<- as.data.frame(t(tmp))
  #colnames(tmp) <- c("ESP_fAC", "ESP_fAN")
  #ESP_goi <- cbind(ESP_goi, tmp)
  save(X2kG_goi, ESP_goi, file="DataSet_Results/ESP_X2kG_goi.RData")
}


### ### ### ### ### ### 
### check with CCDS
### ### ### ### ### ### 
ESP_trunc <- subset(ESP_goi, Impact=="HIGH")
ESP_trunc$CDS <- apply(ESP_trunc[c("Gene", "POS", "REF", "ALT")], 1, function(x) do.call(queryCCDS, as.list(x)))
ESP_trunc <- subset(ESP_trunc, CDS !=0)

X1kG_trunc <- subset(X1kG_goi, Impact=="HIGH" & EFF!="STOP_LOST")
X1kG_trunc$CDS <- apply(X1kG_trunc[c("Gene", "POS", "REF", "ALT")], 1, function(x) do.call(queryCCDS, as.list(x)))
X1kG_trunc <- subset(X1kG_trunc, CDS !=0)
X1kG_trunc <- subset(X1kG_trunc, (is.na(AA_pos) | abs(AA_pos-AALength)/AALength > 0.05 | (Pfam==1 & AftDom==0)))

X2kG_trunc <- subset(X2kG_goi, Impact=="HIGH" & EFF!="STOP_LOST")
X2kG_trunc$CDS <- apply(X2kG_trunc[c("Gene", "POS", "REF", "ALT")], 1, function(x) do.call(queryCCDS, as.list(x)))
X2kG_trunc <- subset(X2kG_trunc, CDS !=0)
X2kG_trunc <- subset(X2kG_trunc, (is.na(AA_pos) | abs(AA_pos-AALength)/AALength > 0.05 | (Pfam==1 & AftDom==0)))
X2kG_trunc <- subset(X2kG_trunc, X2kG_AF < 0.99)
  
### ### ### ### ### ### 
### check variants in X2kG more than 10 times(AF 0.002) but not in X1kG 
### ### ### ### ### ### 
novel_X2kG <- merge(X2kG_trunc[c("var_uid", "X2kG_AF", "VARTYPE")], X1kG_trunc[c("var_uid", "X1kG_AF", "VARTYPE")], by=c("var_uid", "VARTYPE"), all=T)
novel_X2kG <- subset(novel_X2kG, is.na(X1kG_AF) & X2kG_AF > 0.002)
# compare with ESP
novel_X2kG <- arrange(merge(novel_X2kG, ESP_trunc[c("var_uid", "ESP_AF")], by="var_uid", all.x=T), -X2kG_AF)
# compare with truncation calls
novel_X2kG <- arrange(merge(novel_X2kG, hardFilter(trunc_calls, T)@VAR[c("var_uid", "AF", "AC", "AN")], by="var_uid", all.x=T), -X2kG_AF)
# flag variants not in X1kG, ESP or trunc_calls, as suspicious.
novel_X2kG$Suspect <- with(novel_X2kG, is.na(ESP_AF)& is.na(AF))

### ### ### ### ### ### 
### check variants in X1kG more than 10 times(AF 0.002) but not in X2kG 
### ### ### ### ### ### 
miss_X2kG <- merge(X2kG_trunc[c("var_uid", "X2kG_AF", "VARTYPE")], X1kG_trunc[c("var_uid", "X1kG_AF", "X1kG_ERATE","VARTYPE")], by=c("var_uid", "VARTYPE"), all=T)
miss_X2kG <- subset(miss_X2kG, is.na(X2kG_AF) & X1kG_AF > 0.002)
miss_X2kG <- arrange(merge(miss_X2kG, ESP_trunc[c("var_uid", "ESP_AF")], by="var_uid", all.x=T), -X1kG_AF)
miss_X2kG <- arrange(merge(miss_X2kG, hardFilter(trunc_calls, T)@VAR[c("var_uid", "AF", "AC", "AN")], by="var_uid", all.x=T), -X1kG_AF)

both_high<-merge(subset(X2kG_goi, Impact=="HIGH")[c("var_uid", "X2kG_AF", "VARTYPE")], subset(ESP_goi, Impact=="HIGH")[c("var_uid", "ESP_AF", "VARTYPE")], by=c("var_uid", "VARTYPE"), all=T)
both_high$status <- apply(both_high[c("X2kG_AF", "ESP_AF")], 1, function(x) ifelse(!any(is.na(x)), "BOTH", ifelse(is.na(x[1]), "ESP", "X2kG")))
both_high <- droplevels(subset(both_high, status!="BOTH"))
both_high$AF <- apply(both_high[c("X2kG_AF", "ESP_AF")], 1, function(x) x[!is.na(x)])
both_high$VAR <- factor(both_high$VAR=="SNP")
both_high <- subset(both_high, AF< 0.5)
beanplot(AF ~ VARTYPE, data=subset(both_high, status=="ESP"), bw="SJ", log="y", side="both", ll=0.02, col=c("purple", "lightblue", "black"), what=c(0,1,0,1), xlab="Variant Type",
         main="Allele frequence of novel variants in ESP")

require(beanplot)

#p <- ggplot( aes(VARTYPE, ESP_AF), data=ESP_only ) + scale_y_log10(breaks=0.1^seq(0, 5))
#p <- p + geom_jitter(colour="red",position = position_jitter(width = .2, height=.05), alpha=0.8)
#p <- p + theme(panel.background = element_rect(fill = 'white', colour = 'black'))
#p
#ggsave("ESP_novel_AF_by_vartype.png", width=5, height=4)
rm(both_high, ESP_only)
