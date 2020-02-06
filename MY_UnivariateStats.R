#source("src/biostats.R")
library(tidyverse)
#library("vegan")

#Combine all the data to get Long data and make wide
#Get the wide data
dat <- read.csv("KRH_info/LongDat_Normed.csv")

wdat <- dat %>%
  filter(replicate %in% c("AB", "CD", "EF")) %>%
  mutate(SampID_Rep = paste(SampID, replicate, sep = "_")) %>%
  select(Compound.Name, BMISd_Normed, SampID_Rep) %>%
  spread(key = "SampID_Rep", value = "BMISd_Normed")

wdat2 <- wdat 
Samps <- colnames(wdat2)[grepl("B12", colnames(wdat2))] 

#Add Bulk B12 stats
Treat1 <- Samps[grepl("RB12", Samps)]
Treat2 <- Samps[grepl("LB12", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$RvLB12_pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
wdat2$RvLB12_qvalue <- p.adjust(wdat2$RvLB12_pvalue, method = "fdr") #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(RvLB12_FC = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(RB12_ave = rowMeans(wdat2[,Treat1])) %>%
  mutate(LB12_ave = rowMeans(wdat2[,Treat2])) %>%
  mutate(B12Sig = RvLB12_qvalue < 0.05,
         AveSmp = rowMeans(wdat2[,c(Treat1, Treat2)]))

#Add Bulk Light stats
Treat1 <- Samps[grepl("HL", Samps)]
Treat2 <- Samps[grepl("LL", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$HvsLL_pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
wdat2$HvsLL_qvalue <- p.adjust(wdat2$HvsLL_pvalue, method = "fdr") #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(HvsLL_FC = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(HLave = rowMeans(wdat2[,Treat1])) %>%
  mutate(LLave = rowMeans(wdat2[,Treat2])) %>%
  mutate(LightSig = HvsLL_qvalue < 0.05)

#Add B12 stats underHL
Treat1 <- Samps[grepl("RB12HL", Samps)]
Treat2 <- Samps[grepl("LB12HL", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$RvLB12HL_pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
wdat2$RvLB12HL_qvalue <- p.adjust(wdat2$RvLB12HL_pvalue, method = "fdr") #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(RvLB12HL_FC = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(RB12HL_ave = rowMeans(wdat2[,Treat1])) %>%
  mutate(LB12HL_ave = rowMeans(wdat2[,Treat2])) %>%
  mutate(B12HLSig = RvLB12HL_qvalue < 0.05)  


#Add Bulk B12 stats under LL
Treat1 <- Samps[grepl("RB12LL", Samps)]
Treat2 <- Samps[grepl("LB12LL", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$RvLB12LL_pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
wdat2$RvLB12LL_qvalue <- p.adjust(wdat2$RvLB12LL_pvalue, method = "fdr") #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(RvLB12LL_FC = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(RB12LL_ave = rowMeans(wdat2[,Treat1])) %>%
  mutate(LB12LL_ave = rowMeans(wdat2[,Treat2])) %>%
  mutate(B12LLSig = RvLB12LL_qvalue < 0.05)  


#Add  Light stats at RB12
Treat1 <- Samps[grepl("RB12HL", Samps)]
Treat2 <- Samps[grepl("RB12LL", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$HvsLLRB12_pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
wdat2$HvsLLRB12_qvalue <- p.adjust(wdat2$HvsLLRB12_pvalue, method = "fdr") #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(HvsLLRB12_FC = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(HLaveRB12 = rowMeans(wdat2[,Treat1])) %>%
  mutate(LLaveRB12 = rowMeans(wdat2[,Treat2])) %>%
  mutate(LightRB12Sig = HvsLLRB12_qvalue < 0.05)

#Add  Light stats at LB12
Treat1 <- Samps[grepl("LB12HL", Samps)]
Treat2 <- Samps[grepl("LB12LL", Samps)]
Treatsdf <- wdat2[, c(Treat1,Treat2) ]
wdat2$HvsLLLB12_pvalue <- apply(Treatsdf, 1, function(x) {t.test(x[Treat1],x[Treat2])$p.value}) #add a Pvalue for between the two treatments for QC
wdat2$HvsLLLB12_qvalue <- p.adjust(wdat2$HvsLLLB12_pvalue, method = "fdr") #This corrects for false discovery
wdat2 <- wdat2 %>%
  mutate(HvsLLLB12_FC = log2(rowMeans(wdat2[,Treat1])/rowMeans(wdat2[,Treat2])))%>%
  mutate(HLaveLB12 = rowMeans(wdat2[,Treat1])) %>%
  mutate(LLaveLB12 = rowMeans(wdat2[,Treat2])) %>%
  mutate(LightLB12Sig = HvsLLLB12_qvalue < 0.05)


# originally exported as filtered_wide_combined_stats
