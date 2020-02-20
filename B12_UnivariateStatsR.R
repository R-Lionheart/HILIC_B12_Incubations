library(tidyverse)
options(scipen = 999)
BMISd <- read.csv("data_processed/IsoLagran1_0.2_notnormed.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)

# KRH analysis ------------------------------------------------------------
WBMISd <- BMISd %>%
  spread(key = "Replicate.Name", value = "Adjusted.Area")
WBMISd <- WBMISd[complete.cases(WBMISd), ]
mySamps <- colnames(WBMISd)

# Add B12 vs noB12 stats
myTreat1 <- mySamps[grepl("IL1WBT", mySamps)]
myTreat2 <- mySamps[grepl("IL1noBT", mySamps)]
myTreatsdf <- WBMISd[, c(myTreat1, myTreat2)]

WBMISd$WvnoB12_pvalue <- apply(myTreatsdf, 1, function(x) {t.test(x[myTreat1], x[myTreat2])$p.value}) # add a Pvalue for between the two treatments for QC
WBMISd$WvnoB12_qvalue <- p.adjust(WBMISd$WvnoB12_pvalue, method = "fdr") # This corrects for false discovery
WBMISd <- WBMISd %>%
  mutate(WvnoB12_FC = log2(rowMeans(WBMISd[, myTreat1]) / rowMeans(WBMISd[, myTreat2]))) %>%
  mutate(WB12_ave = rowMeans(WBMISd[, myTreat1])) %>%
  mutate(noB12_ave = rowMeans(WBMISd[, myTreat2])) %>%
  mutate(B12vnoB12Sig = WvnoB12_qvalue < 0.1,
         AveSmp = rowMeans(WBMISd[, c(myTreat1, myTreat2)]))

# Add DMB vs DMBnoB12 stats
myTreat1 <- mySamps[grepl("DMB_", mySamps)]
myTreat2 <- mySamps[grepl("DMBnoBT", mySamps)]
myTreatsdf <- WBMISd[, c(myTreat1, myTreat2) ]
WBMISd$WDMBvnoDMB_pvalue <- apply(myTreatsdf, 1, function(x) {t.test(x[myTreat1], x[myTreat2])$p.value}) #add a Pvalue for between the two treatments for QC
WBMISd$WDMBvnoDMB_qvalue <- p.adjust(WBMISd$WDMBvnoDMB_pvalue, method = "fdr") #This corrects for false discovery
WBMISd <- WBMISd %>%
  mutate(WDMBvnoDMB_FC = log2(rowMeans(WBMISd[, myTreat1])/rowMeans(WBMISd[, myTreat2])))%>%
  mutate(WDMBave = rowMeans(WBMISd[,myTreat1])) %>%
  mutate(noDMBave = rowMeans(WBMISd[,myTreat2])) %>%
  mutate(DMBvnoDMBSig = WDMBvnoDMB_qvalue < 0.1)


# Add T0 to Tfinal
myTreat1 <- mySamps[grepl("IT0", mySamps)]
myTreat1 <- myTreat1[grepl("171002", myTreat1)]
myTreat2 <- mySamps[grepl("IL1Control", mySamps)]
myTreatsdf <- WBMISd[, c(myTreat1, myTreat2)]
WBMISd$T0vTF_pvalue <- apply(myTreatsdf, 1, function(x) {t.test(x[myTreat1], x[myTreat2])$p.value}) #add a Pvalue for between the two treatments for QC
WBMISd$T0vTF_qvalue <- p.adjust(WBMISd$T0vTF_pvalue, method = "fdr") #This corrects for false discovery
WBMISd <- WBMISd %>%
  mutate(T0vTF_FC = log2(rowMeans(WBMISd[, myTreat1])/rowMeans(WBMISd[, myTreat2])))%>%
  mutate(T0_ave = rowMeans(WBMISd[,myTreat1])) %>%
  mutate(TF_ave = rowMeans(WBMISd[,myTreat2])) %>%
  mutate(T0vTFSig = T0vTF_qvalue < 0.1)  


#Add T0 to DSW
myTreat1 <- mySamps[grepl("IT0", mySamps)]
myTreat1 <- myTreat1[grepl("171002", myTreat1)]
myTreat2 <- mySamps[grepl("DSW", mySamps)]
myTreatsdf <- WBMISd[, c(myTreat1, myTreat2)]
WBMISd$T0vDSW_pvalue <- apply(myTreatsdf, 1, function(x) {t.test(x[myTreat1],x[myTreat2])$p.value}) #add a Pvalue for between the two treatments for QC
WBMISd$T0vDSW_qvalue <- p.adjust(WBMISd$T0vDSW_pvalue, method = "fdr") #This corrects for false discovery
WBMISd <- WBMISd %>%
  mutate(T0vDSW_FC = log2(rowMeans(WBMISd[,myTreat1])/rowMeans(WBMISd[,myTreat2])))%>%
  mutate(T02_ave = rowMeans(WBMISd[,myTreat1])) %>%
  mutate(DSW_ave = rowMeans(WBMISd[,myTreat2])) %>%
  mutate(T0vDSWSig = T0vDSW_qvalue < 0.1)  

write.csv(WBMISd, "data_processed/WBMISd_wStats.csv")

# originally exported as filtered_wide_combined_stats

