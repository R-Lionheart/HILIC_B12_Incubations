library(tidyverse)
source("src/B12_Functions.R")

# Import all files --------------------------------------------------
filenames <- RemoveCsv(list.files(path = 'data_processed', pattern = "UniVarSig"))

for (i in filenames) {
  filepath <- file.path('data_processed', paste(i,".csv", sep = ""))
  assign(i, read.csv(filepath, stringsAsFactors = FALSE))
}
#######
dsw1_0.2 <- `UniVarSig_Cyclonic_0.2um_IL1DSWvIL1Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig1_0.2_dsw = Significance)

dsw1_5 <- `UniVarSig_Cyclonic_5um_IL1DSWvIL1Control_2020-03-13`%>%
  select(Mass.Feature, Significance) %>%
  rename(Sig1_5_dsw = Significance)

dsw2_0.2 <- `UniVarSig_Anticyclonic_0.2um_IL2DSWvIL2Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig2_0.2_dsw = Significance)

dsw2_5 <- `UniVarSig_Anticyclonic_5um_IL2DSWvIL2Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig2_5_dsw = Significance)

fullDSW <- dsw1_0.2 %>%
  left_join(dsw1_5) %>%
  left_join(dsw2_0.2) %>%
  left_join(dsw2_5) 

########
noBT1_0.2 <- `UniVarSig_Cyclonic_0.2um_IL1noBTvIL1Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig1_0.2_noBT = Significance)

noBT1_5 <- `UniVarSig_Cyclonic_5um_IL1noBTvIL1Control_2020-03-13`%>%
  select(Mass.Feature, Significance) %>%
  rename(Sig1_5_noBT = Significance)

noBT2_0.2 <- `UniVarSig_Anticyclonic_0.2um_IL2noBTvIL2Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig2_0.2_noBT = Significance)

noBT2_5 <- `UniVarSig_Anticyclonic_5um_IL2noBTvIL2Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig2_5_noBT = Significance)

fullBT <- noBT1_0.2 %>%
  left_join(noBT1_5) %>%
  left_join(noBT2_0.2) %>%
  left_join(noBT2_5)


#######
DMB1_0.2 <- `UniVarSig_Cyclonic_0.2um_IL1DMBvIL1Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig1_0.2_DMB = Significance)

DMB1_5 <- `UniVarSig_Cyclonic_5um_IL1DMBvIL1Control_2020-03-13`%>%
  select(Mass.Feature, Significance) %>%
  rename(Sig1_5_DMB = Significance)

DMB2_0.2 <- `UniVarSig_Anticyclonic_0.2um_IL2DMBvIL2Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig2_0.2_DMB = Significance)

DMB2_5 <- `UniVarSig_Anticyclonic_5um_IL2DMBvIL2Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig2_5_DMB = Significance)

fullDMB <- DMB1_0.2 %>%
  left_join(DMB1_5) %>%
  left_join(DMB2_0.2) %>%
  left_join(DMB2_5)

######
DMBnoBT1_0.2 <- `UniVarSig_Cyclonic_0.2um_IL1DMBnoBTvIL1Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig1_0.2_DMBnoBT = Significance)

DMBnoBT1_5 <- `UniVarSig_Cyclonic_5um_IL1DMBnoBTvIL1Control_2020-03-13`%>%
  select(Mass.Feature, Significance) %>%
  rename(Sig1_5_DMBnoBT = Significance)

DMBnoBT2_0.2 <- `UniVarSig_Antiyclonic_0.2um_IL2DMBnoBTvIL2Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig2_0.2_DMBnoBT = Significance)

DMBnoBT2_5 <- `UniVarSig_Antiyclonic_5um_IL2DMBnoBTvIL2Control_2020-03-13` %>%
  select(Mass.Feature, Significance) %>%
  rename(Sig2_5_DMBnoBT = Significance)

fullDMBnoBT <- DMBnoBT1_0.2 %>%
  left_join(DMBnoBT1_5) %>%
  left_join(DMBnoBT2_0.2) %>%
  left_join(DMBnoBT2_5)

####
full.complete <- fullBT %>%
  left_join(fullDMB) %>%
  left_join(fullDSW) %>%
  left_join(fullDMBnoBT)



complete <- read.csv("~/Downloads/AllSigs_Complete.csv", stringsAsFactors = FALSE)

# No cases of compounds that are completely insignificant
Notsigs <- complete %>%
  # filter(is.na(Sig1_0.2_noBT) | Sig1_0.2_noBT == "NotSig") %>%
  # filter(is.na(Sig1_5_noBT) | Sig1_5_noBT == "NotSig") %>%
  # filter(is.na(Sig2_0.2_noBT) | Sig2_0.2_noBT == "NotSig")
  filter(is.na(.[[3]]) | .[[3]] == "NotSig") 

Close_notSigs <- complete %>%
  filter(is.na(.[[3]]) | .[[3]] %in% c("NotSig", "CloseSig"))
  
  
