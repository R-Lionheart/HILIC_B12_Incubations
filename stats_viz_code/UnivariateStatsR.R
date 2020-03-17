library(tidyverse)
options(scipen = 999)


# Upload all required files
BMISd_1_0.2_notnormd <- read.csv("data_processed/IsoLagran1_0.2_notnormd.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)
BMISd_1_5_notnormd <- read.csv("data_processed/IsoLagran1_5_notnormd.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)
BMISd_2_0.2_notnormd <- read.csv("data_processed/IsoLagran2_0.2_notnormd.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)
BMISd_2_5_notnormd <- read.csv("data_processed/IsoLagran2_5_notnormd.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)

# Set user conditions
Condition1 <- "IL1WBT"
Condition2 <- "IL2Control"
SigValue <- "pvalue" # alternative is "qvalue"
file.pattern <- "Cyclonic_5um"
BMISd <- BMISd_2_0.2_notnormd

currentDate <- Sys.Date()

# KRH analysis ------------------------------------------------------------
WBMISd <- BMISd %>%
  spread(key = "Replicate.Name", value = "Adjusted.Area")
WBMISd <- WBMISd[complete.cases(WBMISd), ]
mySamps <- colnames(WBMISd)

# Add Condition1 vs Condition2 stats
myTreat1 <- mySamps[grepl(Condition1, mySamps)]
myTreat2 <- mySamps[grepl(Condition2, mySamps)]
myTreatsdf <- WBMISd[, c(myTreat1, myTreat2)]

WBMISd[, paste(Condition1, "v", Condition2, "_pvalue", sep = "")] <- apply(myTreatsdf, 1, function(x) 
  {t.test(x[myTreat1], x[myTreat2])$p.value}) # add a Pvalue for between the two treatments for QC
WBMISd[, paste(Condition1, "v", Condition2, "_qvalue", sep = "")] <- p.adjust(WBMISd[, ncol(WBMISd)], method = "fdr") 
WBMISd[, paste(Condition1, "v", Condition2, "_FC", sep = "")] <- log2(rowMeans(WBMISd[, myTreat1]) / rowMeans(WBMISd[, myTreat2]))
WBMISd[, paste(Condition1, "_Ave", sep = "")] <- rowMeans(WBMISd[, myTreat1])  
WBMISd[, paste(Condition2, "_Ave", sep = "")] <- rowMeans(WBMISd[, myTreat2])
WBMISd$AveSmp <- rowMeans(WBMISd[, c(myTreat1, myTreat2)])
WBMISd <- WBMISd %>%
  select(Mass.Feature, contains(SigValue), everything()) %>%
  mutate(Significance = ifelse(.[[2]] < 0.1 , "Significant",
                        ifelse(between(.[[2]], 0.1, 0.5), "CloseSig", "NotSig")))

FC_Yaxis <- WBMISd %>%
  select(Mass.Feature, contains("FC")) %>%
  mutate(FC_Yaxis = .[[2]])

dataToPlot <- WBMISd %>%
  left_join(FC_Yaxis) %>%
  select(Mass.Feature, Significance, AveSmp, contains(Condition1), contains(Condition2), contains("FC"))


## Sanity Check
sanitycheck <- BMISd_1_0.2_notnormd %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four")) %>%
  filter(str_detect(SampID, "DSW|Control")) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(myave = mean(Adjusted.Area, na.rm = TRUE))
  
# ggplot(sanitycheck, aes(Mass.Feature, myave, fill = SampID)) +
#   geom_bar(stat = "identity", position = "dodge")


# Condition1 v Condition 2 Significance
a <- ggplot(dataToPlot, aes(x = AveSmp, y = FC_Yaxis, fill = Significance, 
                          label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +
  scale_fill_manual(values = c("lightskyblue3", "grey", "royalblue")) +
  scale_alpha_manual(values = c(1, 0.5, 0.7)) +
  scale_x_log10() +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2)) +
  theme(plot.title = element_text(size = 15),
        legend.position="left",
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        axis.text=element_text(size=8)) +
  labs(x="Average peak size", y=paste("Log2", Condition1, "/", Condition2, sep = "")) +
  theme(legend.position="right") +
  scale_y_continuous(limits=c(-3, 3)) +
  geom_text(data = subset(dataToPlot, Significance == "Significant") ,
            hjust = "inward", nudge_x = 0.05, check_overlap = TRUE)
a
figureFileName <- paste("Figure_", file.pattern, "_", 
                     Condition1, "v", Condition2, "_", 
                     currentDate, ".png", sep = "")

ggsave(filename = figureFileName, plot = a, path = "figures/")

# Return manageable dataset for checking
toDownload <- dataToPlot %>%
  select(Mass.Feature, Significance, AveSmp)

csvFileName <- paste("data_processed/UniVarSig_", file.pattern, "_", 
                     Condition1, "v", Condition2, "_", currentDate, ".csv", sep = "")

write.csv(toDownload, csvFileName, row.names = FALSE)


