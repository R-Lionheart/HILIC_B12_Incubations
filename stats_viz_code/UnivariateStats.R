library(tidyverse)
options(scipen = 999)

# Upload all required files
uploadFiles <- function(myfilepath) {
  # Function for uploading non-standardized files.
  #
  # Returns: dataframe with relevant columns selected.
  uploaded.df <- read.csv(myfilepath, stringsAsFactors = FALSE) %>%
    select(Mass.Feature:Adjusted.Area)
  
  return(uploaded.df)
}

BMISd_1_0.2_notnormd <- uploadFiles("data_processed/IsoLagran1_0.2_notnormd.csv")
BMISd_1_5_notnormd <- uploadFiles("data_processed/IsoLagran1_5_notnormd.csv")
BMISd_2_0.2_notnormd <- uploadFiles("data_processed/IsoLagran2_0.2_notnormd.csv")
BMISd_2_5_notnormd <- uploadFiles("data_processed/IsoLagran2_5_notnormd.csv")

# Set filtering conditions that correspond to the treatments you are comparing.
Condition1 <- "IL1WBT" # Other options: IL1DMBnoBT, IL2WBT, IL1noBt, etc.
Condition2 <- "IL1Control"
SigValue <- "pvalue" # alternative is "qvalue", when using fdr-corrected values.
file.pattern <- "Cyclonic_0.2um" # will be used as a search ID and title for graphs 
SigNumber <- 0.1 # Pvalue cutoff
BMISd <- BMISd_1_0.2_notnormd # Assign correct dataframe for analysis

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

# Add a Pvalue for between the two treatments for QC
WBMISd[, paste(Condition1, "v", Condition2, "_pvalue", sep = "")] <- apply(myTreatsdf, 1, function(x) 
  {t.test(x[myTreat1], x[myTreat2])$p.value}) 
# Add a false-discovery-rate-corrected q value
WBMISd[, paste(Condition1, "v", Condition2, "_qvalue", sep = "")] <- p.adjust(WBMISd[, ncol(WBMISd)], method = "fdr") 
# Calculate fold change: Condition 1 / Condition 2
WBMISd[, paste(Condition1, "v", Condition2, "_FC", sep = "")] <- log2(rowMeans(WBMISd[, myTreat1]) / rowMeans(WBMISd[, myTreat2]))
# Calculate Condition 1 Average
WBMISd[, paste(Condition1, "_Ave", sep = "")] <- rowMeans(WBMISd[, myTreat1])  
# Calculate Condition 2 Average
WBMISd[, paste(Condition2, "_Ave", sep = "")] <- rowMeans(WBMISd[, myTreat2])
# Calculate complete row means
WBMISd$AveSmp <- rowMeans(WBMISd[, c(myTreat1, myTreat2)])
# Organize columns and assign significance
WBMISd <- WBMISd %>%
  select(Mass.Feature, contains(SigValue), everything()) %>%
  mutate(Significance = ifelse(.[[2]] < SigNumber, "Significant",
                        ifelse(between(.[[2]], SigNumber, 0.5), "CloseSig", "NotSig")))

# Adjust fold change axis
FC_Yaxis <- WBMISd %>%
  select(Mass.Feature, contains("FC")) %>%
  mutate(FC_Yaxis = .[[2]]) # Multiply by -1 here to reverse y axis

# Combine for plot.
dataToPlot <- WBMISd %>%
  left_join(FC_Yaxis) %>%
  select(Mass.Feature, Significance, AveSmp, contains(Condition1), contains(Condition2), contains("FC"))


## Sanity Check for fold change ratios
sanitycheck <- BMISd %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four")) %>%
  filter(SampID == Condition1 | SampID == Condition2) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(myave = mean(Adjusted.Area, na.rm = TRUE))
  
ggplot(sanitycheck, aes(Mass.Feature, myave, fill = SampID)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90))


# Condition1 v Condition 2 Significance
SignificancePlot <- ggplot(dataToPlot, aes(x = AveSmp, y = FC_Yaxis, fill = Significance, 
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
  geom_text(data = subset(dataToPlot, Significance == "Significant"),
            hjust = "inward", nudge_x = 0.05, check_overlap = TRUE) +
  geom_text(data = subset(dataToPlot, Significance == "CloseSig"),
            hjust = "inward", nudge_x = 0.05, check_overlap = TRUE,
            color = "lightskyblue3" )
SignificancePlot

figureFileName <- paste("Figure_", file.pattern, "_", 
                     Condition1, "v", Condition2, "_", 
                     currentDate, ".png", sep = "")

ggsave(filename = figureFileName, plot = SignificancePlot, path = "figures/")

# Return manageable dataset for checking
toDownload <- dataToPlot %>%
  select(Mass.Feature, Significance, AveSmp)

csvFileName <- paste("data_processed/UniVarSig_", file.pattern, "_", 
                     Condition1, "v", Condition2, "_", currentDate, ".csv", sep = "")

write.csv(toDownload, csvFileName, row.names = FALSE)


