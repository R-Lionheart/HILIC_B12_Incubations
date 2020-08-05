library(tidyverse)
options(scipen = 999)

source("src/B12_Functions.R")

# Import Chlorophyll ------------------------------------------------------
Chl.pattern <- "IsoLagran" # "ChlAnormd" for chlorophyll, # "IsoLagran" for non-ChlA

filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = regex(Chl.pattern, ignore_case = TRUE)))
for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE) %>%
           select(-X))
}

# Import any non-standardized files ------------------------------------------------------
# Required for both ChlA and non-ChlA situations
dataset.pattern <- "_notstd"

filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = regex(dataset.pattern, ignore_case = TRUE)))
for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE) %>%
           select(-X))
}

# Set filtering conditions that correspond to the treatments you are comparing.
Condition1 <- "IL1IT05um" # Other options: IL1DMBnoBT, IL2WBT, IL1noBt, etc. Include 5um if using those!
Condition2 <- "IL1DSW5um"
SigValue <- "pvalue" # alternative is "qvalue", when using fdr-corrected values.
file.pattern <- "Cyclonic_5um" # will be used as a search ID and title for graphs. Use "wChlA" if ChlA is involved.
SigNumber <- 0.1 # Pvalue cutoff
BMISd <- IsoLagran1_5_notstd # Assign correct dataframe for analysis. Should be non-standardized.

currentDate <- Sys.Date()


# Create table for analysis -----------------------------------------------
if (grepl("wChlA", file.pattern)) {
  WBMISd <- BMISd %>%
    spread(key = "Replicate.Name", value = "Normalized.by.Chla")
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
    mutate(Significance = ifelse(.[[2]] < SigNumber, "Significant", "NotSig"))
  
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
    separate(Replicate.Name, into = c("SampID", "four")) %>%
    filter(SampID == Condition1 | SampID == Condition2) %>%
    group_by(Mass.Feature, SampID) %>%
    mutate(myave = mean(Normalized.by.Chla, na.rm = TRUE))
  
} else {
  
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
    # mutate(Significance = ifelse(.[[2]] < SigNumber, "Significant",
    #                       ifelse(between(.[[2]], SigNumber, 0.5), "CloseSig", "NotSig")))
    mutate(Significance = ifelse(.[[2]] < SigNumber, "Significant", "NotSig"))
  
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
}

# Plotting section --------------------------------------------------------

# Plot non-transformed data
if (grepl("wChlA", file.pattern)) {
  ggplot(sanitycheck, aes(x = reorder(Mass.Feature, -Normalized.by.Chla), 
                          y = Normalized.by.Chla, fill = SampID)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5)) +
    ggtitle(paste(file.pattern, Condition1, "v", Condition2))
} else {
  ggplot(sanitycheck, aes(x = reorder(Mass.Feature, -myave), 
                          y = myave, fill = SampID)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5)) +
    ggtitle(paste(file.pattern, Condition1, "v", Condition2))
}

# Plot log-normalized data
if (grepl("wChlA", file.pattern)) {
  logNormalized <- sanitycheck %>%
    mutate(log.normed = log(Normalized.by.Chla + 1))
} else {
  logNormalized <- sanitycheck %>%
    mutate(log.normed = log(myave + 1))
}

ggplot(logNormalized, aes(x = reorder(Mass.Feature, -log.normed), 
                        y = log.normed, fill = SampID)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5)) +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2, "_logNormalized"))

# Plot square root data
if (grepl("wChlA", file.pattern)) {
  squareRoot <- sanitycheck %>%
    mutate(square.root = sqrt(Normalized.by.Chla))
} else {
  squareRoot <- sanitycheck %>%
    mutate(square.root = sqrt(myave))
}

ggplot(squareRoot, aes(x = reorder(Mass.Feature, -square.root), 
                       y = square.root, fill = SampID)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5)) +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2, "_SquareRoot"))


# Plot fourth root data
if (grepl("wChlA", file.pattern)) {
  fourthRoot <- sanitycheck %>%
    mutate(fourth.root = Normalized.by.Chla ^ 0.25)
} else {
  fourthRoot <- sanitycheck %>%
    mutate(fourth.root = myave ^ 0.25)
}

ggplot(fourthRoot, aes(x = reorder(Mass.Feature, -fourth.root), 
                     y = fourth.root, fill = SampID)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5)) +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2, "_FourthRoot"))


# Condition1 v Condition 2 Significance
SignificancePlot <- ggplot(dataToPlot, aes(x = AveSmp, y = FC_Yaxis, fill = Significance, 
                          label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +
  scale_fill_manual(values = c("grey", "royalblue")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_x_log10() +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2)) +
  theme(plot.title = element_text(size = 25),
        legend.position="left",
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text=element_text(size=10)) +
  labs(x="Average peak size", y=paste("Log2", Condition1, "/", Condition2, sep = "")) +
  theme(legend.position="right") +
  scale_y_continuous(limits=c(-3, 3)) +
  geom_text(data = subset(dataToPlot, Significance == "Significant"),
            hjust = "inward", nudge_x = 0.05, check_overlap = TRUE, size = 6) 
  # geom_text(data = subset(dataToPlot, Significance == "CloseSig"),
  #           hjust = "inward", nudge_x = 0.05, check_overlap = TRUE,
  #           color = "lightskyblue3" )
SignificancePlot

figureFileName <- paste("Figure_", file.pattern, "_", 
                     Condition1, "v", Condition2, "_", 
                     currentDate, ".png", sep = "")

ggsave(filename = figureFileName, plot = SignificancePlot, path = "figures/")



