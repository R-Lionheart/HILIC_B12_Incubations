library(tidyverse)
library(vegan)
options(scipen = 999)

source("src/B12_Functions.R")

# ## Data normalized to Chlorophyll-a

# UploadFiles <- function(myfilepath) {
#   # Function for uploading non-standardized files.
#   #
#   # Returns: dataframe with relevant columns selected.
#   uploaded.df <- read.csv(myfilepath, stringsAsFactors = FALSE) %>%
#     select(Mass.Feature:Adjusted.Area)
#   
#   return(uploaded.df)
# }
# 
# # Upload all required files
# BMISd_1_0.2_notnormd <- UploadFiles("data_processed/IsoLagran1_0.2_notnormd.csv")
# BMISd_1_5_notnormd <- UploadFiles("data_processed/IsoLagran1_5_notnormd.csv")
# BMISd_2_0.2_notnormd <- UploadFiles("data_processed/IsoLagran2_0.2_notnormd.csv")
# BMISd_2_5_notnormd <- UploadFiles("data_processed/IsoLagran2_5_notnormd.csv")
# 
# # Set filtering conditions that correspond to the treatments you are comparing.
Condition1 <- "IL1DSW5um" # Other options: IL1DMBnoBT, IL2WBT, IL1noBt, etc.
Condition2 <- "IL1Control5um"
FilterSize <- "5um"
ChlaEddy <- "IL1"
SigValue <- "pvalue" # alternative is "qvalue", when using fdr-corrected values.
file.pattern <- "Chla Normalized Cyclonic_5um" # will be used as a search ID and title for graphs
SigNumber <- 0.1 # Pvalue cutoff
BMISd <- BMISd_1_5_notnormd # Assign correct dataframe for analysis.

# ## Upload and recode chlorophyll data
# Chlorophyll <- read.csv("data_raw/Dyhrman_MS_Chla.csv", stringsAsFactors = FALSE) %>%
#   select(1:6) %>%
#   rename(Filter.size = Filter..µm.,
#          Chla = Chl.a..µg.L.) %>%
#   mutate(Date = as.character(Date),
#          Eddy = ifelse(str_detect(Date, "13|9"), "IL2", "IL1")) %>%
#   filter(str_detect(Eddy, ChlaEddy)) %>%
#   unite("Replicate.Name", SAMPLE, Eddy, sep = "_", remove = TRUE) %>%
#   filter(Filter.size == 5.0) %>%
#   mutate(Chla = as.numeric(Chla))

# Fix Replicate Name labels -----------------------------------------------
# Chlorophyll_fixed <- Chlorophyll %>%
#   mutate(Replicate.Name = recode(Replicate.Name,
#                                  "C-1_IL1" ="171009_Smp_IL1Control5um_1",
#                                  "C-1_IL2" = "171009_Smp_IL2Control5um_1",
#                                  "C-2_IL1" = "171009_Smp_IL1Control5um_2",
#                                  "C-2_IL2" = "171009_Smp_IL2Control5um_2",
#                                  "C-3_IL1" = "171009_Smp_IL1Control5um_3",
#                                  "C-3_IL2" = "171009_Smp_IL2Control5um_3",
#                                  "ØBT-1_IL1" = "171009_Smp_IL1noBT5um_1",
#                                  "ØBT-1_IL2" = "171009_Smp_IL2noBT5um_1",
#                                  "ØBT-2_IL1" = "171009_Smp_IL1noBT5um_2",
#                                  "ØBT-2_IL2" = "171009_Smp_IL2noBT5um_2",
#                                  "ØBT-3_IL1" = "171009_Smp_IL1noBT5um_3",
#                                  "ØBT-3_IL2" = "171009_Smp_IL2noBT5um_3",
#                                  "DMB-1_IL1" = "171009_Smp_IL1DMB5um_1",
#                                  "DMB-1_IL2" = "171009_Smp_IL2DMB5um_1",
#                                  "DMB-2_IL1" = "171009_Smp_IL1DMB5um_2",
#                                  "DMB-2_IL2" = "171009_Smp_IL2DMB5um_2",
#                                  "DMB-3_IL1" = "171009_Smp_IL1DMB5um_3",
#                                  "DMB-3_IL2" = "171009_Smp_IL2DMB5um_3",
#                                  "DMBØBT-1_IL1" = "171009_Smp_IL1DMBnoBT5um_1",
#                                  "DMBØBT-1_IL2" = "171009_Smp_IL2DMBnoBT5um_1",
#                                  "DMBØBT-2_IL1" = "171009_Smp_IL1DMBnoBT5um_2",
#                                  "DMBØBT-2_IL2" = "171009_Smp_IL2DMBnoBT5um_2",
#                                  "DMBØBT-3_IL1" = "171009_Smp_IL1DMBnoBT5um_3",
#                                  "DMBØBT-3_IL2" = "171009_Smp_IL2DMBnoBT5um_3",
#                                  "DSW-1_IL1" = "171009_Smp_IL1DSW5um_1",
#                                  "DSW-1_IL2" = "171009_Smp_IL2DSW5um_1",
#                                  "DSW-2_IL1" = "171009_Smp_IL1DSW5um_2",
#                                  "DSW-2_IL2" = "171009_Smp_IL2DSW5um_2",
#                                  "DSW-3_IL1" = "171009_Smp_IL1DSW5um_3",
#                                  "DSW-3_IL2" = "171009_Smp_IL2DSW5um_3",
#                                  "IS-1_IL1" = "171009_Smp_IL1IT05um_1",
#                                  "IS-1_IL2" = "171009_Smp_IL2IT05um_1",
#                                  "IS-2_IL1" = "171009_Smp_IL1IT05um_2",
#                                  "IS-2_IL2" = "171009_Smp_IL2IT05um_2",
#                                  "IS-3_IL1" = "171009_Smp_IL1IT05um_3",
#                                  "IS-3_IL2" = "171009_Smp_IL2IT05um_3",
#                                  "WBT-1_IL1" = "171009_Smp_IL1WBT5um_1",
#                                  "WBT-1_IL2" = "171009_Smp_IL2WBT5um_1",
#                                  "WBT-2_IL1" = "171009_Smp_IL1WBT5um_2",
#                                  "WBT-2_IL2" = "171009_Smp_IL2WBT5um_2",
#                                  "WBT-3_IL1" = "171009_Smp_IL1WBT5um_3",
#                                  "WBT-3_IL2" = "171009_Smp_IL2WBT5um_3")) %>%
#   separate(Replicate.Name, into = c("one", "two", "SampID", "replicate")) %>%
#   select(SampID, replicate, Chla) %>%
#   unite(SampID, replicate, col = "Replicate.Name")

# Fix BMISd names to match ------------------------------------------------------------
# BMISd_fixed <- BMISd %>%
#   separate(Replicate.Name, into = c("one", "two", "SampID", "replicate")) %>%
#   select(-c("one", "two")) %>%
#   unite(SampID, replicate, col = "Replicate.Name")
# 
# # Normalize to chlorophyll ------------------------------------------------------------
# complete.set <- Chlorophyll_fixed %>%
#   select(Replicate.Name, Chla) %>%
#   left_join(BMISd_fixed) %>%
#   mutate(Chla = as.numeric(Chla)) %>%
#   mutate(Normalized.by.Chla = Adjusted.Area/Chla) %>%
#   filter(!Chla < 0) %>% # Remove negative Chl-A values
#   select(Mass.Feature, Replicate.Name, Normalized.by.Chla) %>%
#   na.omit()
# 
# ## Standardize
# complete.wide <- ChlMakeWide(complete.set)
# complete.wide[is.na(complete.wide)] <- 1000
# complete.wideT <- t(complete.wide)
# 
# complete.wide.normalizedT <- decostand(complete.wideT, method = "standardize", na.rm = TRUE) 



## Save Chlorophyll files to data_processed/
# write.csv(complete.wide.normalizedT, paste("data_processed/", ChlaEddy, "_", FilterSize, "_ChlA_normd_std.csv", sep = ""))
# write.csv(complete.set, paste("data_processed/", ChlaEddy, "_", FilterSize, "_ChlA_normd_nostd.csv", sep = ""))
# write.csv(Chlorophyll_fixed, "data_processed/Chlorophyll_fixed.csv")


Chl.pattern <- "ChlAnormd"
# # Import Chlorophyll ------------------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = regex(Chl.pattern, ignore_case = TRUE)))

for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

#############################################################################################
# Log normalize the chlorophyll-normalized data ---------------------------

complete.set <- IL1_5um_ChlAnormd_notstd %>%
  select (-X)

complete.set.wide <- ChlMakeWide(complete.set) %>%
  na.omit()

# complete.set.wide <- complete.set %>%
#   tidyr::spread(Replicate.Name, Normalized.by.Chla) 
# complete.set.wide <- column_to_rownames(complete.set.wide, "Mass.Feature")
complete.set.normalized <- decostand(complete.set.wide, method = "log", na.rm = TRUE)

#final.set <- complete.set.normalized[complete.cases(complete.set.normalized), ]
mySamps.normd <- colnames(complete.set.normalized)

# Add Condition 1 v Condition 2 stats -------------------------------------
myTreat1 <- mySamps.normd[grepl(Condition1, mySamps.normd)]
myTreat2 <- mySamps.normd[grepl(Condition2, mySamps.normd)]
myTreatsdf <- final.set[, c(myTreat1, myTreat2)]

# Add a Pvalue for between the two treatments for QC
final.set[, paste(Condition1, "v", Condition2, "_pvalue", sep = "")] <- apply(myTreatsdf, 1, function(x) 
{t.test(x[myTreat1], x[myTreat2])$p.value}) 
# Add a false-discovery-rate-corrected q value
final.set[, paste(Condition1, "v", Condition2, "_qvalue", sep = "")] <- p.adjust(final.set[, ncol(final.set)], method = "fdr") 
# Calculate fold change: Condition 1 / Condition 2
final.set[, paste(Condition1, "v", Condition2, "_FC", sep = "")] <- log2(rowMeans(final.set[, myTreat1]) / rowMeans(final.set[, myTreat2]))
# Calculate Condition 1 Average
final.set[, paste(Condition1, "_Ave", sep = "")] <- rowMeans(final.set[, myTreat1])  
# Calculate Condition 2 Average
final.set[, paste(Condition2, "_Ave", sep = "")] <- rowMeans(final.set[, myTreat2])
# Calculate complete row means
final.set$AveSmp <- rowMeans(final.set[, c(myTreat1, myTreat2)])
# Organize columns and assign significance
final.set <- final.set %>%
  rownames_to_column(var = "Mass.Feature") %>%
  select(Mass.Feature, contains(SigValue), everything()) %>%
  mutate(Significance = ifelse(.[[2]] < SigNumber, "Significant", "NotSig"))

# Adjust fold change axis
FC_Yaxis.normd <- final.set %>%
  select(Mass.Feature, contains("FC")) %>%
  mutate(FC_Yaxis = .[[2]]) # Multiply by -1 here to reverse y axis

# Combine for plot
dataToPlot.normed <- final.set %>%
  left_join(FC_Yaxis.normd) %>%
  select(Mass.Feature, Significance, AveSmp, contains(Condition1), contains(Condition2), contains("FC"))

# Normalized Sanity Check
sanitycheck.normd <- myTreatsdf %>%
  rownames_to_column(var = "Mass.Feature")
sanitycheck.normd <- sanitycheck.normd %>%
  pivot_longer(cols = starts_with("IL"),
               names_to = "SampID",
               values_to = "Area.BMISd.Normd") %>%
  separate(SampID, into = c("SampID", "drop")) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Normalized.Average = mean(Area.BMISd.Normd)) %>%
  select(-c("drop", "Area.BMISd.Normd")) %>%
  unique()

sanitycheck.normd.plot <- ggplot(sanitycheck.normd, aes(x = reorder(Mass.Feature, -Normalized.Average), y = Normalized.Average,
                                  fill = SampID)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 11)) +
  xlab("Log-normalized Area") +
  ggtitle(paste(file.pattern, "Log Normalized", Condition1, "v", Condition2))
sanitycheck.normd.plot

# Log normalized Condition1 v Condition 2 Significance plots -------------------------------------
normdPlot <- ggplot(dataToPlot.normed, aes(x = AveSmp, y = FC_Yaxis, fill = Significance, 
                               label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +
  scale_fill_manual(values = c("grey", "royalblue")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_x_log10() +
  ggtitle(paste("Chl-A Log Normalized", Condition1, "v", Condition2)) +
  theme(plot.title = element_text(size = 25),
        legend.position="left",
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text=element_text(size=10)) +
  labs(x="Average peak size", y=paste("Log2", Condition1, "/", Condition2, sep = "")) +
  theme(legend.position="right") +
  scale_y_continuous(limits=c(-7, 7)) +
  geom_text(data = subset(dataToPlot.normed, Significance == "Significant"),
            hjust = "inward", nudge_x = 0.05, check_overlap = TRUE, size = 6) 
normdPlot

#############################################################################################

# Univariate analysis NO LOG NORMED ------------------------------------------------------------
WBMISd <- complete.set %>%
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
sanitycheck <- complete.set %>%
  separate(Replicate.Name, into = c("SampID", "replicate"), sep = "_") %>%
  filter(SampID == Condition1 | SampID == Condition2) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(myave = mean(Normalized.by.Chla, na.rm = TRUE))

ggplot(sanitycheck, aes(Mass.Feature, myave, fill = SampID)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 11)) +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2))


# Condition1 v Condition 2 Significance
SignificancePlot <- ggplot(dataToPlot, aes(x = AveSmp, y = FC_Yaxis, fill = Significance, 
                                           label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +
  scale_fill_manual(values = c("grey", "royalblue")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_x_log10() +
  ggtitle(paste("Chl-A", Condition1, "v", Condition2)) +
  theme(plot.title = element_text(size = 25),
        legend.position="left",
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text=element_text(size=10)) +
  labs(x="Average peak size", y=paste("Log2", Condition1, "/", Condition2, sep = "")) +
  theme(legend.position="right") +
  scale_y_continuous(limits=c(-7, 7)) +
  geom_text(data = subset(dataToPlot, Significance == "Significant"),
            hjust = "inward", nudge_x = 0.05, check_overlap = TRUE, size = 6) 
SignificancePlot

