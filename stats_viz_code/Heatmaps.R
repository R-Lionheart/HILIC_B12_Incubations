library(tidyverse)
source("src/B12_Functions.R")

file.pattern = "IL"
Norm.Type = "Chl-A Normalized"

filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = file.pattern))
filepath <- file.path("data_processed", paste(filenames, ".csv", sep = ""))

for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE))
}

# Assign normalized file here
data <- IL1_5um_Chla_normd_std

# Rearrange to required format
colnames(data)[1] <- "Replicate.Name"
data2 <- t(data)
colnames(data2) <- as.character(unlist(data2[1,]))
data2 <- as.data.frame(data2[-1,])

data3 <- data.frame(row.names(data2), data2, row.names = NULL)
colnames(data3)[1] <- "Mass.Feature"

for (i in c(2:ncol(data3))) {
  data3[, i] <- as.numeric(as.character(data3[, i]))
}

data3 <- data3 %>%
  mutate(Mass.Feature = as.character(Mass.Feature))

# HEATMAPS for non-Chl standardized data -------------------------------------------------------------------
if (str_detect(file.pattern, "IsoLagran")) {
  
  data4 <- data3 %>%
    pivot_longer(cols = starts_with("X"),
                 names_to = "Replicate.Name",
                 values_to = "Area.BMISd.Normd")
  
  data4$Replicate.Name <- gsub("^.{0,1}", "", data4$Replicate.Name)
  
  heatmap.data <- data4 %>%
    separate(Replicate.Name, into = c("Date", "runtype", "SampID", "replicate")) %>%
    select(Mass.Feature, SampID, Area.BMISd.Normd) %>%
    group_by(Mass.Feature, SampID) %>%
    mutate(Average.Area.Normd = mean(Area.BMISd.Normd)) %>%
    select(Mass.Feature, SampID, Average.Area.Normd) %>%
    unique()

} else if (str_detect(file.pattern, "IL")) {
  
  data4 <- data3 %>%
    pivot_longer(cols = starts_with("IL"),
                 names_to = "Replicate.Name",
                 values_to = "Area.BMISd.Normd")
  
  heatmap.data <- data4 %>%
    separate(Replicate.Name, into = c("SampID", "replicate")) %>%
    select(Mass.Feature, SampID, Area.BMISd.Normd) %>%
    group_by(Mass.Feature, SampID) %>%
    mutate(Average.Area.Normd = mean(Area.BMISd.Normd)) %>%
    select(Mass.Feature, SampID, Average.Area.Normd) %>%
    unique()
}


# Select appropriate heatmap factor levels
# 1_0.2
heatmap.data$SampID <- factor(heatmap.data$SampID,
                              levels = c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW",
                                         "IL1DMBnoBT", "IL1noBT"))
# 1_5
heatmap.data$SampID <- factor(heatmap.data$SampID,
                              levels = c("IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um",
                                         "IL1DMBnoBT5um", "IL1noBT5um"))
# 2_0.2
heatmap.data$SampID <- factor(heatmap.data$SampID,
                              levels = c("IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW",
                                         "IL2DMBnoBT", "IL2noBT"))
# 2_5
heatmap.data$SampID <- factor(heatmap.data$SampID,
                              levels = c("IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", 
                                         "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um")) 

heatmap_1_5 <- ggplot(data = heatmap.data, aes(x = Mass.Feature, y = SampID, fill = Average.Area.Normd)) + 
  geom_tile(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  scale_y_discrete(limits = rev(levels(as.factor(heatmap.data$SampID)))) +
  ggtitle(paste("B12 HILIC Incubation: Cyclonic Eddy 5um", Norm.Type))
print(heatmap_1_5)

require(gridExtra)
grid.arrange(test.heatmap_1_0.2, test.heatmap_1_5, test.heatmap_2_0.2, test.heatmap_2_5, ncol=2)
