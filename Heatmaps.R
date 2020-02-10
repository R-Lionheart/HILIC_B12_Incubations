library(tidyverse)


BMISd <- read.csv("data_processed/IsoLagran_0.2_notnormed.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)

# HEATMAPS -------------------------------------------------------------------
data <- read.csv("data_processed/IsoLagran_0.2_normed.csv", check.names = FALSE, stringsAsFactors = FALSE)
colnames(data)[1] <- "Replicate.Name"
data2t <- t(data)
colnames(data2t) <- as.character(unlist(data2t[1,]))
data2t <- as.data.frame(data2t[-1,])

data3 <- data.frame(row.names(data2t), data2t, row.names = NULL)
colnames(data3)[1] <- "Mass.Feature"

for (i in c(2:ncol(data3))) {
  data3[, i] <- as.numeric(as.character(data3[, i]))
}

data3 <- data3 %>%
  mutate(Mass.Feature = as.character(Mass.Feature))

data4 <- data3 %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "Replicate.Name",
               values_to = "Area.BMISd.Normd")

data4$Replicate.Name <- gsub("^.{0,1}", "", data4$Replicate.Name)

heatmap.data <- data4 %>%
  separate(Replicate.Name, into = c("Date", "runtype", "SampID", "replicate")) %>%
  select(Mass.Feature, SampID, Area.BMISd.Normd) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Average.Area = mean(Area.BMISd.Normd)) %>%
  select(Mass.Feature, SampID, Average.Area) %>%
  unique()

# Select appropriate heatmap factor levels
heatmap.data$SampID <- factor(heatmap.data$SampID,
                              levels = c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW",
                                         "IL1DMBnoBT", "IL1noBT"))

heatmap.data$SampID <- factor(heatmap.data$SampID,
                              levels = c("IL1IT05um", "IL1Control5m", "IL1DMB5um", "IL1WBT5m", "IL1DSW5um",
                                         "IL1DMBnoBT", "IL1noBT"))
heatmap.data$SampID <- factor(heatmap.data$SampID,
                              levels = c("IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW",
                                         "IL2DMBnoBT", "IL2noBT"))

heatmap.data$SampID <- factor(heatmap.data$SampID,
                              levels = c("IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", 
                                         "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um")) 

test.heatmap <- ggplot(data = heatmap.data, aes(x = Mass.Feature, y = SampID, fill = Average.Area)) + 
  geom_tile(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        #axis.ticks.length=unit(.25, "cm"),
        strip.text = element_text(size = 10)) +
  scale_y_discrete(limits = rev(levels(as.factor(heatmap.data$SampID)))) +
  ggtitle("B12 HILIC Incubation: IsoLagran 0.2um") 
print(test.heatmap)

require(gridExtra)
grid.arrange(ncol=2)
