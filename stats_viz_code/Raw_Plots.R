source("src/B12_Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(tidyverse) 
library(vegan)

BMIS.pattern = "BMIS_Output"

# BMIS import and filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

HILIC_all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  select(Mass.Feature, runDate:replicate, Adjusted.Area) %>%
  unite(Replicate.Name, runDate:replicate, sep = "_") %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m")) %>%
  filter(!Mass.Feature == "Inj_vol") %>%
  filter(!str_detect(Mass.Feature, ","))

# Adjust for ILT0 naming issues --------------------------------------------
HILIC_all <- HILIC_all %>%
  mutate(Replicate.Name = recode(Replicate.Name, 
                                 "171002_Smp_IT0_1" ="171002_Smp_IL1IT0_1", 
                                 "171002_Smp_IT0_2" = "171002_Smp_IL1IT0_2",
                                 "171002_Smp_IT0_3" = "171002_Smp_IL1IT0_3",
                                 "171009_Smp_IT05um_1" = "171009_Smp_IL1IT05um_1",
                                 "171009_Smp_IT05um_2" = "171009_Smp_IL1IT05um_2",
                                 "171009_Smp_IT05um_3" = "171009_Smp_IL1IT05um_3",
                                 "171016_Smp_IT0_1" = "171016_Smp_IL2IT0_1",
                                 "171016_Smp_IT0_2" = "171016_Smp_IL2IT0_2",
                                 "171016_Smp_IT0_3" = "171016_Smp_IL2IT0_3",
                                 "171023_Smp_IT05um_1" = "171023_Smp_IL2IT05um_1",
                                 "171023_Smp_IT05um_2" = "171023_Smp_IL2IT05um_2",
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3"))

all.hilics.data <- HILIC_all %>%
  group_by(Mass.Feature) %>%
  mutate(Total.Average = mean(Adjusted.Area))

## Stacked graphs
top.15 <- all.hilics.data %>%
  arrange(desc(Total.Average)) %>%
  select(Mass.Feature) %>%
  unique() %>%
  head(15)
  

stacked.hilics.data <- HILIC_all %>%
  mutate(Full.Total = sum(Adjusted.Area, na.rm = TRUE)) %>%
  filter(Mass.Feature %in% top.15$Mass.Feature) %>%
  #filter(str_detect(Replicate.Name, Treatments)) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four")) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(My.Average = mean(Adjusted.Area, na.rm = TRUE)) %>%
  mutate(Percent.Total = (My.Average / Full.Total)) %>%
  select(Mass.Feature, SampID, My.Average, Percent.Total) %>%
  unique()

# 1_0.2
stacked.hilics.data$SampID <- factor(stacked.hilics.data$SampID, levels = 
                                       c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT", 
                                         "IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT",
                                         "IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um",
                                         "IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))

ggplot(stacked.hilics.data, aes(fill=Mass.Feature, y=My.Average, x=SampID)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, size = 15, vjust = .55),
        legend.text=element_text(size=15)) +
  ggtitle("Top 15 Most Abundant Compounds")


# All HILICS plotted, no filtering
all.hilics.data <- HILIC_all %>%
  group_by(Mass.Feature) %>%
  mutate(Total.Average = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, Total.Average) %>%
  unique()

all.hilics <- ggplot(all.hilics.data, aes(x = reorder(Mass.Feature, -Total.Average), 
                                          y = Total.Average)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(all.hilics)


