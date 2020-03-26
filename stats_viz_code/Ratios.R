library(ggplot2)
library(tidyverse)
## Specific compound ratios

BMISd <- read.csv("data_processed/BMIS_Output_2020-02-20.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature, Adjusted.Area, Run.Cmpd) %>%
  filter(!str_detect(Run.Cmpd, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m")) %>%
  separate(Run.Cmpd, sep = " ", into = c("Replicate.Name"), remove = FALSE) %>%
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
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3")) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), fill = "right", remove = FALSE) %>%
  select(Mass.Feature, SampID, Adjusted.Area) %>%
  drop_na() %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Averages = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, Averages) %>%
  unique() 

#####
Dataset <- BMISd %>%
  filter(str_detect(SampID, "5um")) %>%
  filter(str_detect(SampID, "IL2"))

# Dataset$SampID <- factor(Dataset$SampID,
#                               levels = c("IL1IT05um","IL1Control5um", 
#                                          "IL1DMB5um", "IL1DMBnoBT5um",
#                                          "IL1noBT5um", "IL1WBT5um", 
#                                          "IL1DSW5um"))

# 1_0.2
# Dataset$SampID <- factor(Dataset$SampID,
#                               levels = c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW",
#                                          "IL1DMBnoBT", "IL1noBT"))
# 1_5
# Dataset$SampID <- factor(Dataset$SampID,
#                               levels = c("IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um",
#                                          "IL1DMBnoBT5um", "IL1noBT5um"))
# 2_0.2
# Dataset$SampID <- factor(Dataset$SampID,
#                               levels = c("IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW",
#                                          "IL2DMBnoBT", "IL2noBT"))
# # 2_5
Dataset$SampID <- factor(Dataset$SampID,
                              levels = c("IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um",
                                         "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))
######

# SAM : SAH
SAM.SAH <- Dataset %>%
  filter(Mass.Feature %in% c("Adenosyl Methionine", "Adenosyl Homocysteine")) %>%
  group_by(SampID) %>%
  add_tally() %>%
  filter(!n < 2) %>%
  mutate(Ratios =
          Averages[Mass.Feature == "Adenosyl Methionine"] /
           Averages[Mass.Feature == "Adenosyl Homocysteine"])

Plot2_5 <- ggplot(SAM.SAH, aes(x = SampID, y = Averages, fill = Mass.Feature)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  ggtitle("Anticyclonic Eddy, 5um")
Plot2_5

require(gridExtra)
grid.arrange(Plot1_0.2, Plot1_5, Plot2_0.2, Plot2_0.2, nrow = 2)

library(ggpubr)
ggarrange(Plot1_0.2, Plot1_5, Plot2_0.2, Plot2_5,
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


# Glutamine : Glutamate
glutamine.glutamate <- BMISd %>%
  filter(Mass.Feature %in% c("Glutamine"))

# Glutathione : GSSG
glutathione.gssg <- Dataset %>%
  filter(Mass.Feature %in% c("Glutathione Disulfide"))


