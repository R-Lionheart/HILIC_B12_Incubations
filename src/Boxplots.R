library(gridExtra)
library(patchwork)
library(tidyverse)

options(scipen = 999)
options(digits = 5)
source("src/Functions.R")

replace_nonvalues <- function(x) (gsub(NaN, NA, x))

# Munge full Incubations dataset
Complete.Dataset <- read.csv("data_processed/MSDial_QE_QC_Output_HILIC.csv",
                        stringsAsFactors = FALSE) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m|Process|Std")) %>%
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
  mutate(Area.with.QC = ifelse(str_detect(all.Flags, "Blank.Flag"), 
                              NA, Area.with.QC)) %>%
  mutate(Area.with.QC = ifelse(str_detect(all.Flags, "Blank.Flag"), 
                               NA, Area.with.QC)) %>%
  select(Replicate.Name, Metabolite.Name, Area.with.QC) %>%
  mutate(Size.Fraction = ifelse(str_detect(Replicate.Name, "5um"), "Large.Filter", "Small.Filter"),
         Eddy = ifelse(str_detect(Replicate.Name, "IL1"), "Cyclonic", "Anticyclonic")) %>%
  group_by(Metabolite.Name, Eddy, Size.Fraction) %>%
  mutate(Binned.Group = ifelse(str_detect(Replicate.Name, "DSW"), "DeepSeaWater",
                               ifelse(str_detect(Replicate.Name, "T0"), "TimeZero",
                                      ifelse(str_detect(Replicate.Name, "Control"), "Control",
                                             ifelse(str_detect(Replicate.Name, "DMBnoBT|noBT"), "Nutrients", "NoNutrients"))))) %>%
  mutate(Binned.Group = ifelse(str_detect(Replicate.Name, "5um"), paste(Binned.Group, "LargeFilter", sep = "_"),
                               paste(Binned.Group, "SmallFilter", sep = "_"))) %>%
  mutate(Binned.Group = ifelse(Eddy == "Cyclonic", paste(Binned.Group, "Cyclonic", sep = "_"), 
                               paste(Binned.Group, "Anticyclonic", sep = "_"))) %>%
  ungroup() %>%
  select(Metabolite.Name, Area.with.QC, Binned.Group) %>%
  mutate(Binned.Group = factor(Binned.Group, ordered = TRUE)) %>%
  group_by(Metabolite.Name) %>%
  mutate(CountVals = sum(!is.na(Area.with.QC))) %>%
  filter(CountVals != 0) %>%
  ungroup() %>%
  separate(Binned.Group, c("SampID", "B", "C")) %>%
  unite("Grouping.ID", B:C)


PlotAllCompounds <- function (df) {
  all.plots <- ggplot(df, aes(x = SampID, y = Area.with.QC, fill = SampID)) +
    facet_wrap(~ Grouping.ID, scales = "free") +
    geom_boxplot() +
    scale_fill_grey() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(paste(unique(df[[1]])))
  
  return(all.plots)
}

to.plot <- Complete.Dataset %>%
  group_by(Metabolite.Name) %>%
  group_split()

plotlist <- lapply(to.plot, PlotAllCompounds)

pdf("~/Downloads/B12Incubations_HILIC_AllCompounds.pdf")
lapply(plotlist, print)
dev.off()

