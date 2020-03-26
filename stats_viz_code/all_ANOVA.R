library(tidyverse)
options(scipen = 999)

makeLong <- function(df) {
  colnames(df)[1] <- "Replicate.Name"
  df.long <- pivot_longer(df, -Replicate.Name,
               names_to = "Mass.Feature",
               values_to = "Area.BMISd.Normd") %>%
    select(Mass.Feature, Replicate.Name, Area.BMISd.Normd)
  return(df.long)
} 

# BMISd normd
BMISd_1_0.2_wide <- read.csv("data_processed/IsoLagran1_0.2_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)
BMISd_1_5_wide <- read.csv("data_processed/IsoLagran1_5_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)
BMISd_2_0.2_wide <- read.csv("data_processed/IsoLagran2_0.2_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)
BMISd_2_5_wide <- read.csv("data_processed/IsoLagran2_5_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)

BMISd_1_0.2_long <- makeLong(BMISd_1_0.2_wide)
BMISd_1_5_long <- makeLong(BMISd_1_5_wide)
BMISd_2_0.2_long <- makeLong(BMISd_2_0.2_wide)
BMISd_2_5_long <- makeLong(BMISd_2_5_wide)

BMISD.normd <- BMISd_1_0.2_long %>%
  rbind(BMISd_1_5_long) %>%
  rbind(BMISd_2_0.2_long) %>%
  rbind(BMISd_2_5_long)

rm(list=ls(pattern="long|wide"))


# Complete BMISd not normd
BMISd.notnormd <- read.csv("data_processed/BMIS_Output_2020-03-26.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature, Adjusted.Area, Run.Cmpd) %>%
  filter(!str_detect(Run.Cmpd, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m")) %>%
  separate(Run.Cmpd, sep = " ", into = c("Replicate.Name"), remove = TRUE) %>%
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
  select(Mass.Feature, Replicate.Name, SampID, Adjusted.Area) %>%
  drop_na()


grouped.BMISd <- BMISd.notnormd %>%
  left_join(BMISD.normd) %>%
  select(-Replicate.Name) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Avg.Area = mean(Adjusted.Area)) %>%
  mutate(Count = n()) %>%
  filter(!Count < 3) %>%
  ungroup() %>%
  mutate(SampID = factor(SampID, ordered = TRUE)) 

glimpse(grouped.BMISd)
levels(grouped.BMISd$SampID)

ggplot(grouped.BMISd, aes(x = Mass.Feature, y = Avg.Area, fill = SampID)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_blank())

# Apply ANOVA to dataframe, summarize and check significance
AnovaList <- lapply(split(grouped.BMISd, grouped.BMISd$Mass.Feature), function(i) {
  aov(lm(Adjusted.Area ~ SampID, data = i))
})
AnovaListSummary <- lapply(AnovaList, function(i) {
  summary(i)
})
  
# Summarize ANOVA and create dataframe of significance
AnovaDF <- as.data.frame(do.call(rbind, lapply(AnovaListSummary, function(x) {temp <- unlist(x)})))
colnames(AnovaDF)[9] <- "AnovaP"
AnovaDF$AnovaQ <- p.adjust(AnovaDF$AnovaP, method = "fdr")

AnovaDF <- AnovaDF %>%
  rownames_to_column(var = "Mass.Feature") %>%
  mutate(AnovaSig = ifelse(AnovaQ < 0.1, TRUE, FALSE)) %>%
  select(Mass.Feature, AnovaP, AnovaQ, AnovaSig) %>%
  arrange(Mass.Feature)

toPlot <- grouped.BMISd %>%
  left_join(AnovaDF) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(AveSmp = mean(Area.BMISd.Normd, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, AveSmp, AnovaSig) %>%
  unique() %>%
  arrange(Mass.Feature) %>%
  drop_na() %>%
  group_by(Mass.Feature) %>%
  mutate(TotalAve = mean(AveSmp))


a <- ggplot(toPlot, aes(x = Mass.Feature, y = AveSmp, fill = AnovaSig)) +
  geom_point(size = 2, shape = 21) +  
  scale_fill_manual(values = c("grey", "royalblue4")) +
  ggtitle("TEST") +
  theme(plot.title = element_text(size = 15),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9, angle = 90)) +
  theme(legend.position="right") 
a
