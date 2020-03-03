library(tidyverse)
options(scipen = 999)

complete.BMISd <- read.csv("data_processed/BMIS_Output_2020-02-20.csv", stringsAsFactors = FALSE) %>%
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


grouped.BMISd <- complete.BMISd %>%
  select(-Replicate.Name) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Avg.Area = mean(Adjusted.Area)) %>%
  mutate(Count = n()) %>%
  filter(!Count < 3) %>%
  ungroup() %>%
  mutate(SampID = factor(SampID, ordered = TRUE)) 

glimpse(grouped.BMISd)
levels(grouped.BMISd$SampID)


# Apply ANOVA to dataframe, summarize and check significance
AnovaList <- lapply(split(grouped.BMISd, grouped.BMISd$Mass.Feature), function(i) {
  aov(lm(Adjusted.Area ~ SampID, data = i))
}) 
AnovaListSummary <- lapply(AnovaList, function(i) {
  summary(i)
})


# Apply Tukey to dataframe and check significance
#TukeyList <- lapply(AnovaList, function(x) TukeyHSD(x))

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
  # group_by(Mass.Feature, SampID) %>%
  # mutate(AveSmp = mean(Area.BMISd.Normd, na.rm = TRUE)) %>%
  # select(Mass.Feature, SampID, AveSmp, AnovaSig) %>%
  # unique() %>%
  # left_join(TukeyDF) %>%
  # left_join(Analysis) %>%
  # arrange(Mass.Feature) %>%
  # drop_na() %>%
  # group_by(Mass.Feature) %>%
  # mutate(TotalAve = mean(AveSmp))


a <- ggplot(toPlot, aes(x = Mass.Feature, y = Avg.Area, fill = AnovaSig)) +
  geom_point(size = 3, shape = 21) +  
  scale_fill_manual(values = c("grey", "royalblue4")) +
  ggtitle("TEST") +
  theme(plot.title = element_text(size = 15),
        legend.position="none",
        axis.title.y=element_text(size=9, angle = 90),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9)) +
  theme(legend.position="right") 
a
