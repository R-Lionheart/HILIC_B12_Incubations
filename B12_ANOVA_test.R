library(tidyverse)
options(scipen = 999)

# Import required files
BMISd <- read.csv("data_processed/IsoLagran2_5_notnormed.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)
WBMISd.normed <- read.csv("data_processed/IsoLagran2_5_normed.csv", stringsAsFactors = FALSE, check.names = FALSE)
colnames(WBMISd.normed)[1] <- "Replicate.Name"

# Change all sets to long format
BMISd.normed <- WBMISd.normed %>%
  pivot_longer(-Replicate.Name,
    names_to = "Mass.Feature",
    values_to = "Area.BMISd.Normd") %>%
  select(Mass.Feature, Replicate.Name, Area.BMISd.Normd)


# Combine all data and rearrange
full.BMISd <- BMISd.normed %>%
  left_join(BMISd) %>%
  filter(str_detect(Replicate.Name, "WBT|IL2noBT|IL2Control")) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), remove = FALSE) %>%
  select(Mass.Feature, Replicate.Name, SampID, Area.BMISd.Normd, Adjusted.Area) %>%
  unique()
full.BMISd <- AnovaB12[complete.cases(AnovaB12), ]


# Set up data for ANOVA
AnovaB12 <- full.BMISd %>%
  select(-Replicate.Name) %>%
  mutate(SampID = factor(SampID, ordered = TRUE)) %>%
  arrange(Mass.Feature)
glimpse(AnovaB12) #Use for the future
levels(AnovaB12$SampID)

# Graph normalized areas for reference
ggplot(AnovaB12, aes(x = SampID, y = Area.BMISd.Normd, fill = SampID)) +
  geom_boxplot() +
  facet_wrap(~Mass.Feature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_jitter(shape = 15,
              color = "steelblue",
              position = position_jitter(0.21))

# Split and apply ANOVA function, result is a named list:
AnovaList <- lapply(split(AnovaB12, AnovaB12$Mass.Feature), function(i) {
  anova(lm(Area.BMISd.Normd ~ SampID, data = i))
})

# Turn list to dataframe
AnovaDF <- do.call(rbind, lapply(AnovaList, function(x) {temp <- unlist(x)}))
AnovaDF <- as.data.frame(AnovaDF)
colnames(AnovaDF)[9] <- "Pvalue"

AnovaDF <- AnovaDF %>%
  rownames_to_column(var = "Mass.Feature") %>%
  select(Mass.Feature, Pvalue) %>%
  mutate(Significant = ifelse(Pvalue < 0.1, TRUE, FALSE))

# Join with original data to plot
toPlot <- full.BMISd %>%
  left_join(AnovaDF) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(AveSmp = mean(Area.BMISd.Normd, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, AveSmp, Pvalue, Significant) %>%
  unique()

a <- ggplot(toPlot, aes(x = AveSmp, y = Pvalue, fill = Significant,
                          label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +  
  scale_fill_manual(values = c("grey", "royalblue4")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  ggtitle("TEST") +
  theme(plot.title = element_text(size = 15),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9)) +
  labs(x="Average peak size") +
  labs(x="Average peak size", y=expression(paste(Log[2], "T0/DSW", sep = ""))) +
  theme(legend.position="right") 
  #geom_text(data = toPlot)
a



test_anova_one_way <- aov(Area.BMISd.Normd~SampID, data = test_df)
summary(test_anova_one_way)
summary(test_anova_one_way)[[1]][["Pr(>F)"]]
test <- TukeyHSD(test_anova_one_way)

# two way isn't useful because mass feature isn't grouping?
# test_anova_two_way <- aov(Area.BMISd.Normd~SampID + Mass.Feature, data = test_df)
# summary(test_anova_two_way)





