library(tidyverse)
options(scipen = 999)

####
BMISd <- read.csv("data_processed/IsoLagran1_0.2_notnormed.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)
WBMISd.normed <- read.csv("data_processed/IsoLagran2_5_normed.csv", stringsAsFactors = FALSE)
BMISd.normed <- WBMISd.normed %>%
  pivot_longer(-X,
               names_to = "Mass.Feature",
               values_to = "Area.BMISd.Normd") %>%
  rename(Replicate.Name = X) %>%
  select(Mass.Feature, Replicate.Name, Area.BMISd.Normd)


WBMISd <- BMISd %>%
  spread(key = "Replicate.Name", value = "Adjusted.Area")
WBMISd <- WBMISd[complete.cases(WBMISd), ]
mySamps <- colnames(WBMISd)

AnovaB12 <- BMISd.normed %>%
  filter(str_detect(Replicate.Name, "WBT|IL2noBT|Control")) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), remove = FALSE) %>%
  select(Mass.Feature, Replicate.Name, SampID, Area.BMISd.Normd) %>%
  # group_by(Mass.Feature, SampID) %>%
  # mutate(Average.Adjusted.Area = mean(Adjusted.Area, na.rm = TRUE)) %>%
  #select(Mass.Feature, SampID, Average.Adjusted.Area) %>%
  select(Mass.Feature, Replicate.Name, SampID, Area.BMISd.Normd) %>%
  unique()
AnovaB12 <- AnovaB12[complete.cases(AnovaB12), ]

WAnovaB12 <- AnovaB12 %>%
  pivot_wider(names_from = Mass.Feature,
              values_from = Adjusted.Area)


## my version of online test
test_df <- AnovaB12 %>%
  select(-Replicate.Name) %>%
  #
  filter(str_detect("Adenine", Mass.Feature)) %>%
  #
  mutate(SampID = factor(SampID, ordered = TRUE))
glimpse(test_df)
levels(test_df$SampID)

summary_test_df <- test_df %>%
  group_by(SampID, Mass.Feature) %>%
  summarise(
    count_SampID = n(),
    mean_area = mean(Area.BMISd.Normd, na.rm = TRUE),
    sd_area = sd(Area.BMISd.Normd, na.rm = TRUE)
  )

ggplot(test_df, aes(x = SampID, y = Area.BMISd.Normd, fill = SampID)) +
  geom_boxplot() +
  facet_wrap(~Mass.Feature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_jitter(shape = 15,
              color = "steelblue",
              position = position_jitter(0.21))

test_anova_one_way <- aov(Area.BMISd.Normd~SampID, data = test_df)
summary(test_anova_one_way)
summary(test_anova_one_way)[[1]][["Pr(>F)"]]
TukeyHSD(test_anova_one_way)

# two way isn't useful because mass feature isn't grouping?
# test_anova_two_way <- aov(Area.BMISd.Normd~SampID + Mass.Feature, data = test_df)
# summary(test_anova_two_way)





