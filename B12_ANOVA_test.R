library(dplyr)
PATH <- "https://raw.githubusercontent.com/guru99-edu/R-Programming/master/poisons.csv"
df <- read.csv(PATH) %>%
  select(-X) %>% 
  mutate(poison = factor(poison, ordered = TRUE))
glimpse(df)

levels(df$poison)


df2 <- df %>%
  group_by(poison) %>%
  summarise(
    count_poison = n(),
    mean_time = mean(time, na.rm = TRUE),
    sd_time = sd(time, na.rm = TRUE)
  )

ggplot(df, aes(x = poison, y = time, fill = poison)) +
  geom_boxplot() +
  geom_jitter(shape = 15,
              color = "steelblue",
              position = position_jitter(0.21)) +
  theme_classic()

anova_one_way <- aov(time~poison, data = df)
summary(anova_one_way)

TukeyHSD(anova_one_way)


####
WAnovaB12 <- AnovaB12 %>%
  pivot_wider(names_from = SampID,
              values_from = Average.Adjusted.Area)

AnovaB12_2 <- AnovaB12 %>%
  group_by(Mass.Feature) %>%
  summarise(
    count_Mass.Feature = n(),
    mean_Area = mean(Average.Adjusted.Area, na.rm = TRUE),
    sd_Area = mean(Average.Adjusted.Area, na.rm = TRUE)
  )

ggplot(AnovaB12, aes(x = SampID, y = Average.Adjusted.Area)) +
  geom_boxplot() +
  geom_jitter(shape = 15,
              color = "steelblue",
              position = position_jitter(0.21)) +
  theme_classic()

anova_one_way <- aov(Average.Adjusted.Area~SampID, data = AnovaB12)
summary(anova_one_way)