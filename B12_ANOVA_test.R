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
#############################################################################
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
  # select(Mass.Feature, SampID, Average.Adjusted.Area) %>%
  select(Mass.Feature, Replicate.Name, SampID, Area.BMISd.Normd) %>%
  unique()
AnovaB12 <- AnovaB12[complete.cases(AnovaB12), ]

# WAnovaB12 <- AnovaB12 %>%
#   pivot_wider(names_from = Mass.Feature,
#               values_from = Adjusted.Area)
###############################################################################

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
Norm.Areas <- ggplot(AnovaB12, aes(x = SampID, y = Area.BMISd.Normd, fill = SampID)) +
  geom_boxplot() +
  facet_wrap(~Mass.Feature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_jitter(shape = 15,
              color = "steelblue",
              position = position_jitter(0.21))
#Norm.Areas

# Apply ANOVA to dataframe and check significance
AnovaList <- lapply(split(AnovaB12, AnovaB12$Mass.Feature), function(i) {
  aov(lm(Area.BMISd.Normd ~ SampID, data = i))
})

AnovaDF <- as.data.frame(do.call(rbind, lapply(AnovaList, function(x) {temp <- unlist(x)})))
colnames(AnovaDF)[9] <- "AnovaP"

AnovaDF <- AnovaDF %>%
  rownames_to_column(var = "Mass.Feature") %>%
  mutate(AnovaSig = ifelse(AnovaP < 0.1, TRUE, FALSE)) %>%
  select(Mass.Feature, AnovaP, AnovaSig) %>%
  arrange(Mass.Feature)


# Apply Tukey HSD to AnovaList and check significance
TukeyList <- lapply(AnovaList, function(x) TukeyHSD(x))
TukeyDF <- as.data.frame(do.call(rbind, lapply(TukeyList, function(x) {temp <- unlist(x)}))) %>%
  select(SampID10:12) %>%
  rownames_to_column("Mass.Feature") %>%
  rename(noB12vControl = SampID10,
         WB12vControl = SampID11,
         WB12vnoB12 = SampID12) %>%
  mutate(noB12vControl_Sig = ifelse(noB12vControl < 0.1, TRUE, FALSE),
          WB12vControl_Sig = ifelse(WB12vControl < 0.1, TRUE, FALSE),
          WB12vnoB12_Sig = ifelse(WB12vnoB12 < 0.1, TRUE, FALSE)) %>%
  #mutate(WvnoB12_FC = log2(rowMeans(TukeyDF[, myTreat1]) / rowMeans(WBMISd[, myTreat2])))
  arrange(Mass.Feature)


# Join with original data to plot
toPlot <- full.BMISd %>%
  left_join(AnovaDF) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(AveSmp = mean(Area.BMISd.Normd, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, AveSmp, AnovaSig) %>%
  unique() %>%
  left_join(TukeyDF) %>%
  arrange(Mass.Feature)


a <- ggplot(toPlot, aes(x = AveSmp, y = AnovaP, fill = AnovaSig,
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
  labs(x="Average peak size", y="AnovaP") +
  theme(legend.position="right") 
  #geom_text(data = toPlot)
a

