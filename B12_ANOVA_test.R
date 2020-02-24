library(tidyverse)
options(scipen = 999)

# Import required files
BMISd.long <- read.csv("data_processed/IsoLagran2_5_notnormd.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature:Adjusted.Area)
BMISd.wide <- read.csv("data_processed/IsoLagran2_5_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)
colnames(BMISd.wide)[1] <- "Replicate.Name"
BMISd.wide.notnormd <- BMISd.long %>%
  pivot_wider(names_from = Replicate.Name, 
              values_from = Adjusted.Area)

# Change all sets to long format
BMISd.long.normd <- BMISd.wide %>%
  pivot_longer(-Replicate.Name,
    names_to = "Mass.Feature",
    values_to = "Area.BMISd.Normd") %>%
  select(Mass.Feature, Replicate.Name, Area.BMISd.Normd)

BMISd.wide.normd <- BMISd.long.normd %>%
  spread(key = "Replicate.Name", value = "Area.BMISd.Normd")

# Combine all data and rearrange
full.BMISd <- BMISd.long.normd %>%
  left_join(BMISd.long) %>%
  filter(str_detect(Replicate.Name, "WBT|IL2noBT|IL2Control")) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), remove = FALSE) %>%
  select(Mass.Feature, Replicate.Name, SampID, Area.BMISd.Normd, Adjusted.Area) %>%
  unique()
full.BMISd <- full.BMISd[complete.cases(full.BMISd), ]



# KRH analysis ------------------------------------------------------------
Analysis <- BMISd.wide.notnormd[complete.cases(BMISd.wide.notnormd), ]
mySamps <- colnames(Analysis)

# Add B12 vs noB12 stats
myTreat1 <- mySamps[grepl("IL2WBT", mySamps)]
myTreat2 <- mySamps[grepl("IL2noBT", mySamps)]
myTreat3 <- mySamps[grepl("IL2Control", mySamps)]
myTreatsdf <- Analysis[, c(myTreat1, myTreat2, myTreat3)]

Analysis <- Analysis %>%
  mutate(WvnoB12_FC = log2(rowMeans(Analysis[, myTreat1]) / rowMeans(Analysis[, myTreat2]))) %>%
  mutate(noB12vConv_FC = log2(rowMeans(Analysis[, myTreat2]) / rowMeans(Analysis[, myTreat3]))) %>%
  mutate(WB12vConv_FC = log2(rowMeans(Analysis[, myTreat1]) / rowMeans(Analysis[, myTreat3]))) %>%
  select(matches('Mass|FC'))
  
  
# Set up data for ANOVA
AnovaB12 <- full.BMISd %>%
  select(-Replicate.Name) %>%
  mutate(SampID = factor(SampID, ordered = TRUE)) %>%
  arrange(Mass.Feature) %>%
  ##
  filter(!Mass.Feature == "Hydroxylysine")
  ##
glimpse(AnovaB12) #Use for the future
levels(AnovaB12$SampID)

# Graph normalized areas for reference
Normd.Areas <- ggplot(AnovaB12, aes(x = SampID, y = Area.BMISd.Normd, fill = SampID)) +
  geom_boxplot() +
  facet_wrap(~Mass.Feature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_jitter(shape = 15,
              color = "steelblue",
              position = position_jitter(0.21))
#Normd.Areas

# Apply ANOVA to dataframe, summarize and check significance
AnovaList <- lapply(split(AnovaB12, AnovaB12$Mass.Feature), function(i) {
  aov(lm(Area.BMISd.Normd ~ SampID, data = i))
}) 
AnovaListSummary <- lapply(AnovaList, function(i) {
  summary(i)
})

# Apply Tukey to dataframe and check significance
TukeyList <- lapply(AnovaList, function(x) TukeyHSD(x))

# Summarize ANOVA and create dataframe of significance
AnovaDF <- as.data.frame(do.call(rbind, lapply(AnovaListSummary, function(x) {temp <- unlist(x)})))
colnames(AnovaDF)[9] <- "AnovaP"

AnovaDF <- AnovaDF %>%
  rownames_to_column(var = "Mass.Feature") %>%
  mutate(AnovaSig = ifelse(AnovaP < 0.1, TRUE, FALSE)) %>%
  select(Mass.Feature, AnovaP, AnovaSig) %>%
  arrange(Mass.Feature)

# Summarize Tukey HSD and create dataframe of significance
TukeyDF <- as.data.frame(do.call(rbind, lapply(TukeyList, function(x) {temp <- unlist(x)}))) %>%
  select(SampID10:12) %>%
  rownames_to_column("Mass.Feature") %>%
  rename(noB12vControl = SampID10,
         WB12vControl = SampID11,
         WB12vnoB12 = SampID12) %>%
  mutate(noB12vControl_Sig = ifelse(noB12vControl < 0.1, TRUE, FALSE),
          WB12vControl_Sig = ifelse(WB12vControl < 0.1, TRUE, FALSE),
          WB12vnoB12_Sig = ifelse(WB12vnoB12 < 0.1, TRUE, FALSE)) %>%
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

#mutate(WvnoB12_FC = log2(rowMeans(TukeyDF[, myTreat1]) / rowMeans(WBMISd[, myTreat2])))

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

