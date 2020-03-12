library(tidyverse)
library(gridExtra) 

# Get wide data with stats
mydatAll <- read.csv("data_processed/WBMISd_wStats.csv")

allSigs <- mydatAll %>%
  select(Mass.Feature, contains("Sig"))

# B12
YestoSig <- allSigs %>%
  filter(B12vnoB12Sig == "TRUE") %>%
  mutate(B12SigToPlot = "YES")

NotoSig <- allSigs %>%
  filter(!Mass.Feature %in% YestoSig$Mass.Feature) %>%
  mutate(B12SigToPlot = "NO")

B12SigToPlot <- rbind(YestoSig, NotoSig)

# DMB
YestoSig <- allSigs %>%
  filter(DMBvnoDMBSig == "TRUE") %>%
  mutate(DMBvnoDMBSig = "YES")

NotoSig <- allSigs %>%
  filter(!Mass.Feature %in% YestoSig$Mass.Feature) %>%
  mutate(DMBSigToPlot = "NO")

DMBSigToPlot <- rbind(YestoSig, NotoSig)

# T0 v TFinal
YestoSig <- allSigs %>%
  filter(T0vTFSig == "TRUE") %>%
  mutate(T0vTFSigToPlot = "YES")

NotoSig <- allSigs %>%
  filter(!Mass.Feature %in% YestoSig$Mass.Feature) %>%
  mutate(T0vTFSigToPlot = "NO")

T0vTFSigToPlot <- rbind(YestoSig, NotoSig)

# T0 v Deep Sea Water
YestoSig <- allSigs %>%
  filter(T0vDSWSig == "TRUE") %>%
  mutate(T0vDSWSigToPlot = "YES")

NotoSig <- allSigs %>%
  filter(!Mass.Feature %in% YestoSig$Mass.Feature) %>%
  mutate(T0vDSWSigToPlot = "NO")

T0vDSWSigToPlot <- rbind(YestoSig, NotoSig)

mydatAll <- mydatAll %>% left_join(B12SigToPlot) %>% left_join(DMBSigToPlot) %>% left_join(T0vTFSigToPlot) %>% left_join(T0vDSWSigToPlot)


# B12 vs no B12 significance
a <- ggplot(mydatAll, aes(x = AveSmp, y = -1*(WvnoB12_FC), fill = B12SigToPlot, #alpha = B12SigToPlot
                             label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +  
  scale_fill_manual(values =c("grey", "royalblue4", "lightskyblue3")) +
  scale_alpha_manual(values =c(1, 0.5, 0.7)) +
  scale_x_log10() +
  ggtitle("B12 vs No B12 Significance") +
  theme(plot.title = element_text(size = 15),
        legend.position="left",
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        axis.text=element_text(size=8)) +
  labs(x="Average peak size", y=expression(paste(Log[2], "With/Without Cobalamin", sep = ""))) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(-3, 3)) +
  geom_text(data = subset(mydatAll, B12SigToPlot == "YES" | AveSmp > 100000000 | (-1*(WvnoB12_FC) > 2) | (-1*(WvnoB12_FC) < -2)), 
                          hjust = "inward", nudge_x = 0.05)
a


# Time 0 vs DSW significance
b <- ggplot(mydatAll, aes(x = AveSmp, y = -1*(T0vDSW_FC), fill = T0vDSWSigToPlot, # alpha = T0vDSWSigToPlot, 
                             label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +  
  scale_fill_manual(values = c("grey", "royalblue4")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_x_log10() +
  ggtitle("Time 0 vs DeepSeaWater Significance") +
  theme(plot.title = element_text(size = 15),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9)) +
  labs(x="Average peak size") +
  labs(x="Average peak size",y=expression(paste(Log[2], "T0/DSW", sep = ""))) +
  theme(legend.position="right") +
  scale_y_continuous(limits=c(-6, 6)) +
  geom_text(data = subset(mydatAll, T0vDSWSigToPlot == "YES" | AveSmp > 100000000 | AveSmp < 100000 | (-1*(T0vDSW_FC) > 3)),
            hjust = "inward", nudge_x = 0.05) 
b


require(gridExtra)
grid.arrange(a, b, ncol=2)


