library(tidyverse)
library(gridExtra) 

# Get wide data with stats
mydatAll <- read.csv("data_processed/WBMISd_wStats.csv")

# Make new column that is a combo of BOTH, ONE, or NONE
allSigs <-mydatAll %>%
  select(Mass.Feature, contains("Sig"))

Both <- allSigs %>%
  filter(B12Sig == "TRUE" & DMBSig == "TRUE") %>%
  mutate(B12SigToPlot = "BOTH")

One <- allSigs %>%
  filter(B12Sig == "TRUE" | DMBSig == "TRUE") %>%
  filter(!Mass.Feature%in% Both$Compound.Name) %>%
  mutate(B12SigToPlot = "ONE")

Neither <- allSigs %>%
  filter(!Mass.Feature %in% Both$Mass.Feature) %>%
  filter(!Mass.Feature %in% One$Mass.Feature) %>%
  mutate(B12SigToPlot = "NONE")

B12SigToPlot <- rbind(Neither, One, Both)

Both <- allSigs %>%
  filter(T0vTFSig == "TRUE" & T0vDSWSig) %>%
  mutate(T0vTFSigToPlot = "BOTH")

One <- allSigs %>%
  filter(T0vTFSig == "TRUE" | T0vDSWSig == "TRUE") %>%
  filter(!Mass.Feature %in% Both$Mass.Feature) %>%
  mutate(T0vTFSigToPlot = "ONE")

Neither <- allSigs %>%
  filter(!Mass.Feature %in% Both$Mass.Feature) %>%
  filter(!Mass.Feature %in% One$Mass.Feature) %>%
  mutate(T0vTFSigToPlot = "NONE")

T0vTFSigToPlot <- rbind(Neither, One, Both)

mydatAll <- mydatAll %>% left_join(B12SigToPlot) %>% left_join(T0vTFSigToPlot)


#Tp Plot - B12
# Univariate stats for overall changes
a <- ggplot(mydatAll, aes(x = AveSmp, y = -1*(WvnoB12_FC), fill = B12SigToPlot, alpha = B12SigToPlot,
                             label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +  
  scale_fill_manual(values =c("royalblue4", "grey", "lightskyblue3")) +
  scale_alpha_manual(values =c(1, 0.5, 0.7))+
  scale_x_log10() +
  ggtitle("TEST, B12") +
  theme(plot.title = element_text(face= "italic", size = 9),
        legend.position="none",
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        axis.text=element_text(size=8)) +
  labs(x="Average peak size", y=expression(paste(Log[2], "(Limited/Replete Cobalamin)", sep = ""))) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(-6, 6)) +
  geom_text(hjust = 0, nudge_x = 0.05)
a


#Tp Plot - Light
# Univariate stats for overall changes
b <- ggplot(mydatAll, aes(x = AveSmp, y = -1*(T0vTF_FC), fill = T0vTFSigToPlot, alpha = T0vTFSigToPlot, 
                             label = Mass.Feature)) +
  geom_point(size = 3, shape = 21, stroke=0) +  
  scale_fill_manual(values =c("royalblue4", "grey", "lightskyblue3")) +
  scale_alpha_manual(values =c(1, 0.5, 0.7)) +
  scale_x_log10() +
  ggtitle("TEST, T0vTF") +
  theme(plot.title = element_text(face= "italic"),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9)) +
  labs(x="Average peak size") +
  labs(x="Average peak size",y=expression(paste(Log[2], "(Low/High light)", sep = "")))+
  theme(legend.position="none") +
  scale_y_continuous(limits=c(-6, 6)) +
  geom_text(hjust = 0, nudge_x = 0.05)
b

require(gridExtra)
grid.arrange(a, b, ncol=2)


