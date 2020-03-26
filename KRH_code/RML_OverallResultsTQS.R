#source("src/biostats.R")
library(tidyverse)
#library("vegan")
library(gridExtra) 


# Get wide data with stats
TPTQSdatAll <- read.csv("KRH_info/Filtered_Wide_Combined_wStats.csv")
TPTQSdat <- TPTQSdatAll %>% 
  select(Compound.Name:RB12LL_EF) 


# Make new column that is a combo of BOTH, ONE, or NONE
TPSigs <-TPTQSdatAll %>%
  select(Compound.Name, B12Sig, B12HLSig, B12LLSig, LightSig, LightRB12Sig, LightLB12Sig)
Both <- TPSigs %>%
  filter(B12Sig == "TRUE" | (B12HLSig == "TRUE" & B12LLSig == "TRUE")) %>%
  mutate(B12SigToPlot = "BOTH")
One <- TPSigs %>%
  filter(B12HLSig == "TRUE" | B12LLSig == "TRUE") %>%
  filter(!Compound.Name %in% Both$Compound.Name) %>%
  mutate(B12SigToPlot = "ONE")
Neither <- TPSigs %>%
  filter(!Compound.Name %in% Both$Compound.Name) %>%
  filter(!Compound.Name %in% One$Compound.Name) %>%
  mutate(B12SigToPlot = "NONE")
B12SigToPlot <- rbind(Neither, One, Both)
Both <- TPSigs %>%
  filter(LightSig == "TRUE" | (LightRB12Sig == "TRUE" & LightLB12Sig == "TRUE")) %>%
  mutate(LightSigToPlot = "BOTH")
One <- TPSigs %>%
  filter(LightRB12Sig == "TRUE" | LightLB12Sig == "TRUE") %>%
  filter(!Compound.Name %in% Both$Compound.Name) %>%
  mutate(LightSigToPlot = "ONE")
Neither <- TPSigs %>%
  filter(!Compound.Name %in% Both$Compound.Name) %>%
  filter(!Compound.Name %in% One$Compound.Name) %>%
  mutate(LightSigToPlot = "NONE")
LightSigToPlot <- rbind(Neither, One, Both)
TPTQSdatAll <- TPTQSdatAll %>% left_join(B12SigToPlot) %>% left_join(LightSigToPlot)


#Tp Plot - B12
# Univariate stats for overall changes
a <- ggplot(TPTQSdatAll, aes(x = AveSmp, y = -1*(RvLB12_FC), fill = B12SigToPlot, alpha = B12SigToPlot,
                             text = Compound.Name))+
  geom_point(size = 3, shape = 21, stroke=0)+  
  scale_fill_manual(values =c("royalblue4", "grey", "lightskyblue3")) +
  scale_alpha_manual(values =c(1, 0.5, 0.7))+
  scale_x_log10() +
  ggtitle("Thalassiosira pseudonana") +
  theme(plot.title = element_text(face= "italic", size = 9),
        legend.position="none",
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        axis.text=element_text(size=8))+
  labs(x="Average peak size", y=expression(paste(Log[2], "(Limited/Replete Cobalamin)", sep = ""))) +
  theme(legend.position="none")+
  scale_y_continuous(limits=c(-6, 6))
a


#Tp Plot - Light
# Univariate stats for overall changes
b <- ggplot(TPTQSdatAll, aes(x = AveSmp, y = -1*(HvsLL_FC), fill = LightSigToPlot, alpha = LightSigToPlot, 
                             text = Compound.Name))+
  geom_point(size = 3, shape = 21, stroke=0)+  
  scale_fill_manual(values =c("royalblue4", "grey", "lightskyblue3")) +
  scale_alpha_manual(values =c(1, 0.5, 0.7))+
  scale_x_log10()  +
  theme(plot.title = element_text(face= "italic"),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9))+
  labs(x="Average peak size") +
  labs(x="Average peak size",y=expression(paste(Log[2], "(Low/High light)", sep = "")))+
  theme(legend.position="none") +
  scale_y_continuous(limits=c(-6, 6))
b

require(gridExtra)
grid.arrange(a, b, ncol=2)


