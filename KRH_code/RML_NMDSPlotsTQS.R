source("src/biostats.R")
library(tidyverse)
library("vegan")

# Get wide data with stats

TPTQSdatAll <- read.csv("KRH_version/Filtered_Wide_Combined_wStats.csv")
TPTQSdat <- TPTQSdatAll %>% 
  select(LB12HL_AB:RB12LL_EF) %>% 
  as.data.frame()

row.names(TPTQSdat) <- TPTQSdatAll$Compound.Name
TPTQSdat <- t(TPTQSdat)
TPTreatDat <- read.csv("KRH_version/Treatment_data.csv", header=TRUE, row.names=1)


#Standarize the data
TPExpDat.std <- data.stand(TPTQSdat, method  = 'standardize', margin = 'column', plot = F)

tp.nmds <- metaMDS(TPExpDat.std, distance='euclidean', k=2, autotransform=FALSE, wascores = FALSE, noshare = FALSE, trymax=10000)
pointlocation <- tp.nmds[['points']] %>% as.data.frame() %>% cbind(TPTreatDat)

g<- ggplot(data = pointlocation, aes(x =MDS1, y =  MDS2, fill = Light.Status, 
                                     shape = Cobalamin.Status, group = Supergroup))+
  geom_polygon(fill = NA, color = "black") +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(22,21)) +
  scale_fill_manual(values = c("grey","black")) +   
  ggtitle("Thalassiosira pseudonana") +
  theme(plot.title = element_text(face= "italic", size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text( size = 8),
        axis.text = element_text(size = 8))+
  labs(y="Axis 2") +
  theme(legend.position="none")
g


#Get and write out vector info, write out CSV into the correct folder.-----
vec.info<-envfit(tp.nmds$points, TPExpDat.std, perm=9999)
names <- as.data.frame(vec.info[[1]][[1]])%>%
  mutate(Compound.Name = row.names(.),
         pvalues = vec.info[[1]][[4]]) %>%
  select(Compound.Name, pvalues)

#Get and write out ANOSIM data info
TPLight.anosim <- anosim(TPExpDat.std, TPTreatDat[,4], distance = 'euclidean')
TPLightanosimResults <- c("Tp", "Light", TPLight.anosim$statistic, TPLight.anosim$signif)
TPB12.anosim <- anosim(TPExpDat.std,TPTreatDat[,3], distance = 'euclidean')
TPB12anosimResults <- c("Tp", "Cobalamin", TPB12.anosim$statistic, TPB12.anosim$signif)
AnosimResults <- data.frame(x = TPLightanosimResults,  y = TPB12anosimResults, z= c("Org", "Variable", "ANOSIMStat", "ANOSIMpvalue"))
