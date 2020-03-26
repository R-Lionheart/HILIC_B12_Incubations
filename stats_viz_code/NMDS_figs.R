source("src/B12_Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(tidyverse) 
library(vegan)

# User data
pattern = "BMIS_Output"
percentMissing = 0.5

# Functions
makeWide <- function(df) {
  df.wide <- df %>%
    ungroup() %>%
    tidyr::spread(Replicate.Name, Adjusted.Area) %>%
    as.data.frame()
  
    df.rownames <- df.wide[,-1]
    rownames(df.rownames) <- df.wide[,1]
    
    df.rownames[is.na(df.rownames)] <- NA
    
    #df.noNA <- na.omit(df.rownames)
  
  return(df.rownames)
}
RemoveCsv <- function(full.filepaths) {
  # Remove a .csv file extension and obtain basename from a given list of filepaths.
  #
  # Args
  #   Character strings of filepaths in a directory.
  #
  # Returns
  #   Character strings of file basenames, without a csv extension.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

# BMIS import and filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

HILIC_all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  select(Mass.Feature, runDate:replicate, Adjusted.Area) %>%
  unite(Replicate.Name, runDate:replicate, sep = "_") %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m")) %>%
  filter(!Mass.Feature == "Inj_vol") %>%
  filter(!str_detect(Mass.Feature, ","))

### ADJUST FOR IT0 ISSUES ###
HILIC_fixed <- HILIC_all %>%
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
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3"))


# All HILICS plotted, no filtering
all.hilics.data <- HILIC_fixed %>%
  group_by(Mass.Feature) %>%
  mutate(Total.Average = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, Total.Average) %>%
  unique()
all.hilics <- ggplot(all.hilics.data, aes(x = reorder(Mass.Feature, -Total.Average), 
                                          y = Total.Average)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(all.hilics)

# Filter out compounds that are missing too many peaks --------------------------------------------
HILIC_filtered <- HILIC_fixed %>%
  group_by(Mass.Feature) %>%
  mutate(Missing = sum(is.na(Adjusted.Area))) %>%
  mutate(MF.Count = n()) %>%
  filter(!Missing > (percentMissing*MF.Count)) %>%
  select(-c("Missing", "MF.Count"))

# Separate dataset into groups for analysis -------------------------------------------------------
IsoLagran1 <- HILIC_filtered %>%
  filter(str_detect(Replicate.Name, "IL1"))
IsoLagran2 <- HILIC_filtered %>%
  filter(str_detect(Replicate.Name, "IL2"))
IsoLagran1_0.2 <- IsoLagran1 %>%
  filter(!str_detect(Replicate.Name, "5um"))
write.csv(IsoLagran1_0.2, "data_processed/IsoLagran1_0.2_notnormd.csv")
IsoLagran1_5 <- IsoLagran1 %>%
  filter(str_detect(Replicate.Name, "5um"))
write.csv(IsoLagran1_5, "data_processed/IsoLagran1_5_notnormd.csv")
IsoLagran2_0.2 <- IsoLagran2 %>%
  filter(!str_detect(Replicate.Name, "5um"))
write.csv(IsoLagran2_0.2, "data_processed/IsoLagran2_0.2_notnormd.csv")
IsoLagran2_5 <- IsoLagran2 %>%
  filter(str_detect(Replicate.Name, "5um"))
write.csv(IsoLagran2_5, "data_processed/IsoLagran2_5_notnormd.csv")

rm(list = c("IsoLagran1", "IsoLagran2", "HILIC_fixed"))

# Treatment data info --------------------------------------------------------------------------
# Set eddy and size fraction
Dataset <- IsoLagran2_5
EddySize <- "2, Anticyclonic_5um"

Treatment <- Dataset %>%
  ungroup() %>%
  select(Replicate.Name) %>%
  unique() %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  mutate(Control.Status = ifelse(str_detect(Supergroup, "IT0"),
                                 "Incubation", ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater", 
                                                      ifelse(str_detect(Supergroup, "Control"), "Control", "Treatments")))) %>%
  mutate(Treatment.Status = ifelse(Control.Status == "Control", "Control",
                                   ifelse(Control.Status == "DeepSeaWater", "DeepSeaWater",
                                          ifelse(Control.Status == "Incubation", "TimeZero",
                                                 ifelse(str_detect(Supergroup, "DMBnoBT"), "DMBnoB12",
                                                        ifelse(str_detect(Supergroup, "WBT"), "B12",
                                                               ifelse(str_detect(Supergroup, "DMB"), "DMB", "noB12"))))))) %>%
  select(Replicate.Name, Control.Status, Treatment.Status, Supergroup)

# Transform to wide for analysis -----------------------------------------------------------------
Iso_wide <- makeWide(Dataset)
Iso_wide[is.na(Iso_wide)] <- 1000
Iso_wideT <- t(Iso_wide)

# Standardize + NMDS transformations --------------------------------------------------------------
Iso_wide_normalizedT <- decostand(Iso_wideT, method = "standardize", na.rm = TRUE) 
write.csv(Iso_wide_normalizedT, paste("data_processed/IsoLagran", EddySize, "_normd.csv", sep = ""))

# Visualize scree plot of potential ordination axes
dimcheckMDS(Iso_wide_normalizedT, distance="euclidean", k=10, autotransform=FALSE, trymax=20)


Iso_wide_nmds <- vegan::metaMDS(Iso_wide_normalizedT, distance = "euclidean", 
                                k = 3, autotransform = FALSE, trymax = 100)

# Check stressplots, scree diagrams --------------------------------------------------------------

#stressplot(Iso_wide_nmds, main = paste("Stressplot, Eddy", EddySize, sep = " "))
Iso_wide_nmds$stress # Add a flag if this is high?
#nmds.monte(Iso_wide_normalizedT, distance="euclidean", k=3, autotransform=FALSE, trymax=20)

# Experiments with species scores
envfit(Iso_wide_nmds, Iso_wide_normalizedT) 
scores(Iso_wide_nmds)



plot(Iso_wide_nmds$points,type="n")
text(Iso_wide_nmds,labels=row.names(Iso_wideT))

# isodf <- as.data.frame(Iso_wideT)
# plot(Iso_wide_nmds$points,type="n")
# points(Iso_wide_nmds,cex=Iso_wide$`171023_Smp_IL2Control5um_1`)

Iso_pointlocation <- Iso_wide_nmds[['points']] %>% as.data.frame() %>% cbind(Treatment)

Iso_pointlocation$Treatment.Status <- factor(Iso_pointlocation$Treatment.Status,
                              levels = c("TimeZero", "Control", "DMB", "B12", "DeepSeaWater",
                                         "DMBnoB12", "noB12"))

##############################
# Quick vectors
mysol <- vegan::metaMDS(Iso_wide_normalizedT, distance = "euclidean", 
                        k = 2, autotransform = FALSE, trymax = 100)

myNMDS = data.frame(MDS1 = mysol$points[,1], MDS2 = mysol$points[,2])

myvec.sp<-envfit(mysol$points, Iso_wide_normalizedT, perm=1000)
myvec.sp.df<-as.data.frame(myvec.sp$vectors$arrows*sqrt(myvec.sp$vectors$r))
myvec.sp.df$species<-rownames(myvec.sp.df)

ggplot(data = myNMDS, aes(MDS1, MDS2)) + 
  #geom_point(aes(data = MyMeta, color = MyMeta$amt)) +
  geom_segment(data=myvec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour="grey", inherit_aes=FALSE) + 
  geom_text(data=myvec.sp.df,aes(x=MDS1,y=MDS2,label=species),size=3) + 
  ggtitle(paste("Incubation Experiments: Eddy", EddySize, sep = " "))
  
ordiplot(mysol, choices = c(1, 2), type="text", display="sites",
         xlab="Axis1", ylab="Axis2")
plot(myvec.sp, p.max=.01, col="blue")


mysol.d <- vegdist(Iso_wide_normalizedT, "euclidean")
mysitecl.ward<-hclust(mysol.d,method='ward.D')
mysitecl.class<-cutree(mysitecl.ward, k=4)
groups<-levels(factor(mysitecl.class))
mysite.sc <- scores(mysol)

my.p <- ordiplot(mysite.sc, type="n", main="NMDS combined with clustering")
for (i in 1:length(groups))
{
  points(mysite.sc[mysitecl.class==i,], pch=(14+i), cex=2, col=i+1)
}
text(mysite.sc, row.names(Iso_wideT), pos=4, cex=0.7)
ordicluster(my.p, mysitecl.ward, col="dark grey")

legend("bottomleft", paste("Group", c(1:length(groups))),
       pch=14+c(1:length(groups)), col=1+c(1:length(groups)), pt.cex=2)


ordiplot(Iso_wide_nmds, type="n")
ordihull(Iso_wide_nmds, groups=Treatment$Treatment.Status, draw="polygon",col="grey90",label=F)
orditorp(Iso_wide_nmds, display="sites",
         col=c(rep("red",3), rep("orange",3), rep("green"), 3),
         air=0.01,cex=1.25)


colors=c(rep("red",5),rep("blue",5))
ordiplot(Iso_wide_nmds,type="n")
#Plot convex hulls with colors based on treatment
for(i in unique(Treatment$Treatment.Status)) {
  ordihull(Iso_wide_nmds$point[grep(i,Treatment$Treatment.Status),],draw="polygon",
           groups=Treatment$Treatment.Status[Treatment$Treatment.Status==i],col=colors[grep(i,Treatment$Treatment.Status)],label=F) } 
orditorp(Iso_wide_nmds,display="species",col="red",air=0.01)
orditorp(Iso_wide_nmds,display="sites",col=c(rep("green",5),
                                            rep("blue",5)),air=0.01,cex=1.25)
##############################


# NMDS graph --------------------------------------------------------------
Isograph <- ggplot(data = Iso_pointlocation, aes(x = MDS1, y =  MDS2, 
                                                  shape = Treatment.Status, group = Supergroup)) +
  geom_polygon(fill = NA, color = "black") +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(10, 8, 15, 17, 16, 0, 2)) +
   geom_text(aes(label = Treatment.Status), 
             vjust=-0.25, size = 2.5) +
  ggtitle(paste("Incubation Experiments: Eddy", EddySize, sep = " ")) +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  labs(y = "Axis 2") +
  theme(legend.position = "left")
Isograph


# require(gridExtra)
# grid.arrange(Isograph_1_0.2, Isograph_1_5, Isograph_2_0.2, Isograph_2_5, ncol=2)
