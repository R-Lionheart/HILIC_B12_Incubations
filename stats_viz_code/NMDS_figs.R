source("src/B12_Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(tidyverse) 
library(vegan)

# User data
BMIS.pattern = "BMIS_Output"
Chl.pattern = "ChlA"
percentMissing = 0.5

# Functions
makeNMDS <- function(mydf, hasChlorophyll) {
  
  Treatment <- mydf %>%
    ungroup() %>%
    select(Replicate.Name) %>%
    unique() 
  
  if (str_detect(hasChlorophyll, regex("no", ignore_case = TRUE))) {
    Treatment <- Treatment %>%
      separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) 
    
  } else if (str_detect(hasChlorophyll, regex("yes", ignore_case = TRUE))) {
    Treatment <- Treatment %>%
      separate(Replicate.Name, into = c("Supergroup", "replicate"), remove = FALSE)
  }
  
  Treatment <- Treatment %>%
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
  Iso_wide <- makeWide(mydf)
  Iso_wide[is.na(Iso_wide)] <- 1000
  Iso_wideT <- t(Iso_wide)
  
  # Standardize + distance matrix --------------------------------------------------------------
  df_wide_normalizedT <- decostand(Iso_wideT, method = "standardize", na.rm = TRUE)
  df_dataframe <- as.data.frame(df_wide_normalizedT)
  write.csv(df_wide_normalizedT, paste("data_processed/IsoLagran", EddyInformation, "_normd.csv", sep = ""))
  
  Iso_wide_nmds <- vegan::metaMDS(df_wide_normalizedT, distance = "euclidean", 
                                  k = 3, autotransform = FALSE, trymax = 100, wascores = FALSE)
  # Assign treatment to points ----------------------------------------------
  Iso_pointlocation <- Iso_wide_nmds$points %>% as.data.frame() %>% cbind(Treatment)
  Iso_pointlocation$Treatment.Status <- factor(Iso_pointlocation$Treatment.Status,
                                               levels = c("TimeZero", "Control", "DMB", "B12", "DeepSeaWater",
                                                          "DMBnoB12", "noB12"))
  # Plot NMDS graph ----------------------------------------------
  Isograph <- ggplot() + 
    geom_polygon(data=Iso_pointlocation, aes(x=MDS1, y=MDS2, fill=Treatment.Status, group=Treatment.Status), alpha=0.30) +
    geom_point(data=Iso_pointlocation, aes(x=MDS1,y=MDS2,colour=Treatment.Status),size=4) + 
    geom_text(data=Iso_pointlocation,aes(x=MDS1,y=MDS2,label=Treatment.Status), size=4) +  # add the species labels
    xlim(-20, 10) +
    ggtitle(paste("Incubation Experiments: Eddy", EddyInformation, sep = " ")) 
  print(Isograph)
  ggsave(path = "figures", paste("Incubation Experiments: Eddy", EddyInformation, ".png", sep = ""))
  #ggsave("figures/", paste("Incubation Experiments: Eddy", EddyInformation, ".png", sep = ""))
  
  return(Iso_wide_nmds)
}
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
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

HILIC_all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  select(Mass.Feature, runDate:replicate, Adjusted.Area) %>%
  unite(Replicate.Name, runDate:replicate, sep = "_") %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m")) %>%
  filter(!Mass.Feature == "Inj_vol") %>%
  filter(!str_detect(Mass.Feature, ","))

# Adjust for ILnT0 naming issues --------------------------------------------
HILIC_all <- HILIC_all %>%
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

# Chlorophyll import --------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = regex(Chl.pattern, ignore_case = TRUE)))

for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

# Filter out compounds that are missing too many peaks --------------------------------------------
HILIC_filtered <- HILIC_all %>%
  group_by(Mass.Feature) %>%
  mutate(Missing = sum(is.na(Adjusted.Area))) %>%
  mutate(MF.Count = n()) %>%
  filter(!Missing > (percentMissing*MF.Count)) %>%
  select(-c("Missing", "MF.Count"))

# Separate groups into analysis
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

IL1_5um_ChlA <- IL1_5um_ChlA_normd_nostd %>%
  select(-X) %>%
  rename(Adjusted.Area = Normalized.by.Chla)
IL2_5um_ChlA <- IL2_5um_ChlA_normd_nostd %>%
  select(-X) %>%
  rename(Adjusted.Area = Normalized.by.Chla)

rm(list = c("IsoLagran1", "IsoLagran2", "IL1_5um_ChlA_normd_nostd", "IL2_5um_ChlA_normd_nostd"))

# Set data and run function ------------------------------------------------------------------------
EddyInformation <- "2_Anticyclonic_5um_noChl"
df_wide_normalizedT <- makeNMDS(IsoLagran2_5, hasChlorophyll = "no")

# Experiment with NMDS visualizations ------------------------------------------------------------------------
Iso_wide_nmds <- vegan::metaMDS(df_wide_normalizedT, distance = "euclidean", 
                                k = 3, autotransform = FALSE, trymax = 100, wascores = FALSE)


# Visualize scree plot of potential ordination axes
dimcheckMDS(df_wide_normalizedT, distance="euclidean", k=6, autotransform=FALSE, trymax=20) 
vegan::stressplot(Iso_wide_nmds, main = paste("Stressplot, Eddy", EddyInformation, sep = " "))
# Check stressplots, scree diagrams
#Iso_wide_nmds$stress # Add a flag if this is high?
#nmds.monte(df_wide_normlizedT, distance="euclidean", k=3, autotransform=FALSE, trymax=20)
# Samples in ordinate space -----------------------------------------------
# Plot 2 dimensional NMDS configuration.
# plot(Iso_wide_nmds$points, type="n") # plotting the scores(iso_wide_nmds)
# text(Iso_wide_nmds,labels=row.names(Iso_wideT), cex = 1)
# title(paste("Incubation Experiments: Eddy", EddyInformation, sep = " "))

# See how particular compound changes with location
# iso_df <- as.data.frame(df_wide_normlizedT)
# plot(scores(Iso_wide_nmds), type = "p")
# points(Iso_wide_nmds, cex = iso_df$Ectoine, col = "red")
# title(paste("Incubation Experiments: Eddy", EddyInformation, sep = " "))


# Quick vectors
myNMDS = data.frame(MDS1 = Iso_wide_nmds$points[,1], MDS2 = Iso_wide_nmds$points[,2], 
                    MDS3 = Iso_wide_nmds$points[,3]) #originally had just 2

myvec.sp <- envfit(Iso_wide_nmds$points, df_wide_normlizedT, perm=1000)
myvec.sp.df <- as.data.frame(myvec.sp$vectors$arrows*sqrt(myvec.sp$vectors$r))
myvec.sp.df$species<-rownames(myvec.sp.df)

myNMDS2 <- myNMDS %>%
  rownames_to_column()


# Vectors and ordiplots - messy due to scale differences
# ggplot(data = myNMDS, aes(MDS1, MDS2)) + 
#   #geom_point() + # comment out for neatness
#   geom_segment(data=myvec.sp.df, aes(x=0, xend=MDS1, y=0, yend=MDS2),
#                #arrow = arrow(),
#                colour="grey") +
#   #geom_text(data=myNMDS2, aes(x=MDS1, y=MDS2, label=rowname), size=3) +
#   geom_text(data=myvec.sp.df, aes(x=MDS1, y=MDS2, label=species), size=3) +
#   ggtitle(paste("Incubation Experiments: Eddy", EddyInformation, sep = " "))


# NMDS combined with clustering
# mysol.d <- vegdist(df_wide_normlizedT, "euclidean")
# mysitecl.ward <- hclust(mysol.d,method='ward.D')
# mysitecl.class <- cutree(mysitecl.ward, k=7) # customize k groups here
# groups <- levels(factor(mysitecl.class))
# mysite.sc <- scores(myNMDS)
# 
# my.p <- ordiplot(mysite.sc, type="n", 
#                  main=paste("Incubation Experiments: Eddy", EddyInformation, sep = " "))
# for (i in 1:length(groups))
# {
#   points(mysite.sc[mysitecl.class==i,], pch=(14+i), cex=2, col=i+1)
# }
# text(mysite.sc, row.names(Iso_wideT), pos=4, cex=0.7)
# ordicluster(my.p, mysitecl.ward, col="dark grey")
# 
# legend("bottomleft", paste("Group", c(1:length(groups))),
#        pch=14+c(1:length(groups)), col=1+c(1:length(groups)), pt.cex=2)
# 
# Ordiplot with hulls CURRENTLY NOT WORKING
# ordiplot(Iso_wide_nmds, type="n")
# ordihull(Iso_wide_nmds, groups=Treatment$Treatment.Status, draw="polygon",col="grey90",label=T)
# orditorp(Iso_wide_nmds, display="sites",
#          air=0.01, cex=0.5)
# 
# # Plot convex hulls with colors based on treatment NOT WORKING FOR SAME REASON as above
# colors=c(rep("red",5), rep("blue",5))
# ordiplot(Iso_wide_nmds, type="n")
# for(i in unique(Treatment$Treatment.Status)) {
#   ordihull(Iso_wide_nmds$point[grep(i,Treatment$Treatment.Status),], draw="polygon",
#            groups=Treatment$Treatment.Status[Treatment$Treatment.Status==i],
#            col=colors[grep(i, Treatment$Treatment.Status)], label=T)} 
# #orditorp(Iso_wide_nmds, display="species", col="red", air=0.01)
# orditorp(Iso_wide_nmds, display="sites", air=0.01, cex=1.25)



# Experiments with species scores, not very useful
# envfit(Iso_wide_nmds, df_wide_normlizedT) 
# test <- envfit(df_wide_normlizedT ~ Iso_pointlocation$Supergroup, data = iso_df, perm=1000) #???
# 
# scores(Iso_wide_nmds) # this is the same as the $points call, species scores do not appear to exist
## Stacked graphs
top.15 <- all.hilics.data %>%
  arrange(desc(Total.Average)) %>%
  head(15)
stacked.hilics.data <- HILIC_fixed %>%
  mutate(Full.Total = sum(Adjusted.Area, na.rm = TRUE)) %>%
  filter(Mass.Feature %in% top.15$Mass.Feature) %>%
  #filter(str_detect(Replicate.Name, Treatments)) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four")) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(My.Average = mean(Adjusted.Area, na.rm = TRUE)) %>%
  mutate(Percent.Total = (My.Average / Full.Total)) %>%
  select(Mass.Feature, SampID, My.Average, Percent.Total) %>%
  unique()

# 1_0.2
stacked.hilics.data$SampID <- factor(stacked.hilics.data$SampID, levels = 
                                       c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT", 
                                         "IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT",
                                         "IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um",
                                         "IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))

ggplot(stacked.hilics.data, aes(fill=Mass.Feature, y=My.Average, x=SampID)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, size = 15, vjust = .55),
        legend.text=element_text(size=15)) +
  ggtitle("Top 15 Most Abundant Compounds")

ggplot(stacked.hilics.data, aes(fill=Mass.Feature, y=Percent.Total, x=SampID)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, size = 15, vjust = .55),
        legend.text = element_text(size = 15)) +
  ggtitle("Top 15 Most Abundant Compounds: Percentages")

# All HILICS plotted, no filtering
all.hilics.data <- HILIC_all %>%
  group_by(Mass.Feature) %>%
  mutate(Total.Average = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, Total.Average) %>%
  unique()

all.hilics <- ggplot(all.hilics.data, aes(x = reorder(Mass.Feature, -Total.Average), 
                                          y = Total.Average)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#print(all.hilics)
