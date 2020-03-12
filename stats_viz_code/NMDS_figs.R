source("B12_Functions.R")
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

# Data import and first filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))

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

data4cmap <- HILIC_fixed %>%
  separate(Replicate.Name, into = c("Date", "runtype", "SampID", "replicate")) %>%
  select(Mass.Feature, SampID, Adjusted.Area) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Average.Area.Normed = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, Average.Area.Normed) %>%
  unique()

# All HILICS plotted, no filtering
all.hilics <- ggplot(HILIC_fixed, aes(x = reorder(Mass.Feature, -Adjusted.Area), y = Adjusted.Area)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#print(all.hilics)

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
Dataset <- IsoLagran1_0.2
Dataset <- IsoLagran1_5
Dataset <- IsoLagran2_0.2
Dataset <- IsoLagran2_5

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

# KRH transformations --------------------------------------------------------------
Iso_wide_normalizedT <- decostand(Iso_wideT, method = "standardize", na.rm = TRUE) 
write.csv(Iso_wide_normalizedT, "data_processed/IsoLagran2_0.2_normd.csv")

Iso_wide_nmds <- vegan::metaMDS(Iso_wide_normalizedT, distance = "euclidean", k = 2, autotransform = FALSE, trymax = 100)

stressplot(Iso_wide_nmds)


Iso_pointlocation <- Iso_wide_nmds[['points']] %>% as.data.frame() %>% cbind(Treatment)

Iso_pointlocation$Treatment.Status <- factor(Iso_pointlocation$Treatment.Status,
                              levels = c("TimeZero", "Control", "DMB", "B12", "DeepSeaWater",
                                         "DMBnoB12", "noB12"))

# NMDS graph --------------------------------------------------------------

Isograph_2_5 <- ggplot(data = Iso_pointlocation, aes(x = MDS1, y =  MDS2, 
                                                  shape = Treatment.Status, group = Supergroup)) +
  geom_polygon(fill = NA, color = "black") +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(10, 8, 15, 17, 16, 0, 2)) +
  # geom_text(aes(label = Treatment.Status), 
  #           vjust=-0.25, size = 2.5) +
  ggtitle("IsoLagrangian Eddy 2: 5um") +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))+
  labs(y = "Axis 2") +
  theme(legend.position = "left")
Isograph_2_5


require(gridExtra)
grid.arrange(Isograph_1_0.2, Isograph_1_5, Isograph_2_0.2, Isograph_2_5, ncol=2)
