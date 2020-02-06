source("B12_Functions.R")
source("src/biostats.R")

#library(data.table)
#library(pastecs) # masks dplyr + tidyr
library(plotly)
library(plyr)
library(stringr)
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
  mutate(Replicate.Name1 = plyr::mapvalues(Replicate.Name, c("171002_Smp_IT0_1", "171002_Smp_IT0_2", "171002_Smp_IT0_3"), 
                                  c("171002_Smp_IL1IT0_1", "171002_Smp_IL1IT0_2", "171002_Smp_IL1IT0_3")),
            Replicate.Name2 = plyr::mapvalues(Replicate.Name1, c("171009_Smp_IT05um_1", "171009_Smp_IT05um_2", "171009_Smp_IT05um_3"),
                                  c("171009_Smp_IL1IT05um_1", "171009_Smp_IL1IT05um_2", "171009_Smp_IL1IT05um_3")),
            Replicate.Name3 = plyr::mapvalues(Replicate.Name2, c("171016_Smp_IT0_1", "171016_Smp_IT0_2", "171016_Smp_IT0_3"),
                                  c("171016_Smp_IL2IT0_1", "171016_Smp_IL2IT0_2", "171016_Smp_IL2IT0_3")),
            Replicate.Name4 = plyr::mapvalues(Replicate.Name3, c("171023_Smp_IT05um_1", "171023_Smp_IT05um_2", "171023_Smp_IT05um_3"),
                                  c("171023_Smp_IL2IT05um_1", "171023_Smp_IL2IT05um_2", "171023_Smp_IL2IT05um_3"))) %>%
  select(Mass.Feature, Replicate.Name4, Adjusted.Area) %>%
  dplyr::rename(Replicate.Name = Replicate.Name4)
detach("package:plyr", unload=TRUE)

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
IsoLagran1_5 <- IsoLagran1 %>%
  filter(str_detect(Replicate.Name, "5um"))
write.csv(IsoLagran1_0.2, "data_processed/IsoLagran_0.2_notnormed.csv")

IsoLagran2_0.2 <- IsoLagran2 %>%
  filter(!str_detect(Replicate.Name, "5um"))
IsoLagran2_5 <- IsoLagran2 %>%
  filter(str_detect(Replicate.Name, "5um"))

rm(list = c("IsoLagran1", "IsoLagran2", "HILIC_fixed"))

# Treatment data info --------------------------------------------------------------------------
Dataset <- IsoLagran2_5

Treatment <- Dataset %>%
  ungroup() %>%
  select(Replicate.Name) %>%
  unique() %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  mutate(Control.Status = ifelse(str_detect(Supergroup, "IT0"),
                                 "Incubation", ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater", 
                                                      ifelse(str_detect(Supergroup, "Control"), "Control", "Treatments")))) %>%
  mutate(Treatment.Status = ifelse(Control.Status == "Control", "NoAddedNutrients",
                                  ifelse(Control.Status == "DeepSeaWater", "DeepNutrients",
                                         ifelse(Control.Status == "Incubation", "InSituNutrients",
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
#write.csv(Iso_wide_normalizedT, "data_processed/IsoLagran_0.2_normed.csv")

Iso_wide_nmds <- metaMDS(Iso_wide_normalizedT, distance = "euclidean", k = 2, autotransform = FALSE, trymax = 100)

stressplot(Iso_wide_nmds)


Iso_pointlocation <- Iso_wide_nmds[['points']] %>% as.data.frame() %>% cbind(Treatment)

# NMDS graph --------------------------------------------------------------

Isograph_2_5 <- ggplot(data = Iso_pointlocation, aes(x = MDS1, y =  MDS2, color = Control.Status, 
                                                  shape = Treatment.Status, group = Supergroup)) +
  geom_polygon(fill = NA, color = "black") +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(15, 16, 17, 2, 16, 16, 0)) +
  scale_color_manual(values = c("black","blue3", "chartreuse3", "red")) +   
  ggtitle("IsoLagrangian Eddy 2: 5um") +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text( size = 8),
        axis.text = element_text(size = 8))+
  labs(y = "Axis 2") +
  theme(legend.position = "left")
Isograph_2_5


require(gridExtra)
grid.arrange(Iso_graph1_0.2, Isograph_1_5, Isograph_2_0.2, Isograph_2_5, ncol=2)

############################################################################################
# B12 FOLD CHANGE --------------------------------------------
BTs <- HILIC_grouped %>%
  filter(str_detect(SampID, "noBT|WBT|DSW|Control")) %>%
  filter(!str_detect(SampID, "DMB"))

BTs_0.2_IL1 <- BTs %>%
  filter(!str_detect(SampID, "5um|IL2|700")) 

BTs_0.2_IL2 <- BTs %>%
  filter(!str_detect(SampID, "5um|IL1|700"))

BTs_5_IL1 <- BTs %>%
  filter(!str_detect(SampID, "IL2")) %>%
  filter(str_detect(SampID, "IT0|5um"))

BTs_5_IL2 <- BTs %>%
  filter(!str_detect(SampID, "IL1")) %>%
  filter(str_detect(SampID, "IT0|5um"))


## Add in standard devation for error bars
BTs_0.2_IL1 <- FindStdDev(BTs_0.2_IL1)
BTs_0.2_IL2 <- FindStdDev(BTs_0.2_IL2)
BTs_5_IL1 <- FindStdDev(BTs_5_IL1)
BTs_5_IL2 <- FindStdDev(BTs_5_IL2)

## Bar Plots
BTs_0.2_IL1_plot <- MakeBarPlot(BTs_0.2_IL1, title = "HILIC B12/noB12, 0.2um IL1")
BTs_0.2_IL2_plot <- MakeBarPlot(BTs_0.2_IL2, title = "HILIC B12/noB12, 0.2um IL2")
BTs_5_IL1_plot <- MakeBarPlot(BTs_5_IL1, title = "HILIC B12/noB12, 5um IL1")
BTs_5_IL2_plot <- MakeBarPlot(BTs_5_IL2, title = "HILIC B12/noB12, 5um IL2")

## Facet Plots
BTs_facet_0.2_IL_plot <- MakeFacetGraphs(BTs_0.2_IL1, title = "HILIC B12/noB12, 0.2um IL1",
                                         scale = "free_y")
BTs_facet_0.2_IL2_plot <- MakeFacetGraphs(BTs_0.2_IL2, title = "HILIC B12/noB12, 0.2um IL2",
                                          scale = "free_y")
BTs_facet_5_IL1_plot <- MakeFacetGraphs(BTs_5_IL1, title = "HILIC B12/noB12, 5um IL1",
                                        scale = "free_y")
BTs_facet_5_IL2_plot <- MakeFacetGraphs(BTs_5_IL2, title = "HILIC B12/noB12, 5um IL2",
                                        scale = "free_y")


# DMB FOLD CHANGE --------------------------------------------
DMBs <- HILIC_grouped %>%
  filter(!str_detect(SampID, "WBT|700")) %>%
  filter(!SampID %in% c("IL1noBT", "IT0", "IL1noBT5um", "IL1WBT5um", "IT05um", "IL2noBT", "IL2noBT5um")) 

DMBs_0.2_IL1 <- DMBs %>%
  filter(!str_detect(SampID, "5um|IL2")) 

DMBs_0.2_IL2 <- DMBs %>%
  filter(!str_detect(SampID, "5um|IL1"))

DMBs_5_IL1 <- DMBs %>%
  filter(!str_detect(SampID, "IL2")) %>%
  filter(str_detect(SampID, "IT0|5um"))

DMBs_5_IL2 <- DMBs %>%
  filter(!str_detect(SampID, "IL1")) %>%
  filter(str_detect(SampID, "IT0|5um"))

## Add in standard devation for error bars
DMBs_0.2_IL1 <- FindStdDev(DMBs_0.2_IL1)
DMBs_0.2_IL2 <- FindStdDev(DMBs_0.2_IL2)
DMBs_5_IL1 <- FindStdDev(DMBs_5_IL1)
DMBs_5_IL2 <- FindStdDev(DMBs_5_IL2)

## Plots
DMBs_0.2_IL1_plot <- MakeBarPlot(DMBs_0.2_IL1, title = "HILIC DMB/noDMB, 0.2um IL1")
DMBs_0.2_IL2_plot <- MakeBarPlot(DMBs_0.2_IL2, title = "HILIC DMB/noDMB, 0.2um IL2")
DMBs_5_IL1_plot <- MakeBarPlot(DMBs_5_IL1, title = "HILIC DMB/noDMB, 5um IL1")
DMBs_5_IL2_plot <- MakeBarPlot(DMBs_5_IL2, title = "HILIC DMB/noDMB, 5um IL2")

DMBs_facet_0.2_IL_plot <- MakeFacetGraphs(DMBs_0.2_IL1, title = "HILIC DMB/noDMB, 0.2um IL1",
                                          scale = "free_y")
DMBs_facet_0.2_IL2_plot <- MakeFacetGraphs(DMBs_0.2_IL2, title = "HILIC DMB/noDMB, 0.2um IL2",
                                           scale = "free_y")
DMBs_facet_5_IL1_plot <- MakeFacetGraphs(DMBs_5_IL1, title = "HILIC DMB/noDMB, 5um IL1",
                                         scale = "free_y")
DMBs_facet_5_IL2_plot <- MakeFacetGraphs(DMBs_5_IL2, title = "HILIC DMB/noDMB, 5um IL2",
                                         scale = "free_y")


## Aminos plot
aminos <- HILIC_grouped %>%
  filter(Compound.Type == "Amino Acid") 
aminos <- FindStdDev(aminos)

aminos.5um <- aminos %>%
  filter(!SampID == "IL1Control") %>%
  filter(!SampID == "IL1DSW") %>%
  filter(!SampID == "IL2Control") %>%
  filter(!SampID == "IL2DSW") %>%
  filter(str_detect(SampID, "5um|Control|DSW")) %>%
  arrange(SampID)

aminos.plot.5um <- ggplot(aminos.5um, aes(x = Mass.Feature, y = Area.Ave, fill = SampID)) +
  geom_bar(stat = "identity") + 
  facet_wrap( ~SampID, scales = "fixed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  geom_errorbar(aes(ymin=Area.Ave - Std.dev, ymax=Area.Ave + Std.dev), width=.2,
                position=position_dodge(.9)) +
  ggtitle("HILIC Aminos data")
print(aminos.plot.5um)



