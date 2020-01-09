source("B12_Inc_Functions.R")
source("src/biostats.R")


library(pastecs) # masks dplyr + tidyr
library(stringr)
library(tidyverse) 
library(vegan)

# User data
pattern = "duplicates"
percentMissing = 0.5
badSamps <- c("Sept29QC", "TruePooWeek1", "TruePooWeek2", "TruePooWeek3", "TruePooWeek4", "DSW700m")
currentDate <- Sys.Date()

# Split the dataset into separate eddies

# Data import and first filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))

HILIC_all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  select(Mass.Feature, SampID, Adjusted.Area) %>%
  filter(!SampID %in% badSamps) %>%
  filter(!Mass.Feature == "Inj_vol") %>%
  filter(!str_detect(Mass.Feature, ",")) %>% # Drop internal standards
  group_by(Mass.Feature, SampID) %>%
  mutate(Area.Ave = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, Area.Ave) %>%
  unique()

standards <- read.csv("data_extras/Ingalls_Lab_Standards.csv", stringsAsFactors = FALSE) %>%
  rename(Mass.Feature = Compound.Name) %>%
  filter(Mass.Feature %in% HILIC_all$Mass.Feature)

HILIC_all <- HILIC_all %>%
  left_join(standards %>% select(Mass.Feature, Compound.Type) %>% unique()) %>%
  select(Mass.Feature, SampID, Area.Ave, Compound.Type) %>%
  group_by(Mass.Feature) %>%
  mutate(Missing = sum(is.nan(Area.Ave))) 

# All HILICS plotted, no filtering
plotName <- paste("figures/All_HILICS_", currentDate, ".pdf", sep = "")
pdf(plotName)
all.hilics <- ggplot(HILIC_all, aes(x = reorder(Mass.Feature, -Area.Ave), y = Area.Ave)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10))
print(all.hilics)
dev.off()

# Filter out compounds that are missing too many peaks --------------------------------------------
HILIC_filtered <- HILIC_all %>%
  group_by(Mass.Feature) %>%
  mutate(MF.Count = n()) %>%
  filter(!Missing > (percentMissing*MF.Count)) %>%
  select(-c("Missing", "MF.Count"))


# Separate dataset into groups for analysis -------------------------------------------------------
IsoLagran1 <- HILIC_filtered %>%
  filter(str_detect(SampID, "IL1|IT0"))
IsoLagran2 <- HILIC_filtered %>%
  filter(str_detect(SampID, "IL2|IT0"))

IsoLagran1_0.2 <- IsoLagran1 %>%
  filter(!str_detect(SampID, "5um"))
IsoLagran1_5 <- IsoLagran1 %>%
  filter(str_detect(SampID, "5um"))

IsoLagran2_0.2 <- IsoLagran2 %>%
  filter(!str_detect(SampID, "5um"))
IsoLagran2_5 <- IsoLagran2 %>%
  filter(str_detect(SampID, "5um"))

HILIC_wide_mid <- HILIC_data %>%
  ungroup() %>%
  tidyr::spread(SampID, Area.Ave) %>%
  as.data.frame()

HILIC_wide <- HILIC_wide_mid[,-1]
rownames(HILIC_wide) <- HILIC_wide_mid[,1]

HILIC_wide <- data.frame(HILIC_wide)
HILIC_wide[is.na(HILIC_wide)] <- NA

HILIC_noNA <- na.omit(HILIC_wide)

rm("HILIC_data")
rm("HILIC_wide_mid")

HILIC_wide <- HILIC_grouped %>%
  ungroup() %>%
  select(-Compound.Type) %>%
  tidyr::spread(SampID, Area.Ave)
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



# HEATMAPS -------------------------------------------------------------------
DMBs_IL1 <- DMBs %>%
  ungroup() %>%
  filter(str_detect(SampID, "IL1")) %>%
  mutate(SampID = as.factor(SampID))

DMBs_IL1$SampID <- factor(DMBs_IL1$SampID, 
                         levels = c("IL1Control", "IL1Control5um", "IL1DSW", "IL1DSW5um", 
                                    "IL1DMB", "IL1DMB5um", "IL1DMBnoBT", "IL1DMBnoBT5um"))

DMBs_IL1.heatmap <- ggplot(data = DMBs_IL1, aes(x = Mass.Feature, y = SampID, fill = Area.Ave)) + 
  geom_tile(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.length=unit(.25, "cm"),
        strip.text = element_text(size = 10)) +
  scale_y_discrete(limits = rev(levels(as.factor(DMBs_IL1$SampID))))
  ggtitle("HILIC Full Heatmap") 
print(DMBs_IL1.heatmap)
