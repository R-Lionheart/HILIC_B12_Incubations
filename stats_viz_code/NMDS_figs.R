source("src/B12_Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(tidyverse) 
library(vegan)

### Make NMDS figures

# User data
BMIS.pattern = "BMIS_Output"
Chl.pattern = "ChlA" # Enter this to identify which files are Chlorophyll. Pattern match to your directory.
percentMissing = 0.5

# Functions
MakeNMDS <- function(mydf, hasChlorophyll) {
  
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

# BMIS import and filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

HILIC_all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  select(Mass.Feature, runDate:replicate, Adjusted.Area) %>%
  unite(Replicate.Name, runDate:replicate, sep = "_") %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m")) %>%
  filter(!Mass.Feature == "Inj_vol") %>%
  filter(!str_detect(Mass.Feature, ","))

# Adjust for T0 naming issues --------------------------------------------
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

###################################
Osmolytes <- read.csv("data_processed/Osmolytes_edited.csv", stringsAsFactors = FALSE)

IL15um_osmolytes <- IL1_5um_ChlA %>%
  filter(Mass.Feature %in% Osmolytes$Osmolytes)
  
IL25um_osmolytes <- IL2_5um_ChlA %>%
  filter(Mass.Feature %in% Osmolytes$Osmolytes)

####################################
## SWITCH SCRIPTS TO ANOSIM HERE ##
####################################

# Set data and run function ------------------------------------------------------------------------
EddyInformation <- "1_Cyclonic_5um"
df_wide_normalizedT <- makeNMDS(IsoLagran1_5, hasChlorophyll = "no")


