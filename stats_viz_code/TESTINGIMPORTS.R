source("src/B12_Functions.R")
library(tidyverse)
library(vegan)
options(scipen = 999)

## Import and filtering for stats analysis ##

ChlMakeWide <- function(df) {
  df.wide <- df %>%
    ungroup() %>%
    tidyr::spread(Replicate.Name, Normalized.by.Chla) %>%
    as.data.frame()
  
  df.rownames <- df.wide[,-1]
  rownames(df.rownames) <- df.wide[,1]
  
  df.rownames[is.na(df.rownames)] <- NA
  
  return(df.rownames)  
}
UploadChlorophyll <- function(whichEddy) {
  df <- read.csv("data_raw/Dyhrman_MS_Chla.csv", stringsAsFactors = FALSE) %>%
    select(1:6) %>%
    rename(Filter.size = Filter..µm.,
           Chla = Chl.a..µg.L.) %>%
    mutate(Date = as.character(Date),
           Eddy = ifelse(str_detect(Date, "13|9"), "IL2", "IL1")) %>%
    filter(str_detect(Eddy, whichEddy)) %>%
    unite("Replicate.Name", SAMPLE, Eddy, sep = "_", remove = TRUE) %>%
    filter(Filter.size == 5.0) %>%
    mutate(Chla = as.numeric(Chla))
  
  return(df)
}
UploadFiles <- function(myfilepath) {
  # Function for uploading non-standardized files.
  #
  # Returns: dataframe with relevant columns selected.
  uploaded.df <- read.csv(myfilepath, stringsAsFactors = FALSE) %>%
    select(Mass.Feature:Adjusted.Area)
  
  return(uploaded.df)
}


# User data
BMIS.pattern = "BMIS_Output"
#Chl.pattern = "ChlA" # Enter this to identify which files are Chlorophyll. Pattern match to your directory.
percentMissing = 0.5 # Use as percent cutoff for dropping compounds missing this percent of peaks.

# Set filtering conditions that correspond to the treatments you are comparing.
# Condition1 <- "IL1DSW5um" # Other options: IL1DMBnoBT, IL2WBT, IL1noBt, etc.
# Condition2 <- "IL1Control5um"
# FilterSize <- "5um"
ChlaEddy <- "IL1"
# SigValue <- "pvalue" # alternative is "qvalue", when using fdr-corrected values.
# file.pattern <- "Chla Normalized Cyclonic_5um" # will be used as a search ID and title for graphs 
# SigNumber <- 0.1 # Pvalue cutoff
# BMISd <- BMISd_1_5_notnormd # Assign correct dataframe for analysis.

# BMIS import and filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

all.HILIC <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  select(Mass.Feature, runDate:replicate, Adjusted.Area) %>%
  unite(Replicate.Name, runDate:replicate, sep = "_") %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m"),
         Mass.Feature != "Inj_vol",
         !str_detect(Mass.Feature, ",")) # This line drops Internal Standards. Temporary fix!

# Adjust for T0 naming issues
all.HILIC <- all.HILIC %>%
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

# Filter out compounds that are missing too many peaks --------------------------------------------
HILIC_filtered <- all.HILIC %>%
  group_by(Mass.Feature) %>%
  mutate(Missing = sum(is.na(Adjusted.Area))) %>%
  mutate(MF.Count = n()) %>%
  filter(!Missing > (percentMissing*MF.Count)) %>%
  select(-c("Missing", "MF.Count"))

# Separate groups into analysis and save to data_processed/ subdirectory.
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



## Data normalized to Chlorophyll-a
## Upload and recode chlorophyll data


# # Upload all required files
# BMISd_1_0.2_notnormd <- UploadFiles("data_processed/IsoLagran1_0.2_notnormd.csv")
# BMISd_1_5_notnormd <- UploadFiles("data_processed/IsoLagran1_5_notnormd.csv")
# BMISd_2_0.2_notnormd <- UploadFiles("data_processed/IsoLagran2_0.2_notnormd.csv")
# BMISd_2_5_notnormd <- UploadFiles("data_processed/IsoLagran2_5_notnormd.csv")

Chlorophyll.data <- UploadChlorophyll(whichEddy = ChlaEddy)
BMISd <- IsoLagran2_5

# Fix Replicate Name labels 
Chlorophyll_fixed <- Chlorophyll.data %>%
  mutate(Replicate.Name = recode(Replicate.Name,
                                 "C-1_IL1" ="171009_Smp_IL1Control5um_1",
                                 "C-1_IL2" = "171009_Smp_IL2Control5um_1",
                                 "C-2_IL1" = "171009_Smp_IL1Control5um_2",
                                 "C-2_IL2" = "171009_Smp_IL2Control5um_2",
                                 "C-3_IL1" = "171009_Smp_IL1Control5um_3",
                                 "C-3_IL2" = "171009_Smp_IL2Control5um_3",
                                 "ØBT-1_IL1" = "171009_Smp_IL1noBT5um_1",
                                 "ØBT-1_IL2" = "171009_Smp_IL2noBT5um_1",
                                 "ØBT-2_IL1" = "171009_Smp_IL1noBT5um_2",
                                 "ØBT-2_IL2" = "171009_Smp_IL2noBT5um_2",
                                 "ØBT-3_IL1" = "171009_Smp_IL1noBT5um_3",
                                 "ØBT-3_IL2" = "171009_Smp_IL2noBT5um_3",
                                 "DMB-1_IL1" = "171009_Smp_IL1DMB5um_1",
                                 "DMB-1_IL2" = "171009_Smp_IL2DMB5um_1",
                                 "DMB-2_IL1" = "171009_Smp_IL1DMB5um_2",
                                 "DMB-2_IL2" = "171009_Smp_IL2DMB5um_2",
                                 "DMB-3_IL1" = "171009_Smp_IL1DMB5um_3",
                                 "DMB-3_IL2" = "171009_Smp_IL2DMB5um_3",
                                 "DMBØBT-1_IL1" = "171009_Smp_IL1DMBnoBT5um_1",
                                 "DMBØBT-1_IL2" = "171009_Smp_IL2DMBnoBT5um_1",
                                 "DMBØBT-2_IL1" = "171009_Smp_IL1DMBnoBT5um_2",
                                 "DMBØBT-2_IL2" = "171009_Smp_IL2DMBnoBT5um_2",
                                 "DMBØBT-3_IL1" = "171009_Smp_IL1DMBnoBT5um_3",
                                 "DMBØBT-3_IL2" = "171009_Smp_IL2DMBnoBT5um_3",
                                 "DSW-1_IL1" = "171009_Smp_IL1DSW5um_1",
                                 "DSW-1_IL2" = "171009_Smp_IL2DSW5um_1",
                                 "DSW-2_IL1" = "171009_Smp_IL1DSW5um_2",
                                 "DSW-2_IL2" = "171009_Smp_IL2DSW5um_2",
                                 "DSW-3_IL1" = "171009_Smp_IL1DSW5um_3",
                                 "DSW-3_IL2" = "171009_Smp_IL2DSW5um_3",
                                 "IS-1_IL1" = "171009_Smp_IL1IT05um_1",
                                 "IS-1_IL2" = "171009_Smp_IL2IT05um_1",
                                 "IS-2_IL1" = "171009_Smp_IL1IT05um_2",
                                 "IS-2_IL2" = "171009_Smp_IL2IT05um_2",
                                 "IS-3_IL1" = "171009_Smp_IL1IT05um_3",
                                 "IS-3_IL2" = "171009_Smp_IL2IT05um_3",
                                 "WBT-1_IL1" = "171009_Smp_IL1WBT5um_1",
                                 "WBT-1_IL2" = "171009_Smp_IL2WBT5um_1",
                                 "WBT-2_IL1" = "171009_Smp_IL1WBT5um_2",
                                 "WBT-2_IL2" = "171009_Smp_IL2WBT5um_2",
                                 "WBT-3_IL1" = "171009_Smp_IL1WBT5um_3",
                                 "WBT-3_IL2" = "171009_Smp_IL2WBT5um_3")) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "replicate")) %>%
  select(SampID, replicate, Chla) %>%
  unite(SampID, replicate, col = "Replicate.Name")

# Fix BMISd names to match ------------------------------------------------------------
BMISd_fixed <- BMISd %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "replicate")) %>%
  select(-c("one", "two")) %>%
  unite(SampID, replicate, col = "Replicate.Name")

# Normalize to chlorophyll ------------------------------------------------------------
complete.set <- BMISd_fixed %>%
  left_join(Chlorophyll_fixed %>% select(Replicate.Name, Chla)) %>%
  mutate(Chla = as.numeric(Chla)) %>%
  mutate(Normalized.by.Chla = Adjusted.Area/Chla) %>%
  filter(!Chla < 0) %>% # Remove negative Chl-A values
  select(Mass.Feature, Replicate.Name, Normalized.by.Chla) %>%
  na.omit()

## Standardize
complete.wide <- ChlMakeWide(complete.set)
complete.wide[is.na(complete.wide)] <- 1000
complete.wideT <- t(complete.wide)

complete.wide.normalizedT <- decostand(complete.wideT, method = "standardize", na.rm = TRUE) 

## Save Chlorophyll files to data_processed/
write.csv(complete.wide.normalizedT, paste("data_processed/", ChlaEddy, "_0.5um_ChlA_normd_std.csv", sep = ""))
write.csv(complete.set, paste("data_processed/", ChlaEddy, "_0.5um_ChlA_normd_nostd.csv", sep = ""))
write.csv(Chlorophyll_fixed, "data_processed/Chlorophyll_fixed.csv")





IL1_5um_ChlA <- IL1_5um_ChlA_normd_nostd %>%
  select(-X) %>%
  rename(Adjusted.Area = Normalized.by.Chla)
IL2_5um_ChlA <- IL2_5um_ChlA_normd_nostd %>%
  select(-X) %>%
  rename(Adjusted.Area = Normalized.by.Chla)
rm(list = c("IsoLagran1", "IsoLagran2", "IL1_5um_ChlA_normd_nostd", "IL2_5um_ChlA_normd_nostd"))

# Import Chlorophyll ------------------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = regex(Chl.pattern, ignore_case = TRUE)))

for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}