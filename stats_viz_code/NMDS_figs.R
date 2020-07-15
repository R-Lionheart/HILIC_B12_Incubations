source("src/B12_Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(tidyverse) 
library(vegan)

### Make NMDS figures

# User data
BMIS.pattern = "IsoLagran"
Chl.pattern = "ChlA_normd_notstd" # Enter this to identify which files are Chlorophyll. Pattern match to your directory.
# percentMissing = 0.5

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
  Iso_wide <- MakeWide(mydf)
  Iso_wide[is.na(Iso_wide)] <- 1000
  Iso_wideT <- t(Iso_wide)
  
  # Standardize + distance matrix --------------------------------------------------------------
  df_wide_normalizedT <- decostand(Iso_wideT, method = "standardize", na.rm = TRUE)
  df_dataframe <- as.data.frame(df_wide_normalizedT)
  write.csv(df_wide_normalizedT, paste("data_processed/IsoLagran", EddyInformation, "_std.csv", sep = ""))
  
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
    geom_point(data=Iso_pointlocation, aes(x=MDS1,y=MDS2, colour=Treatment.Status),size=4) + 
    geom_text(data=Iso_pointlocation,aes(x=MDS1,y=MDS2,label=Treatment.Status), size=4) +  # add the species labels
    xlim(-20, 10) +
    ggtitle(paste("Incubation Experiments: Eddy", EddyInformation, sep = " ")) 
  print(Isograph)
  ggsave(path = "figures", paste("Incubation Experiments: Eddy", EddyInformation, ".png", sep = ""))

  return(Iso_wide_nmds)
}

# BMIS import and filtering of unnecessary SampIDs --------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filenames, ".csv", sep = ""))
for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE) %>%
           select(-X))
}

# Chlorophyll import --------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = regex(Chl.pattern, ignore_case = TRUE)))

for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE) %>%
           select(-X) %>%
           rename(Adjusted.Area = Normalized.by.Chla))
}


# Set data and run function ------------------------------------------------------------------------
EddyInformation <- "2_Anticyclonic_5_ChlA"
df_wide_normalizedT <- MakeNMDS(IL2_5um_ChlA_normd_notstd, hasChlorophyll = "yes")


###################################
# Osmolytes <- read.csv("data_processed/Osmolytes_edited.csv", stringsAsFactors = FALSE)
# 
# IL15um_osmolytes <- IL1_0.5um_ChlA_normd_notstd %>%
#   filter(Mass.Feature %in% Osmolytes$Osmolytes)
#   
# IL25um_osmolytes <- IL2_0.5um_ChlA_normd_notstd %>%
#   filter(Mass.Feature %in% Osmolytes$Osmolytes)





