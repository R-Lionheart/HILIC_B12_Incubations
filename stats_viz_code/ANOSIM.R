source("src/B12_Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(indicspecies)
library(tidyverse) 
library(vegan)

## Import your dataset from NMDS_figs.R
EddyInformation <- "1_Cyclonic_5um"
mydf <- IsoLagran1_5
hasChlorophyll <- "no"


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

# Standardize + distance matrix -----------------------------------------------------------------
df_wide_normalizedT <- decostand(Iso_wideT, method = "standardize", na.rm = TRUE)
df_dataframe <- as.data.frame(df_wide_normalizedT) 

ano.Treatment.Status = anosim(df_wide_normalizedT, Treatment$Treatment.Status, 
                              distance = "euclidean", permutations = 9999)
ano.Treatment.Status

ano.Control.Status = anosim(df_wide_normalizedT, Treatment$Control.Status, 
                            distance = "euclidean", permutations = 9999)
ano.Control.Status




# Indicator Species Analysis ----------------------------------------------
which.treatment <- Treatment$Treatment.Status
species.indication = multipatt(df_wide_normalizedT, which.treatment,
                               func = "r.g", control = how(nperm=9999))
summary(species.indication)
individual.species <- as.data.frame(species.indication$sign) %>%
  rownames_to_column() %>%
  mutate(Significant = ifelse(p.value <= 0.05, "Significant", "NotSignificant")) %>%
  rename(Mass.Feature = rowname) %>%
  pivot_longer(cols = s.B12:s.TimeZero, names_to = "SampID")

ggplot(data = individual.species, aes(x=SampID, y=Mass.Feature, 
                                      size=index, color=Significant)) +
  geom_point(alpha=0.5) 

myboxplot <- mydf %>% left_join(Treatment) %>% filter(Mass.Feature %in% c("Betaine", "DMSP"))
ggplot(myboxplot, aes(x = Mass.Feature, y = Adjusted.Area, fill = Treatment.Status)) + 
  geom_boxplot(colour = "black", position = position_dodge(0.5)) +
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), legend.position = "right", 
        axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(face = "bold", size = 12, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.key=element_blank()) + 
  labs(x= "", y = "Relative Abundance (%)", fill = "Time") 


# Experiment with NMDS visualizations ------------------------------------------------------------------------
Iso_wide_nmds <- vegan::metaMDS(df_wide_normalizedT, distance = "euclidean", 
                                k = 3, autotransform = FALSE, trymax = 100, wascores = FALSE)

# Visualize scree plot of potential ordination axes
dimcheckMDS(df_wide_normalizedT, distance="euclidean", k=6, autotransform=FALSE, trymax=20) 
vegan::stressplot(Iso_wide_nmds, main = paste("Stressplot, Eddy", EddyInformation, sep = " "))

# Quick vectors
myNMDS = data.frame(MDS1 = Iso_wide_nmds$points[,1], MDS2 = Iso_wide_nmds$points[,2], 
                    MDS3 = Iso_wide_nmds$points[,3]) #originally had just 2

myvec.sp <- envfit(Iso_wide_nmds$points, df_wide_normlizedT, perm=1000)
myvec.sp.df <- as.data.frame(myvec.sp$vectors$arrows*sqrt(myvec.sp$vectors$r))
myvec.sp.df$species<-rownames(myvec.sp.df)

myNMDS2 <- myNMDS %>%
  rownames_to_column()

