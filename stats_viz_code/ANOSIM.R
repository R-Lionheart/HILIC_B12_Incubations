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
# Indicator species are “A species whose status provides information on the 
# overall condition of the ecosystem and of other species in that ecosystem. 
# They reflect the quality and changes in environmental conditions as well 
# as aspects of community composition.”

#I often use the same groups that I use for the ANOSIM statistical test. 
# This way I can check to see which species are most responsible for the 
# differences in microbial community composition between my groups.

which.treatment <- Treatment$Treatment.Status

# The mulitpatt command results in lists of species that are associated
# with a particular group of samples. If your group has more than 2 categories,
# multipatt will also identify species that are statistically more abundant in combinations of categories.
species.indication = multipatt(df_wide_normalizedT, which.treatment,
                               func = "r.g", control = how(nperm=9999))
summary(species.indication)

# The first list contains the species found significantly more often in the “DSW” grouping. 
# The #sps 3 shows that 3 species were identified as indicators for this group. 
# The first column contains species names, the next column contains the stat value 
# (higher means the OTU is more strongly associated). The p.value column contains 
# the statistical p values for the species association (lower means stronger significance). 
# The final column shows the significance level, which is explained by the Signif. codes
# at the bottom of the output.

individual.species <- as.data.frame(species.indication$sign) %>%
  rownames_to_column() %>%
  mutate(Significant = ifelse(p.value <= 0.05, "Significant", "NotSignificant")) %>%
  rename(Mass.Feature = rowname) %>%
  pivot_longer(cols = s.B12:s.TimeZero, names_to = "SampID")

individual.graph <- individual.species %>%
  group_by(Mass.Feature) %>%
  mutate(GroupName = paste(SampID[value != 0], collapse = "_")) %>%
  mutate(testing = as.character(GroupName)) %>%
  select(-SampID, -value) %>%
  unique()


tbd.layout <- as.data.frame(species.indication$comb)
possible.goodshit <- as.data.frame(species.indication$str)

## Plot types
ind.spec.heatmap <- ggplot(data = individual.species, aes(x = Mass.Feature, y = SampID, 
                                                          fill = stat)) + 
  geom_tile(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10)) 
  #scale_y_discrete(limits = rev(levels(as.factor(heatmap.data$SampID)))) +
print(ind.spec.heatmap)


gg = ggplot(individual.species, aes(x = Mass.Feature, y = stat, fill = SampID)) + 
  geom_boxplot(colour = "black", position = position_dodge(0.5)) +
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), legend.position = "right", 
        axis.text.x = element_text(face = "bold", colour = "black", size = 12,
                                   angle = 90), 
        axis.text.y = element_text(face = "bold", size = 12, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.key=element_blank()) 
gg


ggplot(individual.species, aes(x=stat, y=p.value, size = value, color = SampID)) +
  geom_point(alpha=0.7) 

ggplot(individual.species, aes(x=p.value, y=stat)) +
  geom_point(mapping = aes(color = value), alpha = 1/20) + 
  scale_color_gradient(low="blue", high="orange") +
  coord_polar() +
  facet_wrap(~SampID)

ggplot(individual.graph, aes(Significant, stat, label = Mass.Feature)) +
  geom_point(mapping = aes(color = p.value)) +
  geom_text(position=position_jitter(width=1,height=1)) + 
  facet_wrap(~GroupName)


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

