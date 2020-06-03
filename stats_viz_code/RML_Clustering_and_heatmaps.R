source("src/B12_Functions.R")
source("src/biostats.R")
source("src/coldiss.R")

library(cluster)
library(tidyverse)
library(vegan)
library(viridis)

# TODO: THIS NEEDS TO BE EDITED FOR THE DUPLICATES IN THE HILIC ION MODE STUFF

# User input
file.pattern <- "MSDial_B12"
my.filename <- "All.Data.Raw"

# Imports -----------------------------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_processed", pattern = regex(file.pattern, ignore_case = TRUE)))
for (i in filenames) {
  filepath <- file.path("data_processed", paste(i, ".csv", sep = ""))
  assign(my.filename, read.csv(filepath, stringsAsFactors = FALSE))
}
my.dataframe <- grep(my.filename, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)

# Reassign and prepare data
All.Data <- do.call("data.frame", get(my.dataframe)) %>%
  filter(Dataset == "B12_Incubation",
         !str_detect(Replicate.Name,
           "ProcessBlk|Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m|Std")) %>%
  separate(Replicate.Name, into = c("one", "two", "Sample.ID", "four"), remove = FALSE) %>%
  select(-c("one", "two", "four"),
         Metabolite.name, Area.Value, Replicate.Name, Sample.ID) %>%
  mutate(Experiment = ifelse(str_detect(Sample.ID, "IL1"), "Cyclonic", "Anticyclonic"),
         Area.Value = ifelse(Area.Value == 0, 1, Area.Value))
All.Data <- data.frame(lapply(All.Data, as.character), stringsAsFactors=FALSE)

# Arrange dataframe so each row is a sample and each column is a Mass Feature.
All.Data.Wide <- All.Data %>%
  select(Replicate.Name:Area.Value, Metabolite.name, Experiment) %>%
  group_by(Metabolite.name) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Metabolite.name, values_from = Area.Value)  %>%
  select(-row)

# Identify HILIC ion mode duplicates
Cmpds.with.duplicates <- All.Data %>% 
  select(Replicate.Name, Sample.ID, Experiment, Metabolite.name, Area.Value) %>%
  group_by(Metabolite.name) %>%
  summarise(n = n()) %>%
  filter(n<160) %>% ## THIS NEEDS TO BE EDITED FOR THE DUPLICATES IN THE HILIC ION MODE STUFF
  full_join(All.Data %>% select(Replicate.Name, Sample.ID, Experiment, Metabolite.name, Area.Value) %>%
            select(Metabolite.name, Experiment) %>%
            unique()) %>%
  as.data.frame() %>%
  filter(is.na(n)) %>%
  select(Metabolite.name) %>%
  unique()


spread.by.sample_RML <- All.Data.Wide %>%
  rename(treatment = Sample.ID,
         Sample.ID = Replicate.Name) %>%
  select(-(treatment:Experiment)) %>%
  unique()

result <- spread.by.sample_RML[-1]
row.names(result) <- spread.by.sample_RML$Sample.ID


spread.by.sample_RML <- ChangeClasses(spread.by.sample_RML, 1)


clustering.metabs.TZ_RML <- (decostand(spread.by.sample_RML[,-1], method = 'standardize', 1, na.rm = T))
rownames(clustering.metabs.TZ_RML) <- spread.by.sample_RML$Metabolite.name
row.names.remove <- c("Thiamine monophosphate")

clustering.metabs.TZ_RML <- clustering.metabs.TZ_RML[!(row.names(clustering.metabs.TZ_RML) %in% row.names.remove), ]

# USE CLARA for K-means clustering
## clara for clusters----------
Ave.Silh.Width_RML = c()
k = c()
for(i in 2:15) {
  Metabcl.clara <- clara(clustering.metabs.TZ_RML, k = i, metric = "euclidean",
                         samples = nrow(clustering.metabs.TZ_RML))
  Ave.Silh.Width_RML = c(Ave.Silh.Width_RML, Metabcl.clara$silinfo$avg.width)
  k = c(k, i)
}
plot(k, Ave.Silh.Width_RML, type = "p", xlab = "No. clusters")

## why is this the best number? why 6? I replaced with 17 but not sure if thats right
best.number = which(Ave.Silh.Width_RML == (max(Ave.Silh.Width_RML))) + 1
low.number = which(Ave.Silh.Width_RML[1:7]==(max(Ave.Silh.Width_RML[1:7]))) + 1
if(best.number==low.number){
  best.number = which(Ave.Silh.Width_RML[1:14]==(max(Ave.Silh.Width_RML[7:14]))) + 1
}
Metabcl.clara.high_RML <- clara(clustering.metabs.TZ_RML, k = best.number, metric = "euclidean",
                            samples = nrow(clustering.metabs.TZ_RML), correct.d = FALSE)
Metabcl.clara.low_RML <- clara(clustering.metabs.TZ_RML, k = low.number, metric = "euclidean", 
                           samples = nrow(clustering.metabs.TZ_RML), correct.d = FALSE)

## combine cluster membership with normalized data 
clust.membership_RML <- data.frame(Mass.Feature = names(Metabcl.clara.low_RML$clustering), 
                                   Cluster_low = Metabcl.clara.low_RML$clustering) %>%
  full_join(., data.frame(Mass.Feature = names(Metabcl.clara.high_RML$clustering), 
                          Cluster_best = Metabcl.clara.high_RML$clustering))
heatmap.data_RML = clustering.metabs.TZ_RML %>%
  mutate(Mass.Feature = rownames(clustering.metabs.TZ_RML)) %>%
  full_join(clust.membership_RML) %>%
  gather(Sample.ID, Value, -Mass.Feature, -Cluster_low, -Cluster_best) %>%
  left_join(All.Data %>%
              rename(treatment = Sample.ID,
                     Sample.ID = Replicate.Name) %>%
              dplyr::select(Sample.ID, Experiment, treatment) %>% 
              unique())  %>%
  unique() %>%
  arrange(Cluster_low) %>%
  select(-Sample.ID) %>%
  rename(Sample.ID = treatment)

#####
ggplot() + 
  geom_tile(data = heatmap.data_RML, aes(x = Sample.ID, y = Mass.Feature, fill = Value)) +
  facet_grid(Cluster_best~., 
             scales = "free_y",
             space = "free_y") +
  scale_fill_viridis(option = "viridis")+
  theme(#axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) +
  ggtitle("Best Cluster Number: 11")


ggplot() + 
  geom_tile(data = heatmap.data_RML, aes(x = Sample.ID, y = Mass.Feature, fill = Value)) +
  facet_grid(Cluster_low~., 
             scales = "free_y",
             space = "free_y") +
  scale_fill_viridis(option = "viridis")+
  theme(axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) +
  ggtitle("Low Cluster Number: 2")


### START HERE
heatmap.data.sum_RML <- clust.membership_RML %>%
  left_join(All.Data %>% rename(Mass.Feature = Metabolite.name)) %>%
  unique() %>%
  group_by(Cluster_low) %>%
  summarise(n = n())

ggplot(heatmap.data.sum_RML, aes(x="", y = n)) +
  geom_col(width = 1) +
  facet_wrap(~Cluster_low) +
  coord_polar("y", start = 0) +
  theme_minimal()

heatmap.data.aves_RML <- heatmap.data_RML %>%
  left_join(All.Data) %>%
              #dplyr::select(Sample.ID, Experiment, treatment, Fraction) ) %>%
  group_by(Mass.Feature, Cluster_low, Cluster_best) %>%
  mutate(meanValue = mean(Value))

ggplot() + 
  geom_tile(data = heatmap.data.aves_RML, aes(x = Mass.Feature, y = Replicate.Name, fill = meanValue)) +
  facet_grid(Cluster_best~., 
             scales = "free_y",
             space = "free_y") +
  theme(axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) 

ggplot() + 
  geom_tile(data = heatmap.data.aves_RML, aes(x = Mass.Feature, y = Replicate.Name, fill = meanValue)) +
  facet_grid(Cluster_low~., 
             scales = "free_y",
             space = "free_y") +
  theme(axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) 

# 
# ## ANOVA, MAY NOT BE NECESSARY -----------------
# all.dat <- read.csv("data_processed/IsoLagran1_0.2_notnormd.csv", stringsAsFactors = FALSE) %>%
#   select(-1)
# 
# ## wrangle -------
# mydata.median <- all.dat %>%
#   spread(Mass.Feature, Adjusted.Area) 
# # treatment.key <- all.dat %>%
# #   dplyr::select(SampID, replicate, treatment, Experiment) %>%
# #   unique() %>%
# #   filter(treatment != "0 In Situ",
# #          treatment != "1_T0") %>%
# #   mutate(NP_added = ifelse(grepl("NP", treatment),"yes","no"))
# 
# ## do the same but on all experiments separately ---------
# ## ditch metabolites in few samples
# # cmpds.in.few.smps <- mydata.median %>%
# #   gather(MassFeature, Value, -SampID, -replicate) %>%
# #   filter(!is.na(Value)) %>%
# #   group_by(MassFeature) %>%
# #   summarise(n = n()) %>%
# #   filter(n<25)
# 
# tidy.anova.test.dat.each <- mydata.median %>%
#   gather(Mass.Feature, Value, -Replicate.Name) %>%
#   #filter(!(MassFeature %in% cmpds.in.few.smps$MassFeature)) %>%
#   #full_join(.,treatment.key)%>%
#   #filter(treatment !="Blk") %>%
#   #filter(!is.na(Value)) %>%
#   separate(Replicate.Name, into = c("one", "two", "treatment", "four"), remove = FALSE) %>%
#   select(-c("one", "two", "four")) %>%
#   group_by(Mass.Feature, treatment) %>%
#   do(aov(lm(data = ., Value ~ treatment))) %>%
#   ungroup() %>% ungroup() %>%
#   filter(term==".$treatment") %>%
#   mutate(fdr.p.val = p.adjust(p.value, method = "fdr")) %>%
#   arrange(p.value)
# filter(tidy.anova.test.dat.each, fdr.p.val <= 0.05)
# 
# mass.features.that.change.with.treatments <- tidy.anova.test.dat.each %>%
#   filter(fdr.p.val <= 0.05) %>%
#   group_by(MassFeature) %>%
#   dplyr::select(MassFeature, Experiment) %>%
#   summarise(n = n()) %>%
#   left_join(., tidy.anova.test.dat.each %>%
#               mutate(fdr.p.val = ifelse(fdr.p.val <= 0.05, TRUE, FALSE)) %>%
#               dplyr::select(MassFeature, Experiment, fdr.p.val) %>%
#               spread(Experiment, fdr.p.val)) %>%
#   dplyr::rename(N.Diff = n,
#                 Diff.in.North = `Exp1 - North`,
#                 Diff.in.TZ = `Exp2 - TZ`,
#                 Diff.in.South = `Exp3 - South`)
# 
# all.dat.change <- all.dat %>%
#   filter(MassFeature %in% mass.features.that.change.with.treatments$MassFeature) %>%
#   dplyr::select(MassFeature,Alignment.ID, Average.Rt.min., Average.Mz, 
#                 Metabolite.name, MS.MS.assigned, MS.MS.spectrum,Fraction, FinalBMIS) %>%
#   unique() %>%
#   full_join(mass.features.that.change.with.treatments) 
# 
# 
# ## Do clustering just with comopunds that overlap with the transect -----
# overlap.ids <- read_csv("2020-02-06-Possible.Untargeted.Matches.With.RRExp.and.Transect.with.manual.confirmation.csv")
# confirmed.ids <- overlap.ids %>%
#   filter(Match_y.n == 'y') %>%
#   separate(key2, into = c("MassFeature.RRExp","ID1","ID2","ID3"),sep = " ")
# 
# all.dat.no.is.overlap <- all.dat.no.IS %>%
#   filter(!grepl("MF_", MassFeature)) %>%
#   full_join(., all.dat.no.IS %>%
#               filter(MassFeature %in% confirmed.ids$MassFeature.RRExp))
# 
# wide.all.data <- all.dat.no.is.overlap %>%
#   dplyr::select(MassFeature, WaterVol.Norm.Area, Sample.ID, Experiment, treatment) %>%
#   filter(!grepl("T0",treatment),
#          !grepl("Situ",treatment),
#          !grepl("_Other",Experiment),
#          !is.na(WaterVol.Norm.Area)) %>%
#   mutate(WaterVol.Norm.Area = ifelse(WaterVol.Norm.Area==0,1,WaterVol.Norm.Area)) %>%
#   spread(MassFeature, WaterVol.Norm.Area)
# 
# 
# ## separate by experiment -----------------------------
# exp.wide.data <- split.data.frame(wide.all.data, wide.all.data$Experiment)
# 
# TZ <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 2, my.name = "_subset_transect-MFs_TZ_RRExp")
# ID.TZ <- full_join(TZ,all.dat.no.is.overlap %>% dplyr::select(MassFeature, Metabolite.name) %>% unique()) %>%
#   dplyr::select(MassFeature, Metabolite.name, Cluster_3, Cluster_13) %>% 
#   unique()
# # write_csv(TZ,paste0(Sys.Date(),"_subset_transect-MFs_TZ_ExpMetabolite_Cluster_Membership.csv"))
# North <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 1, my.name = "_subset_transect-MFs_North_RRExp")
# South <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 3, my.name = "_subset_transect-MFs_South_RRExp")


source("src/biostats.R")
source("src/coldiss.R")

library(cluster)
library(tidyverse)
library(vegan)
library(viridis)

# library(readr)
# library(factoextra)
# library(broom)

# Prepare data
all.dat.no.IS_RML <- read.csv("data_processed/MSDial_B12_Transect_combined_2020-05-11.csv",
                              stringsAsFactors = FALSE) %>%
  filter(Dataset == "B12_Incubation") %>%
  separate(Replicate.Name, into = c("one", "two", "Sample.ID", "four"), remove = FALSE) %>%
  select(-c("one", "two", "four")) %>%
  select(Metabolite.name, Area.Value, Replicate.Name, Sample.ID) %>%
  filter(!str_detect(Replicate.Name,
                     "ProcessBlk|Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m|Std")) %>%
  mutate(Experiment = ifelse(str_detect(Sample.ID, "IL1"), "Cyclonic", "Anticyclonic")) %>%
  mutate(Area.Value = ifelse(Area.Value==0,1,Area.Value))

## USE CLARA for K-means clustering
## make a dataframe where each row is a sample and each column is a MF

## read data ------------------

wide.all.data_RML <- all.dat.no.IS_RML %>%
  group_by(Metabolite.name) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Metabolite.name, values_from = Area.Value)  %>%
  select(-row)

## separate by experiment -----------------------------
exp.wide.data_RML <- split.data.frame(wide.all.data_RML, wide.all.data_RML$Experiment)

## z-score by column (by each metabolite)  --------------
#clara.clusters.Angie.Heatmap <- function(exp.wide.data, number=1, my.name = "test") {

metabs.w.few_RML <- wide.all.data_RML %>% # exp.wide.data_RML[[1]] when split by experiment
  gather(Metabolite.name, Area.Value, -Sample.ID, -Experiment, -Replicate.Name) %>%
  filter(!is.na(Area.Value)) %>%
  group_by(Metabolite.name) %>%
  summarise(n = n()) %>%
  filter(n<160) %>% ## THIS NEEDS TO BE EDITED FOR THE DUPLICATES IN THE HILIC ION MODE STUFF
  #full_join(exp.wide.data[[number]] %>% ORIGINAL
  full_join(wide.all.data_RML %>%
              gather(Metabolite.name, Area.Value, -Sample.ID, -Experiment, -Replicate.Name) %>%
              filter(is.na(Area.Value)) %>%
              dplyr::select(Metabolite.name, Experiment) %>%
              unique())

#spread.by.sample <- exp.wide.data[[number]] %>% ORIGINAL
spread.by.sample_RML <- wide.all.data_RML %>%
  rename(treatment = Sample.ID,
         Sample.ID = Replicate.Name) %>%
  gather(Metabolite.name, Area.Value, -Sample.ID, -Experiment, -treatment) %>%
  filter((Metabolite.name %in% metabs.w.few_RML$Metabolite.name)) %>%
  arrange(Sample.ID) %>%
  filter(!is.na(Area.Value)) %>%
  dplyr::select(-Experiment, -treatment) %>%
  spread(Sample.ID, Area.Value) %>%
  arrange(Metabolite.name)
clustering.metabs.TZ_RML <- (decostand(spread.by.sample_RML[,-1], method = 'standardize', 1, na.rm = T))
rownames(clustering.metabs.TZ_RML) <- spread.by.sample_RML$Metabolite.name
row.names.remove <- c("Thiamine monophosphate")

clustering.metabs.TZ_RML <- clustering.metabs.TZ_RML[!(row.names(clustering.metabs.TZ_RML) %in% row.names.remove), ]


## clara for clusters----------
Ave.Silh.Width_RML = c()
k = c()
for(i in 2:15) {
  Metabcl.clara <- clara(clustering.metabs.TZ_RML, k = i, metric = "euclidean",
                         samples = nrow(clustering.metabs.TZ_RML))
  Ave.Silh.Width_RML = c(Ave.Silh.Width_RML, Metabcl.clara$silinfo$avg.width)
  k = c(k, i)
}
plot(k, Ave.Silh.Width_RML, type = "p", xlab = "No. clusters")

## why is this the best number? why 6? I replaced with 17 but not sure if thats right
best.number = which(Ave.Silh.Width_RML == (max(Ave.Silh.Width_RML))) + 1
low.number = which(Ave.Silh.Width_RML[1:7]==(max(Ave.Silh.Width_RML[1:7]))) + 1
if(best.number==low.number){
  best.number = which(Ave.Silh.Width_RML[1:14]==(max(Ave.Silh.Width_RML[7:14]))) + 1
}
Metabcl.clara.high_RML <- clara(clustering.metabs.TZ_RML, k = best.number, metric = "euclidean",
                                samples = nrow(clustering.metabs.TZ_RML), correct.d = FALSE)
Metabcl.clara.low_RML <- clara(clustering.metabs.TZ_RML, k = low.number, metric = "euclidean", 
                               samples = nrow(clustering.metabs.TZ_RML), correct.d = FALSE)

## combine cluster membership with normalized data 
clust.membership_RML <- data.frame(Mass.Feature = names(Metabcl.clara.low_RML$clustering), 
                                   Cluster_low = Metabcl.clara.low_RML$clustering) %>%
  full_join(., data.frame(Mass.Feature = names(Metabcl.clara.high_RML$clustering), 
                          Cluster_best = Metabcl.clara.high_RML$clustering))
heatmap.data_RML = clustering.metabs.TZ_RML %>%
  mutate(Mass.Feature = rownames(clustering.metabs.TZ_RML)) %>%
  full_join(clust.membership_RML) %>%
  gather(Sample.ID, Value, -Mass.Feature, -Cluster_low, -Cluster_best) %>%
  left_join(all.dat.no.IS_RML %>%
              rename(treatment = Sample.ID,
                     Sample.ID = Replicate.Name) %>%
              dplyr::select(Sample.ID, Experiment, treatment) %>% 
              unique())  %>%
  unique() %>%
  arrange(Cluster_low) %>%
  select(-Sample.ID) %>%
  rename(Sample.ID = treatment)

#####
ggplot() + 
  geom_tile(data = heatmap.data_RML, aes(x = Sample.ID, y = Mass.Feature, fill = Value)) +
  facet_grid(Cluster_best~., 
             scales = "free_y",
             space = "free_y") +
  scale_fill_viridis(option = "viridis")+
  theme(#axis.text.y  = element_blank(),
    axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
    strip.text = element_blank()) +
  ggtitle("Best Cluster Number: 11")


ggplot() + 
  geom_tile(data = heatmap.data_RML, aes(x = Sample.ID, y = Mass.Feature, fill = Value)) +
  facet_grid(Cluster_low~., 
             scales = "free_y",
             space = "free_y") +
  scale_fill_viridis(option = "viridis")+
  theme(axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) +
  ggtitle("Low Cluster Number: 2")


heatmap.data.sum <- clust.membership %>%
  left_join(all.dat.no.IS) %>%
  unique() %>%
  group_by(Cluster_3) %>%
  summarise(n = n())

ggplot(heatmap.data.sum, aes(x="", y = n)) +
  geom_col(width = 1) +
  facet_wrap(~Cluster_3) +
  coord_polar("y", start = 0) +
  theme_minimal()

heatmap.data.aves <- heatmap.data %>%
  left_join(all.dat.no.IS) %>%
  #dplyr::select(Sample.ID, Experiment, treatment, Fraction) ) %>%
  group_by(Mass.Feature, Cluster_3, Cluster_13) %>%
  mutate(meanValue = mean(Value))

ggplot() + 
  geom_tile(data = heatmap.data.aves, aes(x = Mass.Feature, y = Replicate.Name, fill = meanValue)) +
  facet_grid(Cluster_13~., 
             scales = "free_y",
             space = "free_y") +
  theme(axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) 

ggplot() + 
  geom_tile(data = heatmap.data.aves, aes(x = Mass.Feature, y = Replicate.Name, fill = meanValue)) +
  facet_grid(Cluster_3~., 
             scales = "free_y",
             space = "free_y") +
  theme(axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) 

heatmap.data


TZ <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 2, my.name = "_TZ_RRExp")
ID.TZ <- full_join(TZ, all.dat %>% dplyr::select(MassFeature, Metabolite.name) %>% unique()) %>%
  dplyr::select(MassFeature, Metabolite.name, Cluster_3, Cluster_13) %>% 
  unique()
# write_csv(TZ,paste0(Sys.Date(),"_TZ_ExpMetabolite_Cluster_Membership.csv"))
North <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 1, my.name = "_North_RRExp")
South <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 3, my.name = "_South_RRExp")


## ANOVA, MAY NOT BE NECESSARY -----------------
all.dat <- read.csv("data_processed/IsoLagran1_0.2_notnormd.csv", stringsAsFactors = FALSE) %>%
  select(-1)

## wrangle -------
mydata.median <- all.dat %>%
  spread(Mass.Feature, Adjusted.Area) 
# treatment.key <- all.dat %>%
#   dplyr::select(SampID, replicate, treatment, Experiment) %>%
#   unique() %>%
#   filter(treatment != "0 In Situ",
#          treatment != "1_T0") %>%
#   mutate(NP_added = ifelse(grepl("NP", treatment),"yes","no"))

## do the same but on all experiments separately ---------
## ditch metabolites in few samples
# cmpds.in.few.smps <- mydata.median %>%
#   gather(MassFeature, Value, -SampID, -replicate) %>%
#   filter(!is.na(Value)) %>%
#   group_by(MassFeature) %>%
#   summarise(n = n()) %>%
#   filter(n<25)

tidy.anova.test.dat.each <- mydata.median %>%
  gather(Mass.Feature, Value, -Replicate.Name) %>%
  #filter(!(MassFeature %in% cmpds.in.few.smps$MassFeature)) %>%
  #full_join(.,treatment.key)%>%
  #filter(treatment !="Blk") %>%
  #filter(!is.na(Value)) %>%
  separate(Replicate.Name, into = c("one", "two", "treatment", "four"), remove = FALSE) %>%
  select(-c("one", "two", "four")) %>%
  group_by(Mass.Feature, treatment) %>%
  do(aov(lm(data = ., Value ~ treatment))) %>%
  ungroup() %>% ungroup() %>%
  filter(term==".$treatment") %>%
  mutate(fdr.p.val = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value)
filter(tidy.anova.test.dat.each, fdr.p.val <= 0.05)

mass.features.that.change.with.treatments <- tidy.anova.test.dat.each %>%
  filter(fdr.p.val <= 0.05) %>%
  group_by(MassFeature) %>%
  dplyr::select(MassFeature, Experiment) %>%
  summarise(n = n()) %>%
  left_join(., tidy.anova.test.dat.each %>%
              mutate(fdr.p.val = ifelse(fdr.p.val <= 0.05, TRUE, FALSE)) %>%
              dplyr::select(MassFeature, Experiment, fdr.p.val) %>%
              spread(Experiment, fdr.p.val)) %>%
  dplyr::rename(N.Diff = n,
                Diff.in.North = `Exp1 - North`,
                Diff.in.TZ = `Exp2 - TZ`,
                Diff.in.South = `Exp3 - South`)

all.dat.change <- all.dat %>%
  filter(MassFeature %in% mass.features.that.change.with.treatments$MassFeature) %>%
  dplyr::select(MassFeature,Alignment.ID, Average.Rt.min., Average.Mz, 
                Metabolite.name, MS.MS.assigned, MS.MS.spectrum,Fraction, FinalBMIS) %>%
  unique() %>%
  full_join(mass.features.that.change.with.treatments) 


## Do clustering just with comopunds that overlap with the transect -----
overlap.ids <- read_csv("2020-02-06-Possible.Untargeted.Matches.With.RRExp.and.Transect.with.manual.confirmation.csv")
confirmed.ids <- overlap.ids %>%
  filter(Match_y.n == 'y') %>%
  separate(key2, into = c("MassFeature.RRExp","ID1","ID2","ID3"),sep = " ")

all.dat.no.is.overlap <- all.dat.no.IS %>%
  filter(!grepl("MF_", MassFeature)) %>%
  full_join(., all.dat.no.IS %>%
              filter(MassFeature %in% confirmed.ids$MassFeature.RRExp))

wide.all.data <- all.dat.no.is.overlap %>%
  dplyr::select(MassFeature, WaterVol.Norm.Area, Sample.ID, Experiment, treatment) %>%
  filter(!grepl("T0",treatment),
         !grepl("Situ",treatment),
         !grepl("_Other",Experiment),
         !is.na(WaterVol.Norm.Area)) %>%
  mutate(WaterVol.Norm.Area = ifelse(WaterVol.Norm.Area==0,1,WaterVol.Norm.Area)) %>%
  spread(MassFeature, WaterVol.Norm.Area)


## separate by experiment -----------------------------
exp.wide.data <- split.data.frame(wide.all.data, wide.all.data$Experiment)

TZ <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 2, my.name = "_subset_transect-MFs_TZ_RRExp")
ID.TZ <- full_join(TZ,all.dat.no.is.overlap %>% dplyr::select(MassFeature, Metabolite.name) %>% unique()) %>%
  dplyr::select(MassFeature, Metabolite.name, Cluster_3, Cluster_13) %>% 
  unique()
# write_csv(TZ,paste0(Sys.Date(),"_subset_transect-MFs_TZ_ExpMetabolite_Cluster_Membership.csv"))
North <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 1, my.name = "_subset_transect-MFs_North_RRExp")
South <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 3, my.name = "_subset_transect-MFs_South_RRExp")
