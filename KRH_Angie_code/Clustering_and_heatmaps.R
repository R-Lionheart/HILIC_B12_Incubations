library(tidyverse)
# library(ggplot2)
# library(tidyr)
# library(stringr)
# library(readr)
# library(dplyr)
library(vegan)
# 
library(cluster)
# library(factoextra)
# library(viridis)
# library(broom)
source("src/biostats.R")
source("src/coldiss.R")
## USE CLARA for K-means clustering
## make a dataframe where each row is a sample and each column is a MF

## read data ------------------
all.dat.no.IS <- read_csv("data_raw/2020-02-13_All_Exp_More_Filtered_MF.csv")

wide.all.data <- all.dat.no.IS %>%
  dplyr::select(MassFeature, WaterVol.Norm.Area, Sample.ID, Experiment, treatment) %>%
  filter(!grepl("T0",treatment),
         !grepl("Situ",treatment),
         !grepl("_Other",Experiment),
         !is.na(WaterVol.Norm.Area)) %>%
  mutate(WaterVol.Norm.Area = ifelse(WaterVol.Norm.Area==0,1,WaterVol.Norm.Area)) %>%
  spread(MassFeature, WaterVol.Norm.Area)


## separate by experiment -----------------------------
exp.wide.data <- split.data.frame(wide.all.data, wide.all.data$Experiment)
## z-score by column (by each metabolite)  --------------

#clara.clusters.Angie.Heatmap <- function(exp.wide.data, number=1, my.name = "test"){
  ## (Just doing TZ experiment right now) 
  
  ## get rid of metabolites only in a few samples 
  #metabs.w.few <- exp.wide.data[[number]] %>% ORIGINAL
  metabs.w.few <- exp.wide.data[[1]] %>%
    gather(MassFeature, WaterVol.Norm.Area, -Sample.ID, -Experiment, -treatment) %>%
    filter(!is.na(WaterVol.Norm.Area)) %>%
    group_by(MassFeature, Experiment) %>%
    summarise(n = n()) %>%
    filter(n<9) %>%
    #full_join(exp.wide.data[[number]] %>% ORIGINAL
    full_join(exp.wide.data[[1]] %>%
                gather(MassFeature, WaterVol.Norm.Area, -Sample.ID, -Experiment, -treatment) %>%
                filter(is.na(WaterVol.Norm.Area)) %>%
                dplyr::select(MassFeature, Experiment) %>%
                unique())
    
  #spread.by.sample <- exp.wide.data[[number]] %>% ORIGINAL
  spread.by.sample <- exp.wide.data[[1]] %>%
    gather(MassFeature, WaterVol.Norm.Area, -Sample.ID, -Experiment, -treatment) %>%
    filter(!(MassFeature %in% metabs.w.few$MassFeature)) %>%
    arrange(Sample.ID) %>%
    dplyr::select(-Experiment, -treatment) %>%
    spread(Sample.ID, WaterVol.Norm.Area) %>%
    arrange(MassFeature)
  clustering.metabs.TZ <- (decostand(spread.by.sample[,-1], method = 'standardize', 1, na.rm = T))
  rownames(clustering.metabs.TZ) <- spread.by.sample$MassFeature
  
  ## clara for clusters----------
  Ave.Silh.Width = c()
  k = c()
  for(i in 2:15){
    Metabcl.clara <- clara(clustering.metabs.TZ, k = i, metric = "euclidean" , samples = nrow(clustering.metabs.TZ))
    Ave.Silh.Width = c(Ave.Silh.Width, Metabcl.clara$silinfo$avg.width)
    k = c(k, i)
  }
  plot(k, Ave.Silh.Width, type = "p", xlab = "No. clusters")
  #jpeg(file=paste0(Sys.Date(),my.name,"_SilhouetteWidthVsNoClusters.jpeg")) ORIGINAL 
  #plot(k, Ave.Silh.Width, type = "p", xlab = "No. clusters") ORIGINAL
  #dev.off() ORIGINAL
  
  best.number =   which(Ave.Silh.Width==(max(Ave.Silh.Width))) +1
  low.number = which(Ave.Silh.Width[1:6]==(max(Ave.Silh.Width[1:6]))) +1
  if(best.number==low.number){
    best.number  = which(Ave.Silh.Width[1:14]==(max(Ave.Silh.Width[7:14]))) +1
  }
  Metabcl.clara.high <- clara(clustering.metabs.TZ, k = best.number, metric = "euclidean" , samples = nrow(clustering.metabs.TZ))
  Metabcl.clara.low <- clara(clustering.metabs.TZ, k = low.number, metric = "euclidean" , samples = nrow(clustering.metabs.TZ))
  
  ## combine cluster membership with normalized data 
  clust.membership <- data.frame(MassFeature = names(Metabcl.clara.low$clustering), Cluster_3 = Metabcl.clara.low$clustering) %>%
    full_join(.,data.frame(MassFeature = names(Metabcl.clara.high$clustering), Cluster_13 = Metabcl.clara.high$clustering) )
  heatmap.data = clustering.metabs.TZ %>%
    mutate(MassFeature = rownames(clustering.metabs.TZ)) %>%
    full_join(clust.membership) %>%
    gather(Sample.ID, Value, -MassFeature, -Cluster_3, -Cluster_13) %>%
    left_join(all.dat.no.IS %>%
                dplyr::select(Sample.ID, Experiment, treatment, Fraction) %>% 
                unique())  %>%
    unique() %>%
    arrange(Cluster_3, Fraction) 
  
  ggplot() + 
    geom_tile(data = heatmap.data, aes(x = Sample.ID, y = MassFeature, fill = Value)) +
    facet_grid(Cluster_13~., 
               scales = "free_y",
               space = "free_y") +
    #scale_fill_viridis(option = "viridis")+
    theme(axis.text.y  = element_blank(),
          axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
          strip.text = element_blank())
  #ggsave(filename = paste0(Sys.Date(),my.name,"_HeatMap_all_replicates_moreGroups.pdf"))
  
  
  ggplot() + 
    geom_tile(data = heatmap.data, aes(x = Sample.ID, y = MassFeature, fill = Value)) +
    facet_grid(Cluster_3~., 
               scales = "free_y",
               space = "free_y") +
    #scale_fill_viridis(option = "viridis")+
    theme(axis.text.y  = element_blank(),
          axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
          strip.text = element_blank())
  #ggsave(filename = paste0(Sys.Date(),my.name,"_HeatMap_all_replicates_fewerGroups.pdf"))
  
  heatmap.data.sum <- clust.membership %>%
    left_join(all.dat.no.IS %>%
                dplyr::select(MassFeature, Fraction)) %>%
    unique() %>%
    group_by(Cluster_3, Fraction) %>%
    summarise(n = n())
  
  ggplot(heatmap.data.sum, aes(x="", y = n, fill = Fraction)) +
    geom_col(width = 1) +
    facet_wrap(~Cluster_3) +
    coord_polar("y", start = 0) +
    theme_minimal()
  
  heatmap.data.aves <- heatmap.data %>%
    left_join(all.dat.no.IS %>%
                dplyr::select(Sample.ID, Experiment, treatment, Fraction) ) %>%
    group_by( Experiment, treatment, MassFeature, Cluster_3, Cluster_13) %>%
    summarise(meanValue = mean(Value))
  
  ggplot() + 
    geom_tile(data = heatmap.data.aves, aes(x = treatment, y = MassFeature, fill = meanValue)) +
    facet_grid(Cluster_13~., 
               scales = "free_y",
               space = "free_y") +
    scale_fill_viridis(option = "viridis")+
    theme(axis.text.y  = element_blank(),
          axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
          strip.text = element_blank()) 
  #ggsave(filename = paste0(Sys.Date(),my.name,"_Heatmap_averages_moreGroups.pdf"))
  
  ggplot() + 
    geom_tile(data = heatmap.data.aves, aes(x = treatment, y = MassFeature, fill = meanValue)) +
    facet_grid(Cluster_3~., 
               scales = "free_y",
               space = "free_y") +
    scale_fill_viridis(option = "viridis")+
    theme(axis.text.y  = element_blank(),
          axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
          strip.text = element_blank()) 
  #ggsave(filename = paste0(Sys.Date(),my.name,"_Heatmap_averages_fewerGroups.pdf"))
  
  heatmap.data
}

TZ <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 2, my.name = "_TZ_RRExp")
ID.TZ <- full_join(TZ,all.dat %>% dplyr::select(MassFeature, Metabolite.name) %>% unique()) %>%
  dplyr::select(MassFeature, Metabolite.name, Cluster_3, Cluster_13) %>% 
  unique()
# write_csv(TZ,paste0(Sys.Date(),"_TZ_ExpMetabolite_Cluster_Membership.csv"))
North <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 1, my.name = "_North_RRExp")
South <- clara.clusters.Angie.Heatmap(exp.wide.data, number = 3, my.name = "_South_RRExp")


## ANOVA -----------------

all.dat <- read_csv("2020-02-13_All_Exp_More_Filtered_MF.csv")

## wrangle -------
mydata.median <- all.dat %>%
  filter(!grepl("DL",MassFeature),
         !grepl("15N",MassFeature),
         !grepl("D3",MassFeature),
         # grepl("MF_",MassFeature),
         type!="Blk",
         grepl("RExp2",SampID)|grepl("RExp3",SampID)|
           grepl("RExp1",SampID),     
         !grepl("T0",SampID),
         # `Sample ID`!="RExp2T0_A"
  ) %>%
  dplyr::select(SampID, replicate,
                MassFeature, WaterVol.Norm.Area) %>%
  spread(MassFeature,WaterVol.Norm.Area ) %>%
  filter(!is.na(SampID))
treatment.key <- all.dat %>%
  dplyr::select(SampID, replicate, treatment, Experiment) %>%
  unique() %>%
  filter(treatment != "0 In Situ",
         treatment != "1_T0") %>%
  mutate(NP_added = ifelse(grepl("NP",treatment),"yes","no"))

## do the same but on all experiments separately ---------
## ditch metabolites in few samples
cmpds.in.few.smps <- mydata.median %>%
  gather(MassFeature, Value, -SampID, -replicate) %>%
  filter(!is.na(Value)) %>%
  group_by(MassFeature) %>%
  summarise(n = n()) %>%
  filter(n<25)

tidy.anova.test.dat.each <- mydata.median %>%
  gather(MassFeature, Value, -SampID, -replicate) %>%
  filter(!(MassFeature %in% cmpds.in.few.smps$MassFeature)) %>%
  full_join(.,treatment.key)%>%
  filter(treatment !="Blk") %>%
  filter(!is.na(Value)) %>%
  group_by(MassFeature, Experiment) %>%
  do(tidy(anova(lm(.$Value ~ .$treatment)))) %>%
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

write_csv(all.dat.change, paste0(Sys.Date(),
                                 "_MFs_ANOVAS_Different_treatments_Experiments_analyzed_separately_0p05.csv"))



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
