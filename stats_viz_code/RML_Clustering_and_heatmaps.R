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

ID.levels <- c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT",    
               "IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT",
               "IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um", "IL1DMBnoBT5um", "IL1noBT5um",
               "IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um", "IL2DMBnoBT5um", "IL2noBT5um")



filenames <- RemoveCsv(list.files(path = "data_processed", pattern = regex(file.pattern, ignore_case = TRUE)))

for (i in filenames) {
  filepath <- file.path("data_processed", paste(i, ".csv", sep = ""))
  assign(my.filename, read.csv(filepath, stringsAsFactors = FALSE))
}

my.dataframe <- grep(my.filename, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
All.Data <- do.call("data.frame", get(my.dataframe)) 


# Tidy and filter data ----------------------------------------------------
All.Data.Filtered <- All.Data %>%
  filter(Dataset == "B12_Incubation") %>%
  separate(Replicate.Name, into = c("one", "two", "Sample.ID", "four"), remove = FALSE) %>%
  select(-c("one", "two", "four")) %>%
  select(Metabolite.name, Area.Value, Replicate.Name, Sample.ID) %>%
  filter(!str_detect(Replicate.Name,
                     "ProcessBlk|Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m|Std")) %>%
  mutate(Experiment = ifelse(str_detect(Sample.ID, "IL1"), "Cyclonic", "Anticyclonic")) %>%
  mutate(Area.Value = ifelse(Area.Value == 0, 1, Area.Value))


# Pivot dataframe to wide -------------------------------------------------
# Make a dataframe where each row is a sample and each column is a MF

Wide.All.Data.Filtered <- All.Data.Filtered %>%
  group_by(Metabolite.name) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Metabolite.name, values_from = Area.Value)  %>%
  select(-row)


# Identify duplicates and filter by replicate counts ----------------------
Metabolite.Counts <- Wide.All.Data.Filtered %>%
  gather(Metabolite.name, Area.Value, -Sample.ID, -Experiment, -Replicate.Name) %>%
  filter(!is.na(Area.Value)) %>%
  group_by(Metabolite.name) %>%
  summarise(n = n()) %>%
  filter(n<160) %>% ## THIS NEEDS TO BE EDITED FOR THE DUPLICATES IN THE HILIC ION MODE STUFF
  full_join(Wide.All.Data.Filtered %>%
              gather(Metabolite.name, Area.Value, -Sample.ID, -Experiment, -Replicate.Name) %>%
              filter(is.na(Area.Value)) %>%
              select(Metabolite.name, Experiment) %>%
              unique())


# Pivot wider and spread by sample ----------------------------------------

Spread.by.Sample <- Wide.All.Data.Filtered %>%
  rename(treatment = Sample.ID,
         Sample.ID = Replicate.Name) %>%
  gather(Metabolite.name, Area.Value, -Sample.ID, -Experiment, -treatment) %>%
  filter((Metabolite.name %in% Metabolite.Counts$Metabolite.name)) %>%
  arrange(Sample.ID) %>%
  filter(!is.na(Area.Value)) %>%
  select(-Experiment, -treatment) %>%
  spread(Sample.ID, Area.Value) %>%
  arrange(Metabolite.name)


# Standardize wide dataset ------------------------------------------------
Standardized.Samples <- (decostand(Spread.by.Sample[,-1], method = 'standardize', 1, na.rm = T))
rownames(Standardized.Samples) <- Spread.by.Sample$Metabolite.name
row.names.remove <- c("Thiamine monophosphate")

Standardized.Samples <- Standardized.Samples[!(row.names(Standardized.Samples) %in% row.names.remove), ]


# Clara for clusters ------------------------------------------------------
Average.Silh.Width = c()
k = c()
for(i in 2:15) {
  Clustered.List <- clara(Standardized.Samples, k = i, metric = "euclidean",
                         samples = nrow(Standardized.Samples))
  Average.Silh.Width = c(Average.Silh.Width, Clustered.List$silinfo$avg.width)
  k = c(k, i)
}
plot(k, Average.Silh.Width, type = "p", xlab = "No. clusters")


# Assign cluster numbers --------------------------------------------------
## why is this the best number? why 6? I replaced with 17 but not sure if thats right
best.number = which(Average.Silh.Width == (max(Average.Silh.Width))) + 1
low.number = which(Average.Silh.Width[1:7]==(max(Average.Silh.Width[1:7]))) + 1
if(best.number==low.number){
  best.number = which(Average.Silh.Width[1:14]==(max(Average.Silh.Width[7:14]))) + 1
}
High.Cluster.Number <- clara(Standardized.Samples, k = best.number, metric = "euclidean",
                                samples = nrow(Standardized.Samples), correct.d = FALSE)
Low.Cluster.Number <- clara(Standardized.Samples, k = low.number, metric = "euclidean", 
                               samples = nrow(Standardized.Samples), correct.d = FALSE)


# Combine cluster membership with normalized data -------------------------
Cluster.Membership <- data.frame(Mass.Feature = names(Low.Cluster.Number$clustering), 
                                   Cluster_low = Low.Cluster.Number$clustering) %>%
  full_join(., data.frame(Mass.Feature = names(High.Cluster.Number$clustering), 
                          Cluster_best = High.Cluster.Number$clustering))


# Adjust layout for heatmap plots -----------------------------------------
Heatmap.Data = Standardized.Samples %>%
  mutate(Mass.Feature = rownames(Standardized.Samples)) %>%
  full_join(Cluster.Membership) %>%
  gather(Sample.ID, Value, -Mass.Feature, -Cluster_low, -Cluster_best) %>%
  left_join(All.Data.Filtered %>%
              rename(treatment = Sample.ID,
                     Sample.ID = Replicate.Name) %>%
              dplyr::select(Sample.ID, Experiment, treatment) %>% 
              unique())  %>%
  unique() %>%
  arrange(Cluster_low) %>%
  select(-Sample.ID) %>%
  rename(Sample.ID = treatment) %>%
  mutate(Sample.ID = factor(Sample.ID, levels = ID.levels))


# Plot data ---------------------------------------------------------------
ggplot() + 
  geom_tile(data = Heatmap.Data, aes(x = Sample.ID, y = Mass.Feature, fill = Value)) +
  facet_grid(Cluster_best~., 
             scales = "free_y",
             space = "free_y") +
  scale_fill_viridis(option = "viridis")+
  theme(#axis.text.y  = element_blank(),
    axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
    strip.text = element_blank()) +
  ggtitle("Best Cluster Number: 11")


ggplot() + 
  geom_tile(data = Heatmap.Data, aes(x = Sample.ID, y = Mass.Feature, fill = Value)) +
  facet_grid(Cluster_low~., 
             scales = "free_y",
             space = "free_y") +
  scale_fill_viridis(option = "viridis")+
  theme(axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) +
  ggtitle("Low Cluster Number: 2")
