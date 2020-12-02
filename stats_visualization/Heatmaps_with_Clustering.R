source("src/B12_Functions.R")
source("src/biostats.R")
source("src/coldiss.R")

library(cluster)
library(tidyverse)
library(vegan)

# User input
file.pattern <- "Time0"
my.filename <- "All.Data.Raw"

ID.levels <- c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT",    
               "IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT",
               "IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um", "IL1DMBnoBT5um", "IL1noBT5um",
               "IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um", "IL2DMBnoBT5um", "IL2noBT5um")


filenames <- RemoveCsv(list.files(path = "data_processed", 
                                  pattern = regex(file.pattern, ignore_case = TRUE)))
for (i in filenames) {
  filepath <- file.path("data_processed", paste(i, ".csv", sep = ""))
  assign(my.filename, read.csv(filepath, stringsAsFactors = FALSE))
}


# Tidy and filter data ----------------------------------------------------
All.Data.Filtered <- All.Data.Raw %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), remove = FALSE) %>%
  select(-c("one", "two", "four")) %>%
  mutate(Experiment = ifelse(str_detect(SampID, "IL1"), "Cyclonic", "Anticyclonic")) %>%
  mutate(Adjusted.Area = ifelse(is.na(Adjusted.Area), 1000, Adjusted.Area))


# Pivot dataframe to wide -------------------------------------------------
# Make a dataframe where each row is a sample and each column is a MF
Wide.All.Data.Filtered <- All.Data.Filtered %>%
  group_by(Mass.Feature) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Mass.Feature, values_from = Adjusted.Area)  %>%
  select(-row)

Spread.by.Sample <- Wide.All.Data.Filtered %>%
  column_to_rownames(var = "Replicate.Name") %>%
  #select(-SampID, -Experiment) %>%
  t() %>%
  as.data.frame()


# Standardize wide dataset ------------------------------------------------
Standardized.Samples <- (decostand(Spread.by.Sample[,-1], method = 'standardize', 1, na.rm = T))


# Clara for clusters ------------------------------------------------------
Average.Silh.Width = c()
k = c()
for(i in 2:15) {
  Clustered.List <- clara(Standardized.Samples, k = i, metric = "euclidean",
                         samples = nrow(Standardized.Samples))
  Average.Silh.Width = c(Average.Silh.Width, Clustered.List$silinfo$avg.width)
  k = c(k, i)
}
plot(k, Average.Silh.Width, type = "b", xlab = "Number of clusters")


# Assign cluster numbers --------------------------------------------------
best.number = which(Average.Silh.Width == (max(Average.Silh.Width))) + 1
low.number = which(Average.Silh.Width[7:14] == (max(Average.Silh.Width[7:14]))) + 1
if(best.number == low.number){
  best.number = which(Average.Silh.Width[1:14] == (max(Average.Silh.Width[7:14]))) + 1
} # If the same, make the best number to a bigger option. In this case, 8 is the best choice above 7.
# Maybe pick 4, or 6 instead.
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
  gather(SampID, Value, -Mass.Feature, -Cluster_low, -Cluster_best) %>%
  left_join(All.Data.Filtered %>%
              rename(treatment = SampID,
                     SampID = Replicate.Name) %>%
              dplyr::select(SampID, Experiment, treatment) %>% 
              unique())  %>%
  unique() %>%
  arrange(Cluster_low) %>%
  select(-SampID) %>%
  rename(SampID = treatment,
         Standardized.Area.Value = Value) %>%
  mutate(SampID = factor(SampID, levels = ID.levels))


# Plot data ---------------------------------------------------------------
ggplot() + 
  geom_tile(data = Heatmap.Data, aes(x = SampID, y = Mass.Feature, fill = Standardized.Area.Value)) +
  facet_grid(Cluster_best~., 
             scales = "free_y",
             space = "free_y") +
  theme(#axis.text.y  = element_blank(),
    axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
    strip.text = element_blank()) +
  ggtitle(paste("Best Cluster Number:", best.number))


ggplot() + 
  geom_tile(data = Heatmap.Data, aes(x = SampID, y = Mass.Feature, fill = Standardized.Area.Value)) +
  facet_grid(Cluster_low~., 
             scales = "free_y",
             space = "free_y") +
  theme(#axis.text.y  = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0, hjust = 0),
        strip.text = element_blank()) +
  ggtitle(paste("Low Cluster Number:", low.number))
