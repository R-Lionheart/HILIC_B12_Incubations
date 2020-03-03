library(cluster)
library(ISLR) # for college dataset
library(progress)
library(readr)
library(Rtsne)
library(tidyverse)
options(scipen = 999)

BMISd <- read.csv("data_processed/BMIS_Output_2020-02-20.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature, Adjusted.Area, Run.Cmpd) %>%
  filter(!str_detect(Run.Cmpd, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m")) %>%
  separate(Run.Cmpd, sep = " ", into = c("Replicate.Name"), remove = FALSE) %>%
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
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3")) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), fill = "right", remove = FALSE) %>%
  select(Mass.Feature, SampID, Adjusted.Area) %>%
  drop_na()
standards <- read.csv("data_extras/Ingalls_Lab_Standards.csv", stringsAsFactors = FALSE) %>%
  rename(Mass.Feature = Compound.Name) 

cluster.test <- BMISd %>%
  left_join(standards) %>%
  replace(is.na(.), "Unknown") %>%
  select(Mass.Feature:Column, m.z) %>%
  mutate(Mass.Feature = as.factor(Mass.Feature),
         SampID = as.factor(SampID),
         Compound.Type = as.factor(Compound.Type),
         Column = as.factor(Column)) %>%
  select(-m.z)

# Second attempt ----------------------------------------------------------
# original
college_clean <- College %>%
  mutate(name = row.names(.),
         accept_rate = Accept/Apps,
         isElite = cut(Top10perc,
                       breaks = c(0, 50, 100),
                       labels = c("Not Elite", "Elite"),
                       include.lowest = TRUE)) %>%
  mutate(isElite = factor(isElite)) %>%
  select(name, accept_rate, Outstate, Enroll,
         Grad.Rate, Private, isElite)
glimpse(college_clean)

# mine
metab_clean <- cluster.test %>%
  mutate(isControl = as.factor(ifelse(str_detect(SampID, "Control"), "Control", "NonControl")),
         eddy = as.factor(ifelse(str_detect(SampID, "IL1"), "Cyclonic", "Anticyclonic"))) 
glimpse(metab_clean)

# original
gower_dist <- daisy(college_clean[, -1], # Remove college name
                    metric = "gower",
                    type = list(logratio = 3))

mygower_dist <- daisy(metab_clean[, -1], # Remove mass.feature
                    metric = "gower",
                    type = list(logratio = 3))

# Check attributes to ensure the correct methods are being used
# (I = interval, N = nominal)
# Note that despite logratio being called, 
# the type remains coded as "I"
# original
summary(gower_dist)
gower_mat <- as.matrix(gower_dist)

# mine
summary(mygower_dist)
mygower_mat <- as.matrix(mygower_dist)

# Output most similar pair
# original
college_clean[
  which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
        arr.ind = TRUE)[1, ], ]
# Most dissimilar pair
college_clean[
  which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
        arr.ind = TRUE)[1, ], ]

# mine
#similar (can this be expanded to include more grouping?)
metab_clean[
  which(mygower_mat == min(mygower_mat[mygower_mat != min(mygower_mat)]),
        arr.ind = TRUE)[1, ], ]
# dissimilar
metab_clean[
  which(mygower_mat == max(mygower_mat[mygower_mat != max(mygower_mat)]),
        arr.ind = TRUE)[1, ], ]

# Calculate silhouette width for many k using PAM, and plot sihouette width (higher is better)
# original
sil_width <- c(NA)
for(i in 2:10){
  pam_fit <- pam(gower_dist,
                 diss = TRUE,
                 k = i)
  sil_width[i] <- pam_fit$silinfo$avg.width
}
plot(1:10, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, sil_width)

# mine
mysil_width <- c(NA)
for(i in 2:10){
  mypam_fit <- pam(mygower_dist,
                 diss = TRUE,
                 k = i)
  mysil_width[i] <- mypam_fit$silinfo$avg.width
}
plot(1:10, mysil_width,
     xlab = "My Number of clusters",
     ylab = "My Silhouette Width")
lines(1:10, mysil_width)


# original
pam_fit <- pam(gower_dist, diss = TRUE, k = 3)
pam_results <- college_clean %>%
  dplyr::select(-name) %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary
college_clean[pam_fit$medoids, ]

# mine
mypam_fit <- pam(mygower_dist, diss = TRUE, k = 4)
mypam_results <- metab_clean %>%
  dplyr::select(-Mass.Feature) %>%
  mutate(cluster = mypam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
mypam_results$the_summary
metab_clean[mypam_fit$medoids, ]

# original
tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)
tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = college_clean$name)

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))

# mine
mytsne_obj <- Rtsne(mygower_dist, is_distance = TRUE)
mytsne_data <- mytsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(mypam_fit$clustering),
         name = metab_clean$Mass.Feature)

ggplot(aes(x = X, y = Y), data = mytsne_data) +
  geom_point(aes(color = cluster))




# First attempt -----------------------------------------------------------
# # Compute Gower distance
# gower_dist <- daisy(cluster.test, metric = "gower")
# gower_mat <- as.matrix(gower_dist)
# # Print most similar clients
# cluster.test[which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]), arr.ind = TRUE)[1, ], ]
# 
# # Print most dissimilar clients
# cluster.test[which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]), arr.ind = TRUE)[1, ], ]
# 
# sil_width <- c(NA)
# for(i in 2:13) {  
#   pam_fit <- pam(gower_dist, diss = TRUE, k = i)  
# sil_width[i] <- pam_fit$silinfo$avg.width  
# }
# 
# plot(1:13, sil_width,
#      xlab = "Number of clusters",
#      ylab = "Silhouette Width")
# lines(1:13, sil_width)
# 
# k <- 13
# pam_fit <- pam(gower_dist, diss = TRUE, k)
# pam_results <- cluster.test %>%
#   mutate(cluster = pam_fit$clustering) %>%
#   group_by(cluster) %>%
#   do(the_summary = summary(.))
# pam_results$the_summary
# 
# 
# tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)
# tsne_data <- tsne_obj$Y %>%
#   data.frame() %>%
#   setNames(c("X", "Y")) %>%
#   mutate(cluster = factor(pam_fit$clustering))
# ggplot(aes(x = X, y = Y), data = tsne_data) +
#   geom_point(aes(color = cluster))


## K means clustering test 
ggplot(cluster.test, aes(x = Mass.Feature, y = Adjusted.Area)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

set.seed(20)
myCluster <- kmeans(cluster.test[, 4], centers = 4)
myCluster 

myCluster$cluster <- as.factor(myCluster$cluster)

ggplot(cluster.test, aes(x = Mass.Feature, y = Adjusted.Area, color = myCluster$cluster)) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point()


ggplot(cluster.test, aes(x = Mass.Feature, y = SampID, color = myCluster$cluster)) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_jitter(shape = 15,
              position = position_jitter(0.21))


ggplot(data = cluster.test, aes(x = Adjusted.Area, y = myCluster$cluster, fill = SampID)) +
  scale_y_discrete(breaks = seq(1, 7, by = 1)) +
  geom_tile() +
  coord_equal() +
  theme_classic()