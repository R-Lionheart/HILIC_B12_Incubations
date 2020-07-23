
# Original online example ------------------------------------------------
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

# original
gower_dist <- daisy(college_clean[, -1], # Remove college name
                    metric = "gower",
                    type = list(logratio = 3))

# Check attributes to ensure the correct methods are being used
# (I = interval, N = nominal)
# Note that despite logratio being called,
# the type remains coded as "I"
# original
summary(gower_dist)
gower_mat <- as.matrix(gower_dist)

# Output most similar pair
# original
college_clean[
  which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
        arr.ind = TRUE)[1, ], ]
# Most dissimilar pair
college_clean[
  which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
        arr.ind = TRUE)[1, ], ]

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

# original
pam_fit <- pam(gower_dist, diss = TRUE, k = 3)
pam_results <- college_clean %>%
  dplyr::select(-name) %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary
college_clean[pam_fit$medoids, ]

# original
tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)
tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = college_clean$name)

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))


## K means clustering test
ggplot(metab_clean, aes(x = Mass.Feature, y = Adjusted.Area)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))


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

# -----------------------------------------------------------
# How to compute CLARA (Clustering Large Applications) in R
# -----------------------------------------------------------

# generate Simulated Data
# Generate 500 objects, divided into 2 clusters.
df <- rbind(cbind(rnorm(200,0,8), rnorm(200,0,8)),
            cbind(rnorm(300,50,8), rnorm(300,50,8)),
            cbind(rnorm(400,100,8), rnorm(400,100,8)))
colnames(df) <- c("x", "y")
rownames(df) <- paste0("S", 1:nrow(df))

# optimal number of clusters
fviz_nbclust(df, clara, method = "silhouette")+
  theme_classic()

# Compute CLARA
clara.res <- clara(df, 3, samples = 50, pamLike = TRUE)

# Print components of clara.res
print(clara.res)

# Add clustering result to the Data
dd <- cbind(df, cluster = clara.res$cluster)
head(dd, n = 4)

# Visualise clusters
fviz_cluster(clara.res,
             palette = c("#00AFBB", "#FC4E07", "#E7B800"), # color palette
             ellipse.type = "t", # Concentration ellipse
             geom = "point", pointsize = 1,
             ggtheme = theme_classic()
)

# -----------------------------------------------------------
# My version
# -----------------------------------------------------------
clara.test <- cluster.test %>%
  select(Adjusted.Area)

# optimal number of clusters
fviz_nbclust(clara.test, clara, method = "silhouette")+
  theme_classic()

# Compute CLARA
myclara.res <- clara(clara.test, 2, samples = 50, pamLike = TRUE)

# Print components of clara.res
print(clara.res)

# Add clustering result to the Data
mydd <- cbind(clara.test, cluster = myclara.res$cluster)
head(mydd, n = 4)

# Visualise clusters
fviz_cluster(myclara.res,
             palette = c("#00AFBB", "#FC4E07"), # color palette
             ellipse.type = "t", # Concentration ellipse
             geom = "point", pointsize = 1,
             ggtheme = theme_classic()
)

# -----------------------------------------------------------
# KRH version
# -----------------------------------------------------------
datclu.clara <- clara(Iso_wide, k = 3, metric = "euclidean", samples = 100, sampsize = nrow(Iso_wide))

fviz_cluster(datclu.clara,
             ellipse.type = "t", # Concentration ellipse
             geom = "point", pointsize = 1,
             ggtheme = theme_classic()
)