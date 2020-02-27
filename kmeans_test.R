# Retention Time Table ----------------------------------------------------
RT.Table <- size.fraction_0.2 %>%
  filter(str_detect(Replicate.Name, "_Std_")) %>%
  filter(str_detect(Metabolite.Name, "_Std")) %>%
  group_by(Metabolite.Name) %>%
  mutate(Mean.RT.Value = mean(RT.Value, na.rm = TRUE),
         Min.RT.Value = min(RT.Value, na.rm = TRUE),
         Max.RT.Value = max(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.Diff = RT.Value - RT.Expected) %>% 
  mutate(RT.Diff.abs = abs(RT.Value - RT.Expected)) %>%
  select(Replicate.Name, Metabolite.Name, RT.Expected, RT.Value, Mean.RT.Value:RT.Diff.abs)
# mutate(Midrange.RT.Diff = (min(RT.Diff) + max(RT.Diff)) / 2) %>%
# mutate(High.Low = ifelse(RT.Diff > Midrange.RT.Diff, "High", "Low")) %>%
# select(Replicate.Name, Metabolite.Name, RT.Expected, RT.Value, Mean.RT.Value:High.Low)


## K means clustering test 
ggplot(RT.Table, aes(RT.Value, Replicate.Name, color = Metabolite.Name)) + 
  geom_point() +
  ggtitle("Standard Retention Time Differences")
ggplot(RT.Table, aes(RT.Diff, Replicate.Name, color = Metabolite.Name)) + 
  geom_point() +
  ggtitle("Expected vs Real Retention Time Differences")

cluster.test <- RT.Table %>%
  arrange(Metabolite.Name)

set.seed(20)
RTCluster <- kmeans(cluster.test[, 8], 2, nstart = 20)
RTCluster 

RTCluster$cluster <- as.factor(RTCluster$cluster)
ggplot(cluster.test, aes(RT.Diff, Replicate.Name, color = RTCluster$cluster)) + 
  geom_point() +
  ggtitle("K-means clustering: RT Value Differences")

cluster.test$cluster <- RTCluster$cluster

# RT differences plot
RT.Table.clustered <- cluster.test %>%
  group_by(Metabolite.Name) %>%
  mutate(High.Low = as.character(ifelse(cluster == 1, "Low", "High"))) %>%
  select(-cluster) %>%
  # TESTING AREA #
  unique() 

RT.Plot <- ggplot(RT.Table.clustered, aes(x = Replicate.Name, y = RT.Diff, fill = High.Low)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~Metabolite.Name, scales = "fixed") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 5),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Retention Time Differences")
print(RT.Plot)

################################################################################
## Constructing Tolerances

Tolerance.Table.High <- RT.Table.clustered %>%
  filter(High.Low == "High") %>%
  unique() %>%
  ## TESTING
  select(Metabolite.Name, RT.Expected, RT.Value) %>%
  #filter(Metabolite.Name == "FA 16:0_Std") %>%
  group_by(Metabolite.Name) %>%
  mutate(Ave.High = mean(RT.Value)) %>%
  mutate(Ave.High.Diff = abs(RT.Expected - Ave.High)) %>%
  select(-RT.Value) %>%
  unique()

Tolerance.Table.Low <- RT.Table.clustered %>%
  filter(High.Low == "Low") %>%
  unique() %>%
  ## TESTING
  select(Metabolite.Name, RT.Expected, RT.Value) %>%
  #filter(Metabolite.Name == "FA 16:0_Std") %>%
  group_by(Metabolite.Name) %>%
  mutate(Ave.Low= mean(RT.Value)) %>%
  mutate(Ave.Low.Diff = abs(RT.Expected - Ave.Low)) %>%
  select(-RT.Value) %>%
  unique()



# Tolerance production ----------------------------------------------------
## in progress
Full.Tolerance.Table <- size.fraction_0.2 %>%
  select(Metabolite.Name, Replicate.Name, RT.Expected, RT.Value) %>%
  #filter(Metabolite.Name == "FA 16:0_Std") %>% 
  filter(str_detect(Metabolite.Name, "_Std")) %>%
  # TESTING AREA #
  left_join(Tolerance.Table.High) %>%
  left_join(Tolerance.Table.Low) %>%
  select(Metabolite.Name, RT.Expected, Ave.High:Ave.Low.Diff) %>%
  unique()

FA16_HighTolerance = unique(Full.Tolerance.Table$Ave.High.Diff)
FA16_LowTolerance = unique(Full.Tolerance.Table$Ave.Low.Diff)
