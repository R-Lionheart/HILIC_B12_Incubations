# Actual script goes here

source("B12_Inc_Functions.R")

# Import all MSDial files --------------------------------------------------

filenames <- RemoveCsv(list.files(path = 'data_raw', pattern = '*.csv'))

for (i in filenames) {
  filepath <- file.path('data_raw', paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

matching.variable <- "hilic"


columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 
                     'SMILES', 'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 
                     'RT.similarity', 'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')

# Set header, filter unknowns ---------------------------------------

runs <- grep(matching.variable, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))


headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs

for (df in seq_along(headers.set)) { 
  headers.set[[df]] <- headers.set[[df]] %>% filter(!Metabolite.name == "Unknown")
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
}

# Change variable classes -------------------------------------------------

classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs


list2env(classes.changed, globalenv())


######
# START HERE #
######




# Rearrange datasets ------------------------------------------------------

# Positive
SN_HILIC.POS_B12.Inc <- SN_HILIC.POS_B12.Inc %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "SN.Value",
    starts_with("X")) %>%
  select(Replicate.Name, SN.Value, everything())

RT_HILIC.POS_B12.Inc <- RT_HILIC.POS_B12.Inc %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "RT.Value",
    starts_with("X")) %>%
  select(Replicate.Name, RT.Value, everything())

Area_HILIC.POS_B12.Inc <- Area_HILIC.POS_B12.Inc %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "Area.Value",
    starts_with("X")) %>%
  select(Replicate.Name, Area.Value, everything())

Mz_HILIC.POS_B12.Inc <- Mz_HILIC.POS_B12.Inc %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "MZ.Value",
    starts_with("X")) %>%
  select(Replicate.Name, MZ.Value, everything())

# Negative
SN_HILIC.NEG_B12.Inc <- SN_HILIC.NEG_B12.Inc %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "SN.Value",
    starts_with("X")) %>%
  select(Replicate.Name, SN.Value, everything())

RT_HILIC.NEG_B12.Inc <- RT_HILIC.NEG_B12.Inc %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "RT.Value",
    starts_with("X")) %>%
  select(Replicate.Name, RT.Value, everything())

Area_HILIC.NEG_B12.Inc <- Area_HILIC.NEG_B12.Inc %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "Area.Value",
    starts_with("X")) %>%
  select(Replicate.Name, Area.Value, everything())

Mz_HILIC.NEG_B12.Inc <- Mz_HILIC.NEG_B12.Inc %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "MZ.Value",
    starts_with("X")) %>%
  select(Replicate.Name, MZ.Value, everything())

# Combine to one dataset --------------------------------------------------
combined.pos <- Area_HILIC.POS_B12.Inc %>%
  left_join(Mz_HILIC.POS_B12.Inc) %>%
  left_join(SN_HILIC.POS_B12.Inc) %>%
  left_join(RT_HILIC.POS_B12.Inc) %>%
  mutate(Column = "HILICPos") %>%
  select(Replicate.Name, Column, Area.Value, MZ.Value, RT.Value, SN.Value, everything())

combined.neg <- Area_HILIC.NEG_B12.Inc %>%
  left_join(Mz_HILIC.NEG_B12.Inc) %>%
  left_join(SN_HILIC.NEG_B12.Inc) %>%
  left_join(RT_HILIC.NEG_B12.Inc) %>%
  mutate(Column = "HILICNeg") %>%
  select(Replicate.Name, Column, Area.Value, MZ.Value, RT.Value, SN.Value, everything())

combined <- rbind(combined.pos, combined.neg) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name)) 

combined$Replicate.Name <- gsub("^.{0,1}", "", combined$Replicate.Name)


currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_combined_", currentDate, ".csv", sep = "")

write.csv(combined, csvFileName, row.names = FALSE)

rm(list = ls())
