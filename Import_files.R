source("B12_Functions.R")

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


# Rearrange datasets ------------------------------------------------------

# Positive
Area.pos <- RearrangeDatasets(Area_HILIC.POS_B12.Inc, parameter = "Area.Value")
Mz.pos   <- RearrangeDatasets(Mz_HILIC.POS_B12.Inc, parameter = "Mz.Value")
RT.pos   <- RearrangeDatasets(RT_HILIC.POS_B12.Inc, parameter = "RT.Value")
SN.pos   <- RearrangeDatasets(SN_HILIC.POS_B12.Inc, parameter = "SN.Value")


# Negative
Area.neg <- RearrangeDatasets(Area_HILIC.NEG_B12.Inc, parameter = "Area.Value")
Mz.neg   <- RearrangeDatasets(Mz_HILIC.NEG_B12.Inc, parameter = "Mz.Value")
RT.neg   <- RearrangeDatasets(RT_HILIC.NEG_B12.Inc, parameter = "RT.Value")
SN.neg   <- RearrangeDatasets(SN_HILIC.NEG_B12.Inc, parameter = "SN.Value")


# Combine to one dataset --------------------------------------------------
combined.pos <- Area.pos %>%
  left_join(Mz.pos) %>%
  left_join(SN.pos) %>%
  left_join(RT.pos) %>%
  mutate(Column = "HILICPos") %>%
  select(Replicate.Name, Column, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

combined.neg <- Area.neg %>%
  left_join(Mz.neg) %>%
  left_join(SN.neg) %>%
  left_join(RT.neg) %>%
  mutate(Column = "HILICNeg") %>%
  select(Replicate.Name, Column, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

combined <- rbind(combined.pos, combined.neg) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name)) 

combined$Replicate.Name <- gsub("^.{0,1}", "", combined$Replicate.Name)


currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_combined_", currentDate, ".csv", sep = "")

write.csv(combined, csvFileName, row.names = FALSE)

rm(list = ls())
