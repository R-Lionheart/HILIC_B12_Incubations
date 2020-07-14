## Osmolyes

source("src/B12_Functions.R")

Raw.Osmolytes <- read.csv("data_extras/Method_Osmolytes.csv", stringsAsFactors = FALSE) %>%
  rename(Osmolytes = 2,
         Possible = 6) %>%
  mutate_all(list(~na_if(.,"")))

Possible <- Raw.Osmolytes %>% 
  select(Possible) %>%
  drop_na() %>%
  slice(-1) %>%
  rbind(Raw.Osmolytes[28:33, 2])

Osmolytes <- Raw.Osmolytes %>%
  select(Osmolytes) %>%
  slice(4:24) 

Osmolytes[nrow(Osmolytes)+7,] <- NA

All.Osmolytes <- cbind(Osmolytes, Possible)

All.Osmolytes$Osmolytes <- capFirst(All.Osmolytes$Osmolytes)
All.Osmolytes$Possible <- capFirst(All.Osmolytes$Possible)

All.Osmolytes[All.Osmolytes == "NANA"] <- NA

write.csv(All.Osmolytes, "data_processed/Osmolytes_edited.csv")