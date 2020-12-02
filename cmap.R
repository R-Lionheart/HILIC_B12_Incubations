## Testing

library(tidyverse)


mydata <- read.csv("Ingalls_B12Incubations_HILIC - data.csv") 

Internal.Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                               stringsAsFactors = FALSE, header = TRUE) %>%
  select(Compound.Name, Compound.Name_SQL) %>%
  filter(Compound.Name %in% mydata$compound) %>%
  rename(compound = Compound.Name)

together <- mydata %>%
  left_join(Internal.Standards, by = "compound") %>%
  mutate(compound = Compound.Name_SQL) %>%
  select(-Compound.Name_SQL) %>%
  mutate(lon = ifelse(lon == 25, -158.5, -157)) %>%
  mutate(value = replace_na(value, ""))


  
write.csv(together, "~/Downloads/HILIC_B12_Incubations_FixedCmpdNames.csv")