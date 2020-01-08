## Function definitions ##

# TODO --------------------------------------------------------
# Organize directory referencing a little easier? Set value ahead of time?
# Nested loop or list comprehension for renaming the uploads. Functions?
# Order function descriptions neatly at beginning of script.
# Cleanly name items in the upload, rather than manually changing it to SN, RZ, etc.
# Figure out a better way to choose columns out of MSDIAL
# Add comments to clarify the process.
# Make function for always downloading the most up to date Ingalls lab standards.
# Figure out a way to preserve the QC parameter values.
# Fix the StandardizeVariables function

library(ggplot2)
library(rlist)
library(stringr)
library(tidyverse)
library(tidyr)
options(scipen=999)

SetHeader <- function(df) {
  # Remove empty or unnecessary lines from machine output, and make column names headers.
  #
  # Args
  #   df: Raw output file from MSDial.
  #
  # Returns
  #   df: modified dataframe with correct headers and no empty lines.
  #
  df <- df[!(is.na(df[1]) | df[1]==""), ]
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  
  return(df)
}

RemoveCsv <- function(full.filepaths) {
  # Remove a .csv file extension and obtain basename from a given list of filepaths.
  #
  # Args
  #   Character strings of filepaths in a directory.
  #
  # Returns
  #   Character strings of file basenames, without a csv extension.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

ChangeClasses <- function(df) {
  # Change specified columns from factors to numeric.
  #
  # Args
  #   df: MSDial dataframe containing sample columns.
  #
  # Returns
  #   df: MSDial dataframe, with specified columns changed to a numeric class. 
  for (i in c(10:ncol(df))) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

IdentifyRunTypes <- function(msdial.file) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   msdial.file: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(msdial.file$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  print(paste("Your runtypes are:", toString(unique(run.type))))
}

TrimWhitespace <- function (x) gsub("^\\s+|\\s+$", "", x)

IdentifyDuplicates <- function(df) {
  # Determine which compounds are detected in both positive and negative HILIC runs.
  # 
  # Args
  #   df: MSDial dataframe, containing all required parameters (MZ, SN, Area, etc),
  #       and modified to long form instead of wide.
  # 
  # Returns
  #   duplicates: Simple dataframe of listed compounds that have been identified as duplicates.
  #
  duplicates <- df %>%
    group_by(Metabolite.name, Replicate.Name) %>%
    mutate(number = 1) %>%
    mutate(ticker = cumsum(number)) %>%
    filter(ticker == 2) %>%
    ungroup() %>%
    select(Metabolite.name) %>%
    unique()
  return(duplicates)
}


# Functions --------------------------------------------

FindStdDev <- function(df) {
  df.first <- df %>%
    group_by(Mass.Feature) %>%
    group_split()
  df.midframe <- lapply(df.first, function(x) mutate(x, Std.dev = sd(Area.Ave, na.rm = TRUE)))
  df.final <- bind_rows(df.midframe)
  
  return(df.final)
}
MakeBarPlot <- function(df, title) {
  df.plot <- ggplot(df, aes(x = Mass.Feature, y = Area.Ave, fill = SampID)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "top",
          axis.ticks.length=unit(.25, "cm"),
          strip.text = element_text(size = 10)) +
    geom_errorbar(aes(ymin=Area.Ave - Std.dev, ymax=Area.Ave + Std.dev), width=.2,
                  position=position_dodge(.9)) +
    ggtitle(title)
  print(df.plot)
}
MakeFacetGraphs <- function (df, title, scale) {
  df.plot <- ggplot(df, aes(x = SampID, y = Area.Ave, fill = SampID)) +
    geom_bar(stat = "identity") +
    facet_wrap( ~Mass.Feature, scales = scale) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
          axis.text.y = element_text(size = 5),
          legend.position = "top",
          axis.ticks.length=unit(.25, "cm"),
          strip.text = element_text(size = 10)) +
    ggtitle(title)
  print(df.plot)
}