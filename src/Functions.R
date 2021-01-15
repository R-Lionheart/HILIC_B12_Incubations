## Function definitions ##

library(RCurl)
library(rlist)
library(tidyverse)

CapitalizeFirstLetter <- function(s) {
  # Convert lowercase first letter of string to uppercase.
  #
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

ChangeClasses <- function(df, start.column) {
  # Change specified columns from factors to numeric.
  #
  # Args
  #   df: MSDial dataframe containing sample columns.
  #   start.column: 
  #
  # Returns
  #   df: MSDial dataframe, with specified columns changed to a numeric class. 
  for (i in c(start.column:ncol(df))) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

ChangeXClasses <- function(df) {
  # Identifies columns starting with X and changes their class to numeric.
  #
  # Args
  #   df: MSDial dataframe reorganized to drop all empty rows at the top.
  #
  # Returns
  #   df: MSDial dataframed with modified sample column classes.
  #
  col.test <- grepl("^X", names(df))
  for (i in which(col.test == TRUE)) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

FixBMISNames <- function(df) {
  # Shortcut for separating Replicate.Name into manageable columns.
  #
  # Args
  #   df: Long format dataframe containing Replicate.Name column.
  #
  # Returns
  #   df.Replicate.Name.split: df with column split.
  #
  df.Replicate.Name.split <- df %>%
    separate(Replicate.Name, into = c("one", "two", "SampID", "replicate")) %>%
    select(-c("one", "two")) %>%
    unite(SampID, replicate, col = "Replicate.Name")
  
  return(df.Replicate.Name.split)
}

IdentifyDuplicates <- function(df) {
  # Determine which compounds are detected in both positive and negative HILIC runs.
  # 
  # Args
  #   df: MSDial or Skyline dataframe in long form.
  # 
  # Returns
  #   duplicates: Simple dataframe of listed compounds that have been identified as duplicates.
  #
  duplicates <- df %>%
    group_by(Metabolite.Name, Replicate.Name) %>%
    mutate(number = 1) %>%
    mutate(ticker = cumsum(number)) %>%
    filter(ticker == 2) %>%
    ungroup() %>%
    select(Metabolite.Name) %>%
    unique()
  return(duplicates)
}

IdentifyRunTypes <- function(df) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   df: Raw output file from Skyline, or other dataframe containing a Replicate.Name
  #       column.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(df$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  print(paste("Your runtypes are:", toString(unique(run.type))))
}

MakeLong <- function(df) {
  # Shortcut for pivoting multiple datasets from wide to long format. 
  #
  # Args
  #   df: Standardized dataframe, mass features as column names and samples as rownames.
  #
  # Returns
  #   df.long: Dataframe, still standardized, but pivoted to long format.
  #
  colnames(df)[1] <- "Replicate.Name"
  df.long <- pivot_longer(df, -Replicate.Name,
                          names_to = "Mass.Feature",
                          values_to = "Area.BMISd.Normd") %>%
    select(Mass.Feature, Replicate.Name, Area.BMISd.Normd)
  
  return(df.long)
} 

MakeWide <- function(df, area.column.name) {
  # Shortcut for pivoting multiple dataframes from long to wide format.
  #
  # Args
  #   df: Long dataframe, containing a Replicate.Name and Area column).
  #   area.column.name: String of the column name that contains information
  #                     pertaining to area.
  #
  # Returns
  #   df.rownames: Original dataframe pivoted wider, with the mass feature column
  #                changed to rownames
  #
  df.wide <- df %>%
    ungroup() %>%
    pivot_wider(names_from = Replicate.Name,
                values_from = all_of(area.column.name)) %>%
    as.data.frame()
  
  df.rownames <- df.wide[,-1]
  rownames(df.rownames) <- df.wide[,1]
  df.rownames[is.na(df.rownames)] <- NA
  
  return(df.rownames)
}

RearrangeDatasets <- function(df, parameter) {
  # Shortcut for altering multiple datasets using the tidyr::gather() function.
  #
  # Args
  #   df: MSDial dataframe with first (n) empty rows removed.
  #   parameter: Table value. This parameter will become the column name when 
  #              changed to long format.
  #
  # Returns
  #   df.long: MSDial dataframe, changed to long format and with a custom-named value column.
  df.long <- df %>%
    pivot_longer(names_to = "Replicate.Name",
                 values_to = "parameter",
                 starts_with("X")) %>%
    select(Replicate.Name, parameter, everything())
  
  names(df.long)[2] <- parameter
  df.long$Replicate.Name <- gsub("^.{0,1}", "", df.long$Replicate.Name)
  
  return(df.long)
}

RemoveCsv <- function(full.filepaths) {
  # Remove a .csv file extension and obtain basename from a given filepath or list
  # of filepaths.
  #
  # Args
  #   Character string of filepaths in a directory.
  #
  # Returns
  #   Character string of file basenames, without a csv extension.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

SetHeader <- function(df) {
  # Remove empty or unnecessary lines from machine output, and make column names headers.
  # Also change the default header "Metabolite.name" to "Metabolite.Name" so as to maintain
  # consistency across analyses.
  #
  # Args
  #   df: Raw output file from MSDial.
  #
  # Returns
  #   df: Modified dataframe with correct headers and no empty lines.
  #
  df <- df[!(is.na(df[1]) | df[1]==""), ]
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  names(df)[names(df) == "Metabolite.name"] <- "Metabolite.Name"
  
  
  return(df)
}

StandardizeMetabolites <- function(df) {
  # Remove any "Ingalls_" prefixes that may be present in a dataframe.
  # Remove "X" prefixes in syntactically correct Replicate Names.
  #
  # Args
  #   df: MSDial dataframe.
  #
  # Returns
  #   df.standardized: Dataframe with above modifications.
  #
  df.standardized <- df %>%
    mutate(Metabolite.Name = ifelse(str_detect(Metabolite.Name, "Ingalls_"), 
                                    sapply(strsplit(Metabolite.Name, "_"), `[`, 2), Metabolite.Name)) 
  
  #df.standardized$Replicate.Name <- gsub("^.{0,1}", "", df.standardized$Replicate.Name)
  
  return(df.standardized)
}


TrimWhitespace <- function (x) gsub("^\\s+|\\s+$", "", x)

UploadFiles <- function(myfilepath) {
  # Shortcut for uploading multiple long-format files.
  # Args
  #   myfilepath: inner-project filepath to dataframe.
  #
  # Returns
  #   uploaded.df: dataframe with relevant columns selected.
  #
  uploaded.df <- read.csv(myfilepath, stringsAsFactors = FALSE) %>%
    select(Mass.Feature:Adjusted.Area)
  
  return(uploaded.df)
}

# FindStdDev <- function(df) {
#   df.first <- df %>%
#     group_by(Mass.Feature) %>%
#     group_split()
#   df.midframe <- lapply(df.first, function(x) mutate(x, Std.dev = sd(Area.Ave, na.rm = TRUE)))
#   df.final <- bind_rows(df.midframe)
#   
#   return(df.final)
# }

# MakeBarPlot <- function(df, title) {
#   df.plot <- ggplot(df, aes(x = Mass.Feature, y = Area.Ave, fill = SampID)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
#           axis.text.y = element_text(size = 10),
#           legend.position = "top",
#           axis.ticks.length=unit(.25, "cm"),
#           strip.text = element_text(size = 10)) +
#     geom_errorbar(aes(ymin=Area.Ave - Std.dev, ymax=Area.Ave + Std.dev), width=.2,
#                   position=position_dodge(.9)) +
#     ggtitle(title)
#   print(df.plot)
# }

# MakeFacetGraphs <- function (df, title, scale) {
#   df.plot <- ggplot(df, aes(x = SampID, y = Area.Ave, fill = SampID)) +
#     geom_bar(stat = "identity") +
#     facet_wrap( ~Mass.Feature, scales = scale) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
#           axis.text.y = element_text(size = 5),
#           legend.position = "top",
#           axis.ticks.length=unit(.25, "cm"),
#           strip.text = element_text(size = 10)) +
#     ggtitle(title)
#   print(df.plot)
# }

