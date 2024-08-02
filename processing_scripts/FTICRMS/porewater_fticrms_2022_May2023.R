## TEMPEST Porewater FTICR
##
## This is a data processing script for TEMPEST, a sub-project of the DOE-funded 
## COMPASS project (https://compass.pnnl.gov/). 
##
## This script imports and processes FTICR data obtained from Formularity 
## 
## Created: November 2022
## Kaizad F. Patel for EXCHANGE
#  Modified by AMP May 2024

# Notes from EMSL:
# Two of the supplicate runs were much lower than the original run. I would discard these r2 data files.
# EUP_60602_TMP_C_POOL_T4_r2_25Sep23_Fir_450SA
# EUP_60602_TMP_SW_SOURCE_HR2_r2_25Sep23_Fir_450SA

# Peak assignments for CHOS, CHOP, and CHON in certain areas of the formularity output have been recommended to be removed by EMSL staff
# need to filter the mass error > 0.5 and < -0.5 


## These functions have been modified from the `fticrrr` package and workflow: 
## https://github.com/kaizadp/fticrrr

################################################## #

# The functions for this script can be found in `Processing_Scripts/fticrrr-functions`
# See function description for more details

##############################
##############################


# 0. load packages --------------------------------------------------------
library(tidyverse)
library(googledrive)



# 1. SET input file paths -------------------------------
#REPORT = read_csv("../TEMPEST-1-porewater/data_do_not_commit/fticrms/ReportforCOMPASS_Nov2023/xtra_COMPASS_60602_2023_Report.csv")

REPORT_path = "https://drive.google.com/drive/folders/1hEmD6CnQlojRYug-3eMHpTFzzBso5-8c"

sample_list <- readxl::read_excel("../TEMPEST-1-porewater/sample lists/Myers-Pigg_TEMPEST_FT_samplelist_Sept2023.xlsx") %>%
  mutate(sample_name = stringr::str_replace(Sample_ID,"pooled","POOL"),
         sample_name = stringr::str_replace(Sample_ID,"Source","SOURCE"),
         Event = stringr::str_extract(Sample_ID, "TMP"),
         Plot = stringr::str_extract(Sample_ID, 'FW|SW|C|ESTUARY'), 
         Grid = stringr::str_extract(Sample_ID, "B4|C3|C6|D5|E3|F4|F6|H3|H6|I5|SOURCE|ESTUARY|POOL|WELL"),
         Timepoint = stringr::str_extract(Sample_ID,"T[0-9]|HR[0-9]"),
         Timepoint = case_when(Timepoint == "HR8" ~ "HR7", #change the estuary HR8 to HR7
                               TRUE ~ Timepoint),
         Pool_Timepoint = stringr::str_extract(sample_name,"[0-9]{8}_\\d{8}|[0-9]{8}-\\d{8}"),
         sample_name = stringr::str_remove(sample_name,"(?<=[0-9]{8})_\\d{8}|(?<=[0-9]{8})-\\d{8}"),
         sample_name = stringr::str_remove(sample_name, "_FW\\d{2}|_SW\\d{2}"))
  
  
import_data = function(directory){
  
  ## a. Create a list of files to download
  files <- 
    drive_ls(directory) %>% 
    filter(grepl("xtra_COMPASS_60602_2023_Report.csv", name))
  
  ## b. Download files to local (don't worry, we'll delete em in a sec)
  lapply(files$id, drive_download, overwrite = TRUE)
  
  dat <- read_csv(files$name)
  
  
  ## c. pull a list of file names
  ## then read all files and combine
  # 
  # filePaths <- files$id
  # dat <- 
  #   do.call(rbind, lapply(filePaths, function(path){
  #     # then add a new column `source` to denote the file name
  #     df <- read.csv(files$name)
  #     #  df <- read.delim(path, skip = 2)
  #     df[["source"]] <- rep(path, nrow(df))
  #     df}))
  # 
  ## d. delete the temporary files
  file.remove(c(files$name))  
  
  ## e. output
  dat
}

## import the raw data files
REPORT = import_data(REPORT_path) %>% bind_rows()

# 2. source the functions --------------------------------------------------------
source("../TEMPEST_Porewater/processing_scripts/FTICRMS/fticrrr_processing_scripts/a-functions_processing.R")

# EUP_60602_TMP_C_POOL_T4_r2_25Sep23_Fir_450SA
# EUP_60602_TMP_SW_SOURCE_HR2_r2_25Sep23_Fir_450SA

report = REPORT %>% select(-TMP_C_POOL_T4_r2_25S, -TMP_SW_SOURCE_HR2_r2_25S) %>%
  #remove the emsl worked up columns at the end:
  select(-(kmd:`cho-ind`), -...142, -class)

reps_suffix = "_(a_r1|b_r1|r1|a_r2|b_r2|r2|r3|2|3)$"

# fticr_meta = make_fticr_meta(report)$meta2

apply_filter_report(report)

fticr_meta = 
  fticr_report %>% 
  # select only the relevant columns for the formula assignments
  dplyr::select(Mass:Candidates) %>% 
  # select only necessary columns
  dplyr::select(Mass, C, H, O, N, S, P, El_comp, Class) %>% 
  # create columns for indices
  dplyr::mutate(AImod = round((1+C-(0.5*O)-S-(0.5*(N+P+H)))/(C-(0.5*O)-S-N-P),4),
                NOSC =  round(4-(((4*C)+H-(3*N)-(2*O)-(2*S))/C),4),
                DBE = round((C-0.5*H+0.5*N+1),0), 
                HC = round(H/C,2),
                OC = round(O/C,2)) %>% 
  # create column/s for formula
  # first, create columns for individual elements
  # then, combine
  dplyr::mutate(formula_c = if_else(C>0,paste0("C",C),as.character(NA)),
                formula_h = if_else(H>0,paste0("H",H),as.character(NA)),
                formula_o = if_else(O>0,paste0("O",O),as.character(NA)),
                formula_n = if_else(N>0,paste0("N",N),as.character(NA)),
                formula_s = if_else(S>0,paste0("S",S),as.character(NA)),
                formula_p = if_else(P>0,paste0("P",P),as.character(NA)),
                formula = paste0(formula_c,formula_h, formula_o, formula_n, formula_s, formula_p),
                formula = str_replace_all(formula,"NA","")) %>% 
  dplyr::select(Mass, formula, El_comp, Class, HC, OC, AImod, DBE, NOSC, C:P)


# subset metadata to include only the formula and columns with numeric data. Required to aggregate the data below. El_comp and Class are effectively removed here.
# Removes molecular formula where the AImod is not presented. If not, future calculations will be interrupted. 
fticr_meta_numeric = 
  fticr_meta %>% 
  dplyr::select(Mass, formula,C,H,O,N,S,P, NOSC, DBE, AImod,HC, OC)%>%
  filter(AImod != "-Inf")%>%
  filter(AImod != "Inf")

fticr_samples <- fticr_report%>% #Removes all data but the Mass and Sample data  
  select(-c(C:Candidates))

# Merges Sample data with filtered metadata. 
fticr_samples <- merge(x=fticr_meta_numeric, y=fticr_samples, by="Mass", all.x=TRUE)

# Simply reorders columns to put the molecular formula in the first column and remove columns other than the elemental formula and the sample data. 
fticr_samples <- fticr_samples%>%
  select(formula,everything())

#aggregate function to remove and merge duplicate rows
meta_agg <- aggregate(.~formula, data=fticr_meta_numeric, mean) #aggregates metadata and takes the mean of duplicate rows with the same molecular formula. 
#This ensures there are no changes in the metadata and only the merging of rows

fticr_samples_agg <- aggregate(.~ formula,data=fticr_samples, sum) #aggregates the metadata and samples to merge rows with the same molecular formula. We keep the metadata in here
#to ensure the metadata and samples are merged consistently. Duplicated molecular formula will have their intensities summed. This will however cause inaccuracies in the 
#metadata. For instance, a merged molecular formula with 6 carbons will show as 12 carbons. This metadata is removed on the next line. 

fticr_samples_agg <- fticr_samples_agg %>% #removing bad metadata to only include sample data in this data frame. 
  select(-c(Mass:OC))


#fticr_data = make_fticr_data(report)$data_samples_blank_corrected
#fticr_blanks = make_fticr_data_intensities(report)$data_blanks

#

# # 3.  Export processed data -----------------------------------------------
# fticr_data %>% write.csv("Data/Processed/EC1_Water_FTICR_L2_20230308.csv", row.names = FALSE)
# fticr_meta %>% write.csv("Data/Processed/EC1_Water_FTICR_meta_L2_20230308.csv", row.names = FALSE)

# 3. Export processed data --------------------------------------------------

L1directory = "https://drive.google.com/drive/folders/1yhukHvW4kCp6mN2jvcqmtq3XA5niKVR3"

fticr_data %>% write.csv("./ec1_water_fticrms_l1.csv", row.names = FALSE)
fticr_meta %>% write.csv("./ec1_water_fticrms_meta_l2.csv", row.names = FALSE)

L2directory = "https://drive.google.com/drive/folders/1M-ASGuRoKqswiKbUWylWzoAyUmMPm367"


drive_upload(media = "./ec1_water_fticrms_l1.csv", name= "ec1_water_fticrms_l1.csv", path = L1directory)
drive_upload(media = "./ec1_water_fticrms_meta_l2.csv", name= "ec1_water_fticrms_meta_l2.csv", path = L2directory)

file.remove("./ec1_water_fticrms_l1.csv")
file.remove("./ec1_water_fticrms_meta_l2.csv")

# 4. Check with Metadata for missing:

import_l1_data = function(directory){
  
  ## a. Create a list of files to download
  files <- 
    drive_ls(directory) %>% 
    filter(name == "ec1_water_fticrms_L1.csv")
  
  ## b. Download files to local (don't worry, we'll delete em in a sec)
  lapply(files$id, drive_download, overwrite = TRUE)
  
  dat <- read.csv(files$name)
  
  ## d. delete the temporary files
  file.remove(c(files$name))  
  
  ## e. output
  dat
}



## import the L2 data files
fticr_data_l1 = import_l1_data(L1directory)

fticr_sample_list = colnames(fticr_data_l1)

fticr_sample_kits = as_tibble(fticr_sample_list) %>% 
  filter(stringr::str_detect(value, "K[0-9]{3}")) %>%
  mutate(kit_id= stringr::str_extract(value, "K[0-9]{3}"),
         campaign = "EC1",
         transect_location = "water",
         data_collected = TRUE) %>%
  dplyr::select(campaign, kit_id, transect_location, data_collected)


source("./Processing_Scripts/Metadata_kit_list.R")

metadata_collected %>%
  filter(sample_method == "bottle_1l") -> meta_filter

fticr_sample_kits %>%
  full_join(meta_filter, by = c("campaign", "kit_id", "transect_location")) %>%
  mutate(data_collected = case_when(is.na(data_collected) ~ FALSE,
                                    TRUE ~ data_collected),
         note = case_when(kit_id %in% c("K001","K007") ~ "kit compromised",
                          kit_id %in% c("K027") ~ "sample compromised",
                          kit_id %in% c("K044") ~ "sample not sent for analysis")
         ) %>%
  rename(sample_analyzed = data_collected) %>%
  select(campaign, kit_id, transect_location,sample_analyzed, note) %>%
  arrange(kit_id)-> fticr_sample_list

View(fticr_sample_list)

#save sample list
fticr_sample_list %>% write.csv("./ec1_water_fticrms_sample_list.csv", row.names = FALSE)
drive_upload(media = "./ec1_water_fticrms_sample_list.csv", name= "ec1_water_fticrms_sample_list.csv", path = L2directory)
file.remove("./ec1_water_fticrms_sample_list.csv")

#save L2 fticrms data 
fticr_data %>% select(-K001, -BLANK) %>% write.csv("./ec1_water_fticrms_l2.csv", row.names = FALSE)
drive_upload(media = "./ec1_water_fticrms_l2.csv", name= "ec1_water_fticrms_L2.csv", path = L2directory)
file.remove("./ec1_water_fticrms_l2.csv")

