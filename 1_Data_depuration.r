#############################################################################################
##################  Script for data cleaning from Lopez-Reyes et al. 2023 ##################
#############################################################################################
##### --------------------------------------------------------------------------------- #####
#############################################################################################
#################################### LOAD THE PACKAGES ######################################
#############################################################################################
library(CoordinateCleaner)
library(sp)
library(rgdal)
library(dplyr)
library(raster)
library(occ)
library(dismo)
library(data.table)
library(maptools)
library(rgeos)
library(ENMeval)
library(stringr)
library(kuenm)
library(raster)
library(rangemap)
library(smoothr)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################# LOAD THE RECORDS AND LAYERS ###############################
#############################################################################################
var <- list.files("E:/7_Doctorado/Cap_1/2_Layers/No_cor",
                     pattern = ".tif",
                     full.names = TRUE)
stack <- raster::stack(var)
names(stack) = paste0("bio_", c(0, 10, 13, 14, 15, 2, 7), "_10m")

r <- raster(var[1])
r2 <- raster(var[2])

gbif <- read.csv("E:/7_Doctorado/Cap_1/1_Records/Raw/Sce_gbif.csv")
nat <- read.csv("E:/7_Doctorado/Cap_1/1_Records/Raw/Sce_nat.csv")
#############################################################################################
################################# CLEAN THE PRESENCE RECORDS ################################
#############################################################################################
##### --------------------------------------------------------------------------------- #####
############# GBIF #############
gbif_1 <- dplyr::select(gbif, speciesKey, decimalLatitude, decimalLongitude, species)
gbif_2 <- na.omit(gbif_1)
colnames(gbif_2) <- c("speciesKey", "Lat", "Long", "Species")
gbif_3 <- gbif_2[c("Species", "Long", "Lat", "speciesKey")]
gbif_env <- raster::extract(stack, gbif_3[,c("Long","Lat")])
gbif_f_occ_env <- data.frame(gbif_3, gbif_env)
gbif_occ_env <- na.omit(gbif_f_occ_env)
gbif_data_clean <- (gbif_occ_env[,c("Species","Long","Lat", "speciesKey")])
gbif_specieslist <- split(gbif_data_clean, gbif_data_clean$speciesKey)
gbif_allNames <- names(gbif_specieslist)
############# SAVE THE INDIVIDUAL RECORDS TO MERGE SUBSPECIES #############
for (i in 1:length(gbif_specieslist)){
setwd("E:/7_Doctorado/Cap_1/1_Records/Clean/Gbif")
spp <- gbif_specieslist[[i]]
spp$Species <- paste0(spp[2,1])
spp$Species <- gsub(" ", "_", spp$Species)
name <- paste0(spp[2,1])
name <- gsub(" ", "_", name)
write.csv(spp, paste0(name, ".csv"), row.names = FALSE)}
############# LOAD THE RECORDS AND COMBINE IN ONLY ONE #############
setwd("E:/7_Doctorado/Cap_1/1_Records/Clean/Gbif")
files_gbif <- list.files (pattern = ".csv")
temp_gbif <- lapply(files_gbif, fread, sep=",")
data_gbif <- rbindlist(temp_gbif, fill = TRUE)
data_gbif <- dplyr::select(data_gbif, Species, Long, Lat)
write.csv(data_gbif, file="all_gbif.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####
############# NATURALISTA #############
nat_1 <- dplyr::select(nat, latitude, longitude, scientific_name)
nat_2 <- na.omit(nat_1)
colnames(nat_2) <- c("Lat", "Long", "Species")
nat_3 <- nat_2[c("Species", "Long", "Lat")]
nat_env <- raster::extract(stack, nat_3[,c("Long","Lat")])
nat_f_occ_env <- data.frame(nat_3, nat_env)
nat_occ_env <- na.omit(nat_f_occ_env)
nat_data_clean <- (nat_occ_env[,c("Species","Long","Lat")])
nat_data_clean$Species <- gsub(" ", "_", nat_data_clean$Species)
nat_data_clean$Species <- word(nat_data_clean$Species, 1, 2, sep="_")
nat_specieslist <- split(nat_data_clean, nat_data_clean$Species)
nat_allNames <- names(nat_specieslist)
############# SAVE THE INDIVIDUAL RECORDS TO MERGE SUBSPECIES #############
for (i in 1:length(nat_allNames)){
setwd("E:/7_Doctorado/Cap_1/1_Records/Clean/Nat")
spp <- nat_specieslist[[i]]
spp$Species <- paste0(spp[2,1])
name <- paste0(spp[2,1])
name <- gsub(" ", "_", name)
write.csv(spp, paste0(name, ".csv"), row.names = FALSE)}
############# LOAD THE RECORDS AND COMBINE IN ONLY ONE #############
setwd("E:/7_Doctorado/Cap_1/1_Records/Clean/Nat")
files_nat <- list.files (pattern = ".csv")
temp_nat = lapply(files_nat, fread, sep=",")
data_nat <- rbindlist(temp_nat)
write.csv(data_nat,file="all_nat.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####
##### COMBINE AND CLEAN GBIF AND NATURALISTA RECORDS, THEN A SPLIT THEM PER SPECIES #####
all_nat <- read.csv("E:/7_Doctorado/Cap_1/1_Records/Clean/Nat/all_nat.csv")
all_gbif <- read.csv("E:/7_Doctorado/Cap_1/1_Records/Clean/gbif/all_gbif.csv")
joint <- rbind.data.frame(all_gbif, all_nat)
##### SPLIT THE MERGE RECORDS FOR EVERY SPECIES #####
all_specieslist <- split(joint, joint$Species)
all_allNames <- names(all_specieslist)
##### REMOVE DUPLICATE RECORDS AND GEOGRAPHIC OUTLIERS #####
for (i in 1:length(all_specieslist)){
all_occ <- all_specieslist[[i]]
all_xy <- all_occ[c(2:3)]
all_samp <- gridSample(all_xy, r, n = 1)
all_samp1 <- gridSample(all_samp, r2, n = 1)
all_samp1["Species"] <- NA
all_samp1$Species <- all_occ[1,1]
all_samp2 <- all_samp1[c("Species","Long","Lat")]
all_species_name <- all_occ[1,1]
all_genus <- substr(all_species_name, 1, 10)
all_epit <- substr(all_species_name, 12, nchar(all_species_name))
all_names <- paste0(all_genus, "_", all_epit)
clean <- cc_outl(all_samp2, lon = "Long", lat = "Lat", species = "Species", method = "quantile",
  mltpl = 2, value = "clean", sampling_thresh = 0, verbose = TRUE, min_occs = 3,
  thinning = FALSE)
##### SAVE SPECIES WITH MORE THAN TWO RECORDS (FOR THE HULLS MODELS) #####
setwd("E:/7_Doctorado/Cap_1/1_Records/Joint")
if (length(clean[,1]) > 2){
write.csv(clean, paste0(all_names, ".csv"), row.names = FALSE)}}
##### --------------------------------------------------------------------------------- #####
############# SAVE SPECIES WITH AT LEAST 5 RECORDS FOR ENM #############
setwd("E:/7_Doctorado/Cap_1/1_Records/Joint")
all_10 <- list.files(pattern = ".csv")
for (i in 1:length(all_10)){
setwd("E:/7_Doctorado/Cap_1/1_Records/Joint")
data <- read.csv(all_10[i])
name <- data[1,1]
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
if (length(data[,1]) > 4){
write.csv(data, paste0(name, ".csv"), row.names = FALSE)}}
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### END OF THE CODE ######################################
#############################################################################################
##### --------------------------------------------------------------------------------- #####










