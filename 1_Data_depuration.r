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
######################## CREATE THE CALIBRATION AND TRANSFERING REGIONS #####################
#############################################################################################
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
data <- list.files (pattern = ".csv")
############# CALIBRATION (M) #############
for(i in 1:length(data)){
spp <- read.csv(data[[i]], header=T, sep=",")
spp_name <- substr(data[i],1,nchar(data[i])-4)
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/M1")
chull1 <- rangemap_hull(occurrences = spp, hull_type = "convex", buffer_distance = 50000, 
						extent_of_occurrence = FALSE, area_of_occupancy = FALSE,
						save_shp = TRUE, name = spp_name,
						overwrite = TRUE)
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
next}
#### BEFORE YOU CAN RUN THE NEX CODE YOU NEED TO REMOVE THE UNIQUE RECORDS FROM THE M1 FOLDER ####
mydir <- "E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/M1/"
delfiles <- dir(path = mydir, pattern="unique_records")
file.remove(file.path(mydir, delfiles)) 
############# CORRECT THE HULLS COLUMNS AND REMOVE ISLANDS  #############
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/M1")
chull <-list.files(pattern=".shp")
for(i in 1:length(chull)){
mod_1 <- shapefile(chull[[i]])
name <- substr(chull[i],1,nchar(chull[i])-4)
mod_1@data
mod_1$species = name
names(mod_1@data)[names(mod_1@data)=="species"] <- "BINOMIAL"
mod_1$areakm2 <- area(mod_1)
max <- max(mod_1$areakm2)/2
mod_1 <- drop_crumbs(mod_1, threshold = max, drop_empty = TRUE)
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/M")
writeOGR(mod_1, dsn=getwd(), layer = name, driver="ESRI Shapefile")
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/M1")}
############# AFTER THIS YOU CAN DELETE THE M1 FOLDER #############
##### --------------------------------------------------------------------------------- #####
############# PROJECTION (G) #############
jjm <- readOGR("E:/7_Doctorado/Cap_1/2_Layers/Eco", "eco_NA")
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
press <- list.files(pattern = ".csv")
for(i in 1:length(press)){
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
occ <- read.csv(press[[i]])
occ1<-occ[c(2,3)]
occ2 <- SpatialPoints(occ1)
plot(jjm, axes=TRUE)
points(occ2, col="red")
CRS.new<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")  # este es el que tiene el shape de WWF
proj4string(occ2)<-CRS.new
jjm_subset <- jjm[occ2, ]
plot(jjm_subset, axes=TRUE)
points(occ2, col="red")
name <- occ[1,1]
name <- gsub(" ", "_", name)
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases"))
writeOGR(jjm_subset, "G", paste(press[[i]]), driver="ESRI Shapefile")
next}
##### --------------------------------------------------------------------------------- #####
#############################################################################################
############################## CREATING DATA FOR KUENM (SWD) ################################
#############################################################################################
############# CREATE THE FOLDERS FOR EACH SPECIES #############
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
sp <- list.files(pattern = ".csv")
##### CREATE ONE FOLDER PER SPECIES #####
for (i in 1:length(sp)){
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
data <- read.csv(sp[i])
name <- data[1,1]
name <- gsub(" ", "_", name)
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species")
dir.create(name)
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name))
dir.create("Background")
dir.create("M_var")
dir.create("G_var")
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name, "/M_var"))
dir.create("Set_1")
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name, "/G_var"))
dir.create("Set_1")
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name, "/G_var/Set_1"))
dir.create("current")
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name))}
##### --------------------------------------------------------------------------------- #####
############# CROP THE LAYERS FOR M AND G REGIONS #############
############# MS #############
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/M")
M <-list.files(pattern=".shp")
datum <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
for (i in 1:length(M)){
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/M")
Ms <-shapefile(M[[i]])
Ms@proj4string<- datum
cor_C <- crop(stack,Ms)
mas_C <- mask(cor_C,Ms)
name <- substr(M[i],1,nchar(M[i])-4)
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name, "/M_var/Set_1"))
name_var <- names(stack)
writeRaster(mas_C, filename=paste0(name_var), bylayer = T, format = "ascii")
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/M")
next}
##### --------------------------------------------------------------------------------- #####
############# GS #############
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/G")
G <-list.files(pattern=".shp")
datum <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
for (i in 1:length(G)){
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/G")
Gs <-shapefile(G[[i]])
Gs@proj4string<- datum
cor_C <- crop(stack, Gs)
mas_C <- mask(cor_C, Gs)
name <- substr(G[i], 1, nchar(G[i])-8)
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name, "/G_var/Set_1/current"))
name_var <- names(stack)
writeRaster(mas_C, filename = paste0(name_var), bylayer = TRUE, format = "ascii")
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/G")
next}
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################## PREPARE RECORDS FOR SWD ##################################
#############################################################################################
############# JOINT #############
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
all_10 <- list.files(pattern = ".csv")
for (i in 1:length(all_10)){
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
data <- read.csv(all_10[i])
name <- data[1,1]
path <- paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name)
setwd(path)

path <- paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name, "/M_var/Set_1")
var <- list.files(path,
					pattern = ".asc",
					full.names = TRUE)
					
stack <- raster::stack(var)

prep <- prepare_swd(data, species = "Species", longitude = "Long",
                    latitude = "Lat", raster.layers = stack,
                    sample.size = 10000)

if (length(data[,1]) > 4){
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name))
write.csv(prep$joint, paste0("occ_joint.csv"), row.names = FALSE)
setwd(paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name, "/Background"))
write.csv(prep$background, paste0("Set_1.csv"), row.names = FALSE)}}
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### TEST AND TRAIN #######################################
#############################################################################################
phts <- list.files("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species",
					full.names = TRUE)
base <- list.files("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint",
					pattern = ".csv")

var <- list.files("E:/7_Doctorado/Cap_1/2_Layers/No_cor",
                     pattern = ".tif",
                     full.names = TRUE)
stack <- raster::stack(var)
names(stack) = paste0("bio_", c(0, 10, 13, 14, 15, 2, 7), "_10m")

for (i in 1:length(base)){
g <- read.csv(paste0(phts[i], "/occ_joint.csv"))
name <- g[1,1]
name <- gsub(" ", "_", name)
if (length(g[,1]) > 25){
name <- substr(base[i],1,nchar(base[i])-4)
occ <- g[c(2:3)]
lay_swd <- g[c(4:10)]
bg.coords <- randomPoints(r, 10000)
chk1.pts <- get.checkerboard1(occ, stack@layers, bg.coords, aggregation.factor = 2, gridSampleN = 10000)
occ["chk1"] <- chk1.pts$occs.grp
occ <- cbind(occ, g)
occ <- occ[c("Species", "Long", "Lat", "chk1", "bio_0_10m", "bio_10_10m", "bio_13_10m",
"bio_14_10m", "bio_15_10m", "bio_2_10m", "bio_7_10m")]
ck1 <- occ[occ$chk1=="1",]
ck2 <- occ[occ$chk1=="2",]
ck1 <- ck1[c("Species", "Long", "Lat", "chk1", "bio_0_10m", "bio_10_10m", "bio_13_10m",
"bio_14_10m", "bio_15_10m", "bio_2_10m", "bio_7_10m")]
ck2 <- ck2[c("Species", "Long", "Lat", "chk1", "bio_0_10m", "bio_10_10m", "bio_13_10m",
"bio_14_10m", "bio_15_10m", "bio_2_10m", "bio_7_10m")]
path <- paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name)
setwd(path)
write.csv(ck2,paste0("occ_test.csv"), row.names = F)
write.csv(ck1,paste0("occ_train.csv"), row.names = F)
}
else{
name <- substr(base[i],1,nchar(base[i])-4)
samp <- sample(nrow(g), round(0.6 * nrow(g))) 
ck2 <- g[samp,]
ck1 <- g[-samp,]
path <- paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", name)
setwd(path)
write.csv(ck1, paste0("occ_test.csv"), row.names = F)
write.csv(ck2, paste0("occ_train.csv"), row.names = F)}}
##### --------------------------------------------------------------------------------- #####
##### END OF THE CODE #####










