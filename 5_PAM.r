#############################################################################################
##################################### LOAD THE PACKAGES #####################################
#############################################################################################
library(raster)
library(rgdal)
library(sp)
library(wesanderson)
library(dplyr)
library(viridis)
library(letsR)
library(sf)
library(ape)
library(caper)
library(ggplot2)
library(patchwork)
library(data.table)
library(ggpubr)
library(climateStability)
library(tibble)
library(terra)
#############################################################################################
################################# LOAD THE NECESSARY DATA ###################################
#############################################################################################
##### --------------------------------------------------------------------------------- #####
#############################################################################################
#################################### ENVIRONMENTAL LAYERS ###################################
#############################################################################################
setwd("G:/7_Doctorado/Cap_1/2_Layers/Pred_layer")
vars <- list.files("G:/7_Doctorado/Cap_1/2_Layers/Pred_layer",
						pattern = ".tif",
						full.names = TRUE)
stck <- stack(vars)
names(stck) = paste0(c("CCV", "Elevation", "CEH_min", "CEH_max", "CVH"))
##### I NEED TO PROJECT THE ENVIRONMENTAL LAYERS IN THE SAME PROJECTION OF THE PAM #####
stck_proj <- projectRaster(stck, crs = CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"))
##### --------------------------------------------------------------------------------- #####
#############################################################################################
######################################### PHYLOGENIES #######################################
#############################################################################################
leache <- read.nexus("C:/Users/Kevin/Downloads/Filogenias/leache.nex")
tonini <- read.nexus("C:/Users/Kevin/Downloads/Filogenias/tonini.nex")
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################# MODELS (POLYGONS) #########################################
#############################################################################################
setwd("G:/7_Doctorado/Cap_1/3_Models/ahull")
ahull <-list.files(pattern=".shp")
setwd("G:/7_Doctorado/Cap_1/3_Models/chull")
chull <-list.files(pattern=".shp")
setwd("G:/7_Doctorado/Cap_1/3_Models/kuenm")
sdm <-list.files(pattern = ".shp")
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################# TRANSFORM THE POLYGONS ####################################
#############################################################################################
##### AHULL ##### 
setwd("G:/7_Doctorado/Cap_1/3_Models/ahull")
mod_ahull <- lapply(ahull, st_read) ##### READ THE SHAPEFILES #####
com_ahull <- do.call(what = sf:::rbind.sf, args = mod_ahull) ##### COMBINE THEM #####
mult_ahull <- st_cast(com_ahull, "MULTIPOLYGON") ##### CREATE A MULTIPOLYGON OBJECT #####
st_crs(mult_ahull) = 4326 ##### ASSINGS A COORDINATE SYSTEM (NOT PROJECTED YET) #####
dist_ahull <- as_Spatial(mult_ahull, cast = TRUE, ##### MAKING IT SPATIAL #####
						 IDs = paste0("BINOMIAL", seq_along(mult_ahull)))
dist_ahull <- spTransform(dist_ahull, CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m")) ##### PROJECT THE OBJECT #####		 
##### CHULL ##### 
setwd("G:/7_Doctorado/Cap_1/3_Models/chull")
mod_chull <- lapply(chull, st_read) ##### READ THE SHAPEFILES #####
com_chull <- do.call(what = sf:::rbind.sf, args = mod_chull) ##### COMBINE THEM #####
mult_chull <- st_cast(com_chull, "MULTIPOLYGON") ##### CREATE A MULTIPOLYGON OBJECT #####
st_crs(mult_chull) = 4326 ##### ASSINGS A COORDINATE SYSTEM (NOT PROJECTED YET) #####
dist_chull <- as_Spatial(mult_chull, cast = TRUE, ##### MAKING IT SPATIAL #####
						 IDs = paste0("BINOMIAL", seq_along(mult_chull)))
dist_chull <- spTransform(dist_chull, CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m")) ##### PROJECT THE OBJECT #####	
#####  SDM ##### 
setwd("G:/7_Doctorado/Cap_1/3_Models/kuenm")
mod_sdm <- lapply(sdm, st_read) ##### READ THE SHAPEFILES #####
com_sdm <- do.call(what = sf:::rbind.sf, args = mod_sdm) ##### COMBINE THEM #####
mult_sdm <- st_cast(com_sdm, "MULTIPOLYGON") ##### CREATE A MULTIPOLYGON OBJECT #####
st_crs(mult_sdm) = 4326 ##### ASSINGS A COORDINATE SYSTEM (NOT PROJECTED YET) #####
dist_sdm <- as_Spatial(mult_sdm, cast = TRUE, ##### MAKING IT SPATIAL #####
						 IDs = paste0("BINOMIAL", seq_along(mult_sdm)))
dist_sdm <- spTransform(dist_sdm, CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m")) ##### PROJECT THE OBJECT #####	
##### --------------------------------------------------------------------------------- #####						 
#############################################################################################
################################# CONSTRUCT THE PAMS ########################################
#############################################################################################
##### GET THE RESOLUTION AND THE EXTENT (NECCESARY FOR THE "lets.presab" function) #####
ahull_bbox <- bbox(dist_ahull)
ahull_long_min <- ahull_bbox[1]
ahull_long_max <- ahull_bbox[3]
ahull_lat_min <- ahull_bbox[2]
ahull_lat_max <- ahull_bbox[4]

chull_bbox <- bbox(dist_chull)
chull_long_min <- chull_bbox[1]
chull_long_max <- chull_bbox[3]
chull_lat_min <- chull_bbox[2]
chull_lat_max <- chull_bbox[4]

sdm_bbox <- bbox(dist_sdm)
sdm_long_min <- sdm_bbox[1]
sdm_long_max <- sdm_bbox[3]
sdm_lat_min <- sdm_bbox[2]
sdm_lat_max <- sdm_bbox[4]
##### CONSTRUCT THE PAM FOR ALL THE METHODS (ALPHA-HULL, CONVEXHULL, SDM) ##### 
PAM_ahull <- lets.presab(dist_ahull, count = T, xmn = ahull_long_min, xmx = ahull_long_max, ymn = ahull_lat_min, 
						ymx = ahull_lat_max, resol = 55500, crs = "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m", 
						crs.grid = "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m", remove.cells = TRUE)
PAM_chull <- lets.presab(dist_chull, count = T, xmn = chull_long_min, xmx = chull_long_max, ymn = chull_lat_min, 
						ymx = chull_lat_max, resol = 55500, crs = "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m", 
						crs.grid = "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m", remove.cells = TRUE)
PAM_sdm <- lets.presab(dist_sdm, count = T, xmn = sdm_long_min, xmx = sdm_long_max, ymn = sdm_lat_min, 
						ymx = sdm_lat_max, resol = 55500, crs = "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m", 
						crs.grid = "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m", remove.cells = TRUE)					  
####### BECAUSE THE PAM IS PROJECTED, WE NEED TO MATCH THE LAYERS WITH EXTENT AND RESOLUTION OF THE PAM #######
Richness_mask <- crop(PAM_ahull$Richness_Raster, dist_ahull)
Richness_mask <- mask(PAM_ahull$Richness_Raster, dist_ahull)
stck_proj1 <- crop(stck_proj, Richness_mask)
stck_proj1 <- resample(stck_proj1, Richness_mask, method = "bilinear")
stck_proj1 <- mask(stck_proj1, Richness_mask)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################### DATA NECESSARY FOR THE SPECIES APPROACH ANALYSIS #####################
#############################################################################################
####### RANGE SIZE PER SPECIES #######
setwd("G:/7_Doctorado/Cap_1/4_data/RS")
####### AHULL #######
RS_ahull <- lets.rangesize(PAM_ahull, units = "squaremeter")
RS_ahull <- RS_ahull / 1000000000000
RS_ahull1 <- enframe(RS_ahull, name = "Species", value = "Range_size")
RS_ahull1 <- as.data.frame(RS_ahull1)
write.csv(RS_ahull1, "rs_ahull.csv", row.names = FALSE)
RS_ahull1$Range_size <- log(RS_ahull1$Range_size)
write.csv(RS_ahull1, "rs_ahull_log.csv", row.names = FALSE)
####### CHULL #######
RS_chull <- lets.rangesize(PAM_chull, units = "squaremeter")
RS_chull <- RS_chull / 1000000000000
RS_chull1 <- enframe(RS_chull, name = "Species", value = "Range_size")
RS_chull1 <- as.data.frame(RS_chull1)
write.csv(RS_chull1, "rs_chull.csv", row.names = FALSE)
RS_chull1$Range_size <- log(RS_chull1$Range_size)
write.csv(RS_chull1, "rs_chull_log.csv", row.names = FALSE)
####### SDM #######
RS_sdm <- lets.rangesize(PAM_sdm, units = "squaremeter")
RS_sdm <- RS_sdm / 1000000000000
RS_sdm1 <- enframe(RS_sdm, name = "Species", value = "Range_size")
RS_sdm1 <- as.data.frame(RS_sdm1)
write.csv(RS_sdm1, "rs_sdm.csv", row.names = FALSE)
RS_sdm1$Range_size <- log(RS_sdm1$Range_size)
write.csv(RS_sdm1, "rs_sdm_log.csv", row.names = FALSE)
####### MIDPOINTS PER SPECIES #######
####### AHULL #######
mp_ahull <- lets.midpoint(PAM_ahull, planar = FALSE)
colnames(mp_ahull) <- c("Species", "Long", "Lat")
rsmp_ahull <- cbind(mp_ahull, RS_ahull1[,2])
####### CHULL #######
mp_chull <- lets.midpoint(PAM_chull, planar = FALSE)
colnames(mp_chull) <- c("Species", "Long", "Lat")
rsmp_chull <- cbind(mp_chull, RS_chull1[,2])
####### SDM #######
mp_sdm <- lets.midpoint(PAM_sdm, planar = FALSE)
colnames(mp_sdm) <- c("Species", "Long", "Lat")
rsmp_sdm <- cbind(mp_sdm, RS_sdm1[,2])
####### SAVE MY INDIVIDUALS RESULTS #######
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(rsmp_ahull, "rs_mp_ahull.csv", row.names = FALSE)
write.csv(rsmp_chull, "rs_mp_chull.csv", row.names = FALSE)
write.csv(rsmp_sdm, "rs_mp_sdm.csv", row.names = FALSE)
####### CREATE A SINGLE CSV WITH THE RANGE SIZE AND MIDPOINTS VALUES FOR ALL THE METHODS #######
all_mprs <- merge(rsmp_ahull, rsmp_chull, by = "Species", all = T)
all_met_RS <- merge(all_mprs, rsmp_sdm, by = "Species", all = T)
all_met_RS <- all_met_RS[complete.cases(all_met_RS), ]
colnames(all_met_RS) <- c("Species", "Long_ahull",  "Lat_ahull", "RS_ahull", "Long_chull",  "Lat_chull", "RS_chull", 
                          "Long_sdm",  "Lat_sdm", "RS_sdm")
write.csv(all_met_RS, "all_met_RS", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
######## I NEED TO COMPARE THE PREDICTIONS BETWEEN METHODS (RANGE SIZE AND MIDPOINT) ######## 
#############################################################################################
##### RANGE SIZE VS RANGE ZISE FOR ALL METHODS #####
all_met_RS$RS_ahull <- (all_met_RS$RS_ahull - min(all_met_RS$RS_ahull)) / (max(all_met_RS$RS_ahull) - min(all_met_RS$RS_ahull))
all_met_RS$RS_chull <- (all_met_RS$RS_chull - min(all_met_RS$RS_chull)) / (max(all_met_RS$RS_chull) - min(all_met_RS$RS_chull))
all_met_RS$RS_sdm <- (all_met_RS$RS_sdm - min(all_met_RS$RS_sdm)) / (max(all_met_RS$RS_sdm) - min(all_met_RS$RS_sdm))

a_c_hull_rs <- ggplot(all_met_RS, aes(RS_ahull, RS_chull)) +
  geom_point(color = "black", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Range Size - Alpha Hull", y = "Range Size - Convex Hull") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("number of species = 103") +  geom_abline(intercept = 0, slope = 1) +
  theme_classic()  

a_sdm_rs <- ggplot(all_met_RS, aes(RS_ahull, RS_sdm)) +
  geom_point(color = "black", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Range Size - Alpha Hull", y = "Range Size - SDM") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("number of species = 95") +  geom_abline(intercept = 0, slope = 1) +
  theme_classic()  

c_sdm_rs <- ggplot(all_met_RS, aes(RS_chull, RS_sdm)) +
  geom_point(color = "black", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Range Size - Convex Hull", y = "Range Size - SDM") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("number of species = 95") +  geom_abline(intercept = 0, slope = 1) +
  theme_classic()
##### MIDPOINTS VS MIDPOINTS FOR ALL METHODS #####
all_met_RS$Lat_ahull <- (all_met_RS$Lat_ahull - min(all_met_RS$Lat_ahull)) / (max(all_met_RS$Lat_ahull) - min(all_met_RS$Lat_ahull))
all_met_RS$Lat_chull <- (all_met_RS$Lat_chull - min(all_met_RS$Lat_chull)) / (max(all_met_RS$Lat_chull) - min(all_met_RS$Lat_chull))
all_met_RS$Lat_sdm <- (all_met_RS$Lat_sdm - min(all_met_RS$Lat_sdm)) / (max(all_met_RS$Lat_sdm) - min(all_met_RS$Lat_sdm))

a_c_hull_mp <- ggplot(all_met_RS, aes(Lat_ahull, Lat_chull)) +
  geom_point(color = "black", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Midpoints - Alpha Hull", y = "Midpoints - Convex Hull") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("") +  geom_abline(intercept = 0, slope = 1) +
  theme_classic()  

a_sdm_mp <- ggplot(all_met_RS, aes(Lat_ahull, Lat_sdm)) +
  geom_point(color = "black", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Midpoints - Alpha Hull", y = "Midpoints - SDM") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("") +  geom_abline(intercept = 0, slope = 1) +
  theme_classic()  

c_sdm_mp <- ggplot(all_met_RS, aes(Lat_chull, Lat_sdm)) +
  geom_point(color = "black", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Midpoints - Convex Hull", y = "Midpoints - SDM") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("") +  geom_abline(intercept = 0, slope = 1) +
  theme_classic()
ggarrange(a_c_hull_rs, a_sdm_rs, c_sdm_rs, a_c_hull_mp, a_sdm_mp, c_sdm_mp, ncol = 3, nrow = 2)
setwd("G:/7_Doctorado/Cap_1/5_figures")
ggsave("Fig_1.tiff", width = 10, height = 6, dpi = 1000) 
##### --------------------------------------------------------------------------------- #####
#############################################################################################
### THE NEXT ANALYSIS ARE PERFORMED ONLY WITH ALPHAHULL (METHOD WITH LESS OVERPREDICTION) ###
#############################################################################################
##### --------------------------------------------------------------------------------- #####
#############################################################################################
## CALCULATE THE MEAN VALUE FOR EVERY VARIABLE INSIDE THE GEOGRAPHIC RANGE OF EACH SPECIES ##
#############################################################################################
setwd("G:/7_Doctorado/Cap_1/3_Models/ahull")
dist <- list.files(pattern=".shp")
##### CLIMATE CHANGE VELOCITY (CCV) EXPRESSED IN METERS PER YEAR #####
ccv_list <- list()
for (i in 1:length(dist_ahull)){
ccv <- stck_proj1[[1]]
pol_dist <- shapefile(dist[i])
pol_dist <- spTransform(pol_dist, CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"))
ccv_dist <- crop(ccv, pol_dist)
ccv_dist <- mask(ccv_dist, pol_dist)
val_ccv <- cellStats(ccv_dist, "mean")
ccv_list[i] <- val_ccv}
ahull_CCV <- unlist(ccv_list)
##### ELEVATION #####
Elev_list <- list()
for (i in 1:length(dist_ahull)){
elev <- stck_proj1[[2]]
pol_dist <- shapefile(dist[i])
pol_dist <- spTransform(pol_dist, CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"))
elev_dist <- crop(elev, pol_dist)
elev_dist <- mask(elev, pol_dist)
val_elev <- cellStats(elev_dist, "mean")
Elev_list[i] <- val_elev}
ahull_elev <- unlist(Elev_list)
##### MINIMUN AND MAXIMUN MONTLY TEMPERATURE FOR THE CLIMATE EXTREMES HYPOTHESIS (CHE) #####
temp_max_list <- list()
temp_min_list <- list()
for (i in 1:length(dist_ahull)){
temp_min <- stck_proj1[[3]]
temp_max <- stck_proj1[[4]]
pol_dist <- shapefile(dist[i])
pol_dist <- spTransform(pol_dist, CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"))
temp_min_dist <- crop(temp_min, pol_dist)
temp_min_dist <- mask(temp_min, pol_dist)
temp_max_dist <- crop(temp_max, pol_dist)
temp_max_dist <- mask(temp_max, pol_dist)
val_temp_min <- cellStats(temp_min_dist, "mean")
val_temp_max <- cellStats(temp_max_dist, "mean")
temp_min_list[i] <- val_temp_min
temp_max_list[i] <- val_temp_max}
temp_min_list <- unlist(temp_min_list)
temp_max_list <- unlist(temp_max_list)
ahull_CEH <- cbind(temp_min_list, temp_max_list)
##### ANUAL TEMPERATURE RANGE FOR THE CLIMATE VARIABILITY HYPOTHESIS (CVH) #####
temp_rg_list <- list()
for (i in 1:length(dist_ahull)){
temp_rng <- stck_proj1[[5]]
pol_dist <- shapefile(dist[i])
pol_dist <- spTransform(pol_dist, CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"))
temp_rng_dist <- crop(temp_rng, pol_dist)
temp_rng_dist <- mask(temp_rng_dist, pol_dist)
val_temp_rg <- cellStats(temp_rng_dist, "mean")
temp_rg_list[i] <- val_temp_rg}
ahull_CVH <- unlist(temp_rg_list)
##### THEN I COMBINE THE LIST OF EACH LAYER WITH THE RANGE SIZE AND COORDINATES OF EACH SPECIES #####
sp_vars_all <- cbind(rsmp_ahull, ahull_CCV, ahull_elev, ahull_CEH, ahull_CVH)
colnames(sp_vars_all) <- c("Species", "Long", "Lat", "Range_size", "CCV", "Elevation", "CEH_min", "CEH_max",
                           "CVH")
##### SAVE MY RESULTS #####
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(sp_vars_all, "all_data_specieslevel.csv", row.names = FALSE) 
###### SPECIES LEVEL #####
##### CREATE THE COMPARATIVE DATA #####
comp_ton <- comparative.data(tonini, data = sp_vars_all, names.col = "Species", vcv = T)
comp_lea <- comparative.data(leache, data = sp_vars_all, names.col = "Species", vcv = T)
####### GET THE VARIABLES FOR THE SPECIES LEVEL FOR EACH DATA SET (ALL, TONINI, LEACHE) #######
##### GET THE NAMES OF THE SPECIES MATCHING WITH THE PHYLOGENIES #####
ton_names <- as.data.frame(comp_ton$phy$tip.label)
colnames(ton_names) <- c("Species")
lea_names <- as.data.frame(comp_lea$phy$tip.label)
colnames(lea_names) <- c("Species")
##### GET THE VALUES ONLY FOR THE SPECIES MATCHING WITH THE PHYLOGENIES #####
sp_vars_all <- sp_vars_all
sp_vars_ton <- merge(ton_names, sp_vars_all, by = "Species", all = FALSE)
sp_vars_lea <- merge(lea_names, sp_vars_all, by = "Species", all = FALSE)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(sp_vars_ton, "ton_data_specieslevel.csv", row.names = FALSE)
write.csv(sp_vars_lea, "lea_data_specieslevel.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################# ASSEMBLAGES APPROACH ######################################
#############################################################################################
##### CALCULATES RANGE SIZE PER SITE (MEDIAN) #####
RS_sites <- lets.maplizer(PAM_ahull, log(RS_ahull), rownames(RS_ahull), func = median, ras = TRUE) ## log 
RS_matrix <- RS_sites$Matrix

ap_df <- as.data.frame(ap)
##### PLOT RASTER ######
library(ggplot2)
library(reshape2)
library(tidyverse)
world <- map_data("world")
world

ggplot(ap_df, aes(x, y, fill = layer)) +
  geom_tile(colour="black") +
  scale_fill_gradientn(colours = hcl.colors(3, palette = "Fall", 0.8)) +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.1
  )  + 
  theme_bw()

##### ADD THE ENVIRONMENTAL VARIABLES TO THE PAM #####
vars_PAM <- lets.addvar(PAM_ahull, stck_proj1, fun = median, onlyvar = FALSE)
PAM_onlyvars <- vars_PAM[,106:110] ##### SELECT THE COLUMNS ONLY WITH ENVIRONMENTAL LAYERS #####
##### JOINT THE RS WITH THE FULL MATRIX (SITES AND ENVIRONMENTAL LAYERS) #####
sites_vars_all <- cbind(RS_matrix, PAM_onlyvars)
sites_vars_all <- as.data.frame(sites_vars_all) ##### MAKE IT A DATAFRAME OBJECT #####
colnames(sites_vars_all) <- c("Long", "Lat", "Range_size", "CCV", "Elevation", "CEH_min", "CEH_max", "CVH")
##### SAVE THE MATRIX #####
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(sites_vars_all, "all_data_siteslevel.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
####### GET THE VARIABLES FOR THE SITES LEVEL FOR EACH DATA SET (ALL, TONINI, LEACHE) #######
##### I NEED TO SUBSET THE PAMS #####
PAM_ton <- lets.subsetPAM(PAM_ahull, comp_ton$phy$tip.label, remove.cells = FALSE)
PAM_lea <- lets.subsetPAM(PAM_ahull, comp_lea$phy$tip.label, remove.cells = FALSE)
##### GET THE VALUES OF EACH VARIABLES AT SITES LEVEL #####
##### ALL #####
sites_vars_all <- sites_vars_all
##### TONINI #####
vars_PAM_ton <- lets.addvar(PAM_ton, stck_proj1, fun = median, onlyvar = TRUE)
RS_ton_sp <- lets.rangesize(PAM_ton, units = "squaremeter")
RS_ton_sp <- RS_ton_sp/1000000000000
RS_ton_sites <- lets.maplizer(PAM_ton, log(RS_ton_sp), rownames(RS_ton_sp), func = median, ras = TRUE) ##log
RS_ton_sites_mtx <- RS_ton_sites$Matrix
sites_vars_ton <- cbind(RS_ton_sites_mtx, vars_PAM_ton)
sites_vars_ton <- as.data.frame(sites_vars_ton)
colnames(sites_vars_ton) <- colnames(sites_vars_all)
##### LEACHE #####
vars_PAM_lea <- lets.addvar(PAM_lea, stck_proj1, fun = median, onlyvar = TRUE)
RS_lea_sp <- lets.rangesize(PAM_lea, units = "squaremeter")
RS_lea_sp <- RS_lea_sp/1000000000000
RS_lea_sites <- lets.maplizer(PAM_lea, log(RS_lea_sp), rownames(RS_lea_sp), func = median, ras = TRUE) ##log
RS_lea_sites_mtx <- RS_lea_sites$Matrix
sites_vars_lea <- cbind(RS_lea_sites_mtx, vars_PAM_lea)
sites_vars_lea <- as.data.frame(sites_vars_lea)
colnames(sites_vars_lea) <- colnames(sites_vars_all)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(sites_vars_ton, "ton_data_siteslevel.csv", row.names = FALSE)
write.csv(sites_vars_lea, "lea_data_siteslevel.csv", row.names = FALSE)
##### CHECK OF THERE ARE ANY DIFFERENCES IN THE NUMBER OF CELLS FOR TONINI AND LEACHE VS ALL THE SPECIES #####
cells_ton <- PAM_ahull$P[,comp_ton$phy$tip.label]
cells_ton_1 <- which(rowSums(cells_ton) == 0)
cells_lea <- PAM_ahull$P[,comp_lea$phy$tip.label]
cells_lea_1 <- which(rowSums(cells_lea) == 0)
##### THE SAME WITH THE SUBSET PAMS #####
cells_ton_2 <- which(rowSums(PAM_ton$P) == 0)
cells_lea_2 <- which(rowSums(PAM_lea$P) == 0)
save.image("G:/7_Doctorado/Cap_1/7_Codes/Rdata/3_PAM_and_layers.RData")
##### --------------------------------------------------------------------------------- #####
#############################################################################################
####################################### END OF THE CODE #####################################
#############################################################################################
