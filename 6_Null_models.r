
##### LOAD THE LIBRARIES #####
library(letsR)
library(raster)
library(rangemodelR)
library(spdep)
library(maptools)
library(ggplot2)
library(tibble)
library(caper)
library(broom)
#set.seed(08101994)
##### I NEED TO LOAD THEN PAM ######
load("E:/7_Doctorado/Cap_1/7_Codes/Rdata/5_PAM.RData")
##### GET RASTER #####
rast_full <- PAM_ahull$Richness_Raster * 0 ##### PAM_ahull es mi objeto tipo matriz de presencia ausencia contenido en el RData #####

##### CROP THE RASTER TO THE EXTENT OF ANALYSIS (ALL THE POLYGONS MERGED) #####
r_mask <- crop(rast_full, dist_ahull) ###### "dist_ahull" es un multipoligono que contiene el area de distribución de todas las 
                                      ##### especies, el cuál cargamos con el RData ######
r_mask <- mask(rast_full, dist_ahull)

##### CONVERTING FROM RASTER TO POLYGON ####
grids_shape <- rasterToPolygons(r_mask)

##### NOW I HAVE TO PLOT THE GRIDS #####
plot(grids_shape)

##### ASSING ID TO THE GRIDS #####
slot(grids_shape,"data") <- cbind("ID" = 1:length(grids_shape), slot(grids_shape, "data"))

##### CREATE A PAM IN THE GRID ######
projection(dist_ahull) <- projection(grids_shape)
PAM_grid <- lets.presab.grid(dist_ahull, grids_shape, "ID")

##### CREATE NEIGHBOUR LIST FOR THE RANGEMODELR SIMULATIONS #####
neighs <- poly2nb(grids_shape)

##### SAVE THE COORDINATES BASED CON RASTER VALUES != NA ##### 
grids_cells <- which(!is.na(values(r_mask)))
grids_ncell <- length(which(!is.na(values(r_mask))))
gridcoords <- xyFromCell(object = r_mask, cell = grids_cells)
colnames(gridcoords) <- c("Long", "Lat")

#### RUN THE FUNCTION "rangemod2d" ######
null_models <- rangemod2d(PAM_grid$PAM, grids_shape, "ID", nb = neighs, rsize = "observed",
                       var = NULL, reps = 1000, degen = TRUE)

##### ADDING COORDINATES TO THE MATRIX #####
grid_PAM <- cbind(gridcoords, PAM_grid$PAM)

##### PLOT RANDOM SPECIES RANGE #####
sim_ranges_1 <- which(null_models$degenerate.matrices[[3]][,64] == 1)
length(sim_ranges_1)

##### TRANSFORM THE ONE SPECIES IN RASTER ######
list_rast <- list()
for(i in 1:100){
sim_ranges_1 <- which(null_models$degenerate.matrices[[i]][,64] == 1)
null_model_raster_i <- rasterize(grid_PAM[sim_ranges_1, 1:2], r_mask, 1)

list_rast[i] <- null_model_raster_i

}
###### PLOT ######
plot(null_model_raster, col = c("white", "red"), legend = FALSE)
plot(grids_shape, add = T)

stck_rast <- stack(list_rast)
out_dir <- "G:/7_Doctorado/Cap_1/5_figures"

grid_sf <- st_as_sf(grids_shape)

anim <- tm_shape(grid_sf) + tm_polygons() + tm_shape(stck_rast) + tm_raster(palette = "red3", legend.show = FALSE) + tm_facets(nrow = 1, ncol = 1)

tmap_animation(anim, "anim_file.gif")

###### SAVE THE SIMULATED RANGES ######
sce_ranges <- vector("list", length(null_models$degenerate.matrices))
for (i in 1:length(null_models$degenerate.matrices)){
  ##### GET THE PAM #####
  sce_PAM <- null_models$degenerate.matrices[[i]]
  ###### GET EACH SPECIES SIMULATED RANGES #####
  sce_ranges[[i]] <- vector("list", length(ncol(sce_PAM)))
  for (j in 1:ncol(sce_PAM)){
    print(paste(i, '_', j, sep = ""))
    sim_ranges_i <- which(sce_PAM[, j] == 1)
    sce_PAM_df <- as.data.frame(sce_PAM)
    names_species <- colnames(sce_PAM_df[j])
    range_raster <- rasterize(grid_PAM[sim_ranges_i, 1:2], r_mask, 1)
    raster_to_polygon <- rasterToPolygons(range_raster) 
    raster_to_polygon_2 <- unionSpatialPolygons(raster_to_polygon, ID = rep(1, times = length(raster_to_polygon)), avoidGEOS = F)
    IDs <- sapply(slot(raster_to_polygon_2, "polygons"), function(x) slot(x, "ID"))
    data_frame <- data.frame(j)
    colnames(data_frame) <- c("IDs")
    raster_to_polygon_3 <- SpatialPolygonsDataFrame(raster_to_polygon_2, data_frame)
    raster_to_polygon_3$BINOMIAL <- names_species
    sce_ranges[[i]][[j]] <- raster_to_polygon_3}}

###### SAVE 1,000 SIMULATIONS #######
save.image("E:/7_Doctorado/Cap_1/7_Codes/Rdata/7_1000_sims.RData")

load("E:/7_Doctorado/Cap_1/7_Codes/Rdata/7_1000_sims.RData")
##### GET THE NAMES FOR THE SETS OF DATA (TONINI AND LEACHE) #######
leache <- read.nexus("C:/Users/Tsunnyk/Downloads/Filogenias/leache.nex")
tonini <- read.nexus("C:/Users/Tsunnyk/Downloads/Filogenias/tonini.nex")
ton_names <- as.data.frame(comp_ton$phy$tip.label)
colnames(ton_names) <- c("Species")
lea_names <- as.data.frame(comp_lea$phy$tip.label)
colnames(lea_names) <- c("Species")

###### CREATE AN EMPTY DATA FRAME #####
mydata <- list()
vars_list <- list()
vars_list_2 <- list()
mydata_sites_level <- list()

##### GET THE VARIABLES ######
for(i in 1:100){
  sim_i <- sce_ranges[[i]]
  sim_i_comb <- do.call(bind, sim_i) ##### COMBINE THE POLYGONS #####
  spatial_sim_i <- spTransform(sim_i_comb, CRS("+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m")) ##### PROJECT THE MULTIPOLYGON #####		 
  ###### RUN FOR EVERY SPECIES INSIDE THE SIMULATION ######
  sim_bbox <- bbox(spatial_sim_i)
  sim_long_min <- ahull_bbox[1]
  sim_long_max <- ahull_bbox[3]
  sim_lat_min <- ahull_bbox[2]
  sim_lat_max <- ahull_bbox[4]
  
  PAM_sim_i <- lets.presab(spatial_sim_i, count = TRUE, xmn = sim_long_min, xmx = sim_long_max, ymn = sim_lat_min, 
                           ymx = sim_lat_max, resol = 55500, crs = "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m", 
                           crs.grid = "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m", remove.cells = TRUE)
  
  ####### SUBSET THE PAMS #######
  PAM_sim_ton <- lets.subsetPAM(PAM_sim_i, comp_ton$phy$tip.label, remove.cells = TRUE)
  PAM_sim_lea <- lets.subsetPAM(PAM_sim_i, comp_lea$phy$tip.label, remove.cells = TRUE)
  
  ##### GET THE VALUES OF EACH VARIABLES AT SITES LEVEL #####
  ####### ALL #######
  ###### CREATE A MATRIX CONTAINING THE ENVIRONMENTAL LAYERS ###### 
  PAM_sim_env_i <- lets.addvar(PAM_sim_i, stck_proj1, fun = mean)
  vars <- PAM_sim_env_i[, -c(1:105)]
  coordinates_sites <- as.data.frame(PAM_sim_env_i[,c(1,2)])
  colnames(coordinates_sites) <- c("Long", "Lat")
  
  ###### CALCULATES NULL RANGES SIZES (SITES LEVEL) #####
  Null_RS_sites <- lets.maplizer(PAM_sim_i, log(RS_ahull), rownames(RS_ahull), func = median, ras = TRUE)
  Null_RS_sites_df <- as.data.frame(Null_RS_sites$Matrix)
  Null_all_sites_df <- cbind(coordinates_sites, Null_RS_sites_df[,3], vars)
  colnames(Null_all_sites_df) <- c("Long", "Lat", "Range_size", "CCV", "Elevation", "CEH_min", "CEH_max", "CVH")
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All")
  write.csv(Null_all_sites_df, paste0("RS_sites_all_sim_", i, ".csv"), row.names = FALSE)
  
  ####### TONINI #######
  ###### CREATE A MATRIX CONTAINING THE ENVIRONMENTAL LAYERS ###### 
  PAM_sim_env_ton <- lets.addvar(PAM_sim_ton, stck_proj1, fun = mean)
  vars_ton <- PAM_sim_env_ton[, -c(1:94)]
  coordinates_sites_ton <- as.data.frame(PAM_sim_env_ton[,c(1,2)])
  colnames(coordinates_sites_ton) <- c("Long", "Lat")
  
  ###### CALCULATES NULL RANGES SIZES (SITES LEVEL) #####
  Null_RS_sites_ton <- lets.maplizer(PAM_sim_ton, log(RS_ton_sp), rownames(RS_ton_sp), func = median, ras = TRUE)
  Null_RS_sites_df_ton <- as.data.frame(Null_RS_sites_ton$Matrix)
  Null_ton_sites_df <- cbind(coordinates_sites_ton, Null_RS_sites_df_ton[,3], vars_ton)
  colnames(Null_ton_sites_df) <- c("Long", "Lat", "Range_size", "CCV", "Elevation", "CEH_min", "CEH_max", "CVH")
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton")
  write.csv(Null_ton_sites_df, paste0("RS_sites_ton_sim_", i, ".csv"), row.names = FALSE)
  
  ####### LEACHE #######
  ###### CREATE A MATRIX CONTAINING THE ENVIRONMENTAL LAYERS ###### 
  PAM_sim_env_lea <- lets.addvar(PAM_sim_lea, stck_proj1, fun = mean)
  vars_lea <- PAM_sim_env_lea[, -c(1:83)]
  coordinates_sites_lea <- as.data.frame(PAM_sim_env_lea[,c(1,2)])
  colnames(coordinates_sites_lea) <- c("Long", "Lat")
  
  ###### CALCULATES NULL RANGES SIZES (SITES LEVEL) #####
  Null_RS_sites_lea <- lets.maplizer(PAM_sim_lea, log(RS_lea_sp), rownames(RS_lea_sp), func = median, ras = TRUE)
  Null_RS_sites_df_lea <- as.data.frame(Null_RS_sites_lea$Matrix)
  Null_lea_sites_df <- cbind(coordinates_sites_lea, Null_RS_sites_df_lea[,3], vars_lea)
  colnames(Null_lea_sites_df) <- c("Long", "Lat", "Range_size", "CCV", "Elevation", "CEH_min", "CEH_max", "CVH")
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea")
  write.csv(Null_lea_sites_df, paste0("RS_sites_lea_sim_", i, ".csv"), row.names = FALSE)
  
  
  ###### CALCULATES ENVIRONMENTAL MEANS (SPECIES LEVEL) ###### 
  species <- PAM_sim_env_i[, -c(1, 2, 106, 107, 108, 109, 110)]
  for (j in 1:103){
    
    PAM_vars_sub <- subset(vars, species[, j] == 1)
    
    vars_mean <- colMeans(PAM_vars_sub, na.rm = FALSE)
    
    vars_list[[j]] <- vars_mean
    
  }
  
  ####### GET DATA FRAME WITH THE ENVIRONMENTAL VALUES PER SPECIES #####
  vars_unlist <- data.frame(matrix(unlist(vars_list), nrow = length(vars_list), byrow = TRUE))
  colnames(vars_unlist) <- c("CCV", "Elevation", "CEH_min", "CEH_max", "CVH")
  
  ###### MIDPOINTS ######
  midpoint_i <- lets.midpoint(PAM_sim_i, planar = FALSE)
  midpoint_only <- as.data.frame(midpoint_i)
  colnames(midpoint_only) <- c("Species", "Long", "Lat")
  
  ####### RANGE SIZE #####
  RS_sim_sp <- RS_ahull1
  Null_all_sp_df <- cbind(midpoint_only, RS_ahull1[,2], vars_unlist)
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/All")
  write.csv(Null_all_sp_df, paste0("RS_sp_all_sim_", i, ".csv"), row.names = FALSE)
  
  Null_ton_sp_df <- merge(ton_names,  Null_all_sp_df, by = "Species", all = FALSE)
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton")
  write.csv(Null_ton_sp_df, paste0("RS_sp_ton_sim_", i, ".csv"), row.names = FALSE)
  
  Null_lea_sp_df <- merge(lea_names,  Null_all_sp_df, by = "Species", all = FALSE)
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea")
  write.csv(Null_lea_sp_df, paste0("RS_sp_lea_sim_", i, ".csv"), row.names = FALSE)
}

save.image("E:/7_Doctorado/Cap_1/7_Codes/Rdata/6_Null_models.RData")
####### OLS AND PGLS #####
load("E:/7_Doctorado/Cap_1/7_Codes/Rdata/6_Null_models.RData")

###############################################################################################
########################### PERFORM THE OLS AT SPECIES LEVEL ################################## 
###############################################################################################
###### LOAD THE PACKAGES #######
library(letsR)
library(raster)
library(rangemodelR)
library(spdep)
library(maptools)
library(ggplot2)
library(tibble)
library(MuMIn)
library(tidyr)

########### ALL ##############
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_all <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/All", 
                           pattern = ".csv")
list_ols_sp_all <- list()
####### OLS - SPECIES LEVEL #####
for (i in 1:length(null_list_sp_all)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/All")
  null_data_sp_all <- read.csv(null_list_sp_all[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_all$Long <- scale(null_data_sp_all$Long, center = T)
  null_data_sp_all$Lat <- scale(null_data_sp_all$Lat, center = T)
  null_data_sp_all$Range_size <- scale(null_data_sp_all$Range_size, center = T)
  null_data_sp_all$CCV <- scale(null_data_sp_all$CCV, center = T)
  null_data_sp_all$Elevation  <- scale(null_data_sp_all$Elevation, center = T)
  null_data_sp_all$CEH_min <- scale(null_data_sp_all$CEH_min, center = T)
  null_data_sp_all$CVH <- scale(null_data_sp_all$CVH, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sp_null_all <- lm(Range_size ~ CCV + Elevation + CVH + CEH_min, null_data_sp_all)
  ###### DREDGE ######
  ols_sp_null_all_dg <- dredge(ols_sp_null_all, rank = "AIC", m.lim = c(1,4))
  ###### MODELS WITH <= .95 #######
  ols_sp_null_all_sum <- summary(model.avg(ols_sp_null_all_dg, subset = cumsum(weight) <= .95))
  ###### DATAFRAME WITH THE COEFFICIENTS #######
  ols_sp_null_all_sum_df <- as.data.frame(ols_sp_null_all_sum$coefmat.full)
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sp_null_all_est <- ols_sp_null_all_sum_df[1]
  
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  ols_sp_null_all_bm <- get.models(ols_sp_null_all_dg, 1)[[1]]
  ##### GET THEIR COEFFICIENTS #######
  ols_sp_null_all_bm_s  <- as.data.frame(broom::glance(ols_sp_null_all_bm))
  ols_sp_null_all_bm_s_1 <- ols_sp_null_all_bm_s[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sp_null_all_est <- tibble::rownames_to_column(ols_sp_null_all_est, "Hypothesis")
  ols_sp_null_all_wider <- ols_sp_null_all_est %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sp_null_all_wider_df <- as.data.frame(ols_sp_null_all_wider[,-1])
  ols_sp_null_all_wider_df
  
  final_ols_df_sp_all <- cbind(Sim_number, ols_sp_null_all_wider_df, ols_sp_null_all_bm_s_1)
  
  list_ols_sp_all[[i]] <- final_ols_df_sp_all
}

final_ols_null_sp_all <- as.data.frame(do.call(rbind, list_ols_sp_all))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sp_all, "final_ols_null_sp_all.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

########### TONINI ##############
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_ton <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton", 
                               pattern = ".csv")
list_ols_sp_ton <- list()

for (i in 1:length(null_list_sp_ton)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton")
  null_data_sp_ton <- read.csv(null_list_sp_ton[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_ton$Long <- scale(null_data_sp_ton$Long, center = T)
  null_data_sp_ton$Lat <- scale(null_data_sp_ton$Lat, center = T)
  null_data_sp_ton$Range_size <- scale(null_data_sp_ton$Range_size, center = T)
  null_data_sp_ton$CCV <- scale(null_data_sp_ton$CCV, center = T)
  null_data_sp_ton$Elevation  <- scale(null_data_sp_ton$Elevation, center = T)
  null_data_sp_ton$CEH_min <- scale(null_data_sp_ton$CEH_min, center = T)
  null_data_sp_ton$CVH <- scale(null_data_sp_ton$CVH, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sp_null_ton <- lm(Range_size ~ CCV + Elevation + CVH + CEH_min, null_data_sp_ton)
  ###### DREDGE ######
  ols_sp_null_ton_dg <- dredge(ols_sp_null_ton, rank = "AIC", m.lim = c(1,4))
  ###### MODELS WITH <= .95 #######
  ols_sp_null_ton_sum <- summary(model.avg(ols_sp_null_ton_dg, subset = cumsum(weight) <= .95))
  ###### DATAFRAME WITH THE COEFFICIENTS #######
  ols_sp_null_ton_sum_df <- as.data.frame(ols_sp_null_ton_sum$coefmat.full)
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sp_null_ton_est <- ols_sp_null_ton_sum_df[1]
  
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  ols_sp_null_ton_bm <- get.models(ols_sp_null_ton_dg, 1)[[1]]
  ##### GET THEIR COEFFICIENTS #######
  ols_sp_null_ton_bm_s  <- as.data.frame(broom::glance(ols_sp_null_ton_bm))
  ols_sp_null_ton_bm_s_1 <- ols_sp_null_ton_bm_s[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sp_null_ton_est <- tibble::rownames_to_column(ols_sp_null_ton_est, "Hypothesis")
  ols_sp_null_ton_wider <- ols_sp_null_ton_est %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sp_null_ton_wider_df <- as.data.frame(ols_sp_null_ton_wider[,-1])
  ols_sp_null_ton_wider_df
  
  final_ols_df_sp_ton <- cbind(Sim_number, ols_sp_null_ton_wider_df, ols_sp_null_ton_bm_s_1)
  
  list_ols_sp_ton[[i]] <- final_ols_df_sp_ton
}


final_ols_null_sp_ton <- as.data.frame(do.call(rbind, list_ols_sp_ton))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sp_ton, "final_ols_null_sp_ton.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

########### LEACHE ##############
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_lea <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea", 
                               pattern = ".csv")
list_ols_sp_lea <- list()
####### OLS - SPECIES LEVEL #####
for (i in 1:length(null_list_sp_lea)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea")
  null_data_sp_lea <- read.csv(null_list_sp_lea[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_lea$Long <- scale(null_data_sp_lea$Long, center = T)
  null_data_sp_lea$Lat <- scale(null_data_sp_lea$Lat, center = T)
  null_data_sp_lea$Range_size <- scale(null_data_sp_lea$Range_size, center = T)
  null_data_sp_lea$CCV <- scale(null_data_sp_lea$CCV, center = T)
  null_data_sp_lea$Elevation  <- scale(null_data_sp_lea$Elevation, center = T)
  null_data_sp_lea$CEH_min <- scale(null_data_sp_lea$CEH_min, center = T)
  null_data_sp_lea$CVH <- scale(null_data_sp_lea$CVH, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sp_null_lea <- lm(Range_size ~ CCV + Elevation + CVH + CEH_min, null_data_sp_lea)
  ###### DREDGE ######
  ols_sp_null_lea_dg <- dredge(ols_sp_null_lea, rank = "AIC", m.lim = c(1,4))
  ###### MODELS WITH <= .95 #######
  ols_sp_null_lea_sum <- summary(model.avg(ols_sp_null_lea_dg, subset = cumsum(weight) <= .95))
  ###### DATAFRAME WITH THE COEFFICIENTS #######
  ols_sp_null_lea_sum_df <- as.data.frame(ols_sp_null_lea_sum$coefmat.full)
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sp_null_lea_est <- ols_sp_null_lea_sum_df[1]
  
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  ols_sp_null_lea_bm <- get.models(ols_sp_null_lea_dg, 1)[[1]]
  ##### GET THEIR COEFFICIENTS #######
  ols_sp_null_lea_bm_s  <- as.data.frame(broom::glance(ols_sp_null_lea_bm))
  ols_sp_null_lea_bm_s_1 <- ols_sp_null_lea_bm_s[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sp_null_lea_est <- tibble::rownames_to_column(ols_sp_null_lea_est, "Hypothesis")
  ols_sp_null_lea_wider <- ols_sp_null_lea_est %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sp_null_lea_wider_df <- as.data.frame(ols_sp_null_lea_wider[,-1])
  ols_sp_null_lea_wider_df
  
  final_ols_df_sp_lea <- cbind(Sim_number, ols_sp_null_lea_wider_df, ols_sp_null_lea_bm_s_1)
  
  list_ols_sp_lea[[i]] <- final_ols_df_sp_lea
}


final_ols_null_sp_lea <- as.data.frame(do.call(rbind, list_ols_sp_lea))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sp_lea, "final_ols_null_sp_lea.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

###############################################################################################
#################################### PERFORM THE PGLS #########################################
###############################################################################################
########### TONINI ##############
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_ton <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton", 
                               pattern = ".csv")
list_pgls_sp_ton <- list()
####### OLS #####
for (i in 1:length(null_list_sp_ton)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton")
  null_data_sp_ton <- read.csv(null_list_sp_ton[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_ton$Long <- scale(null_data_sp_ton$Long, center = T)
  null_data_sp_ton$Lat <- scale(null_data_sp_ton$Lat, center = T)
  null_data_sp_ton$Range_size <- scale(null_data_sp_ton$Range_size, center = T)
  null_data_sp_ton$CCV <- scale(null_data_sp_ton$CCV, center = T)
  null_data_sp_ton$Elevation  <- scale(null_data_sp_ton$Elevation, center = T)
  null_data_sp_ton$CEH_min <- scale(null_data_sp_ton$CEH_min, center = T)
  null_data_sp_ton$CVH <- scale(null_data_sp_ton$CVH, center = T)
  
  comp_null_ton <- comparative.data(tonini, data = null_data_sp_ton, names.col = "Species", vcv = T)
  ##### PERFORM THE PGLS ANALYSIS #####
  options(na.action = "na.fail")
  ###### PGLS #######
  pgls_sp_null_ton <- pgls(Range_size ~ CCV + Elevation + CEH_min + CVH, data = comp_null_ton, lambda = "ML")
  pgls_sp_null_ton_dg <- dredge(pgls_sp_null_ton, rank = "AIC", m.lim = c(1,4))
  pgls_sp_null_ton_dg
  pgls_sp_null_ton_sum <- summary(model.avg(pgls_sp_null_ton_dg, subset = cumsum(weight) <= .95))
  pgls_sp_null_ton_sum <- as.data.frame(pgls_sp_null_ton_sum$coefmat.full)
  pgls_sp_null_ton_est <- pgls_sp_null_ton_sum[1]
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  pgls_sp_null_ton_bm <- summary(get.models(pgls_sp_null_ton_dg, 1)[[1]])
  pgls_sp_null_ton_bm
  
  ##### GET THEIR COEFFICIENTS #######
  pgls_ton_null_r2 <- pgls_sp_null_ton_bm$adj.r.squared
  pgls_ton_null_p <- pgls_sp_null_ton_bm$coefficients[2,4]
  pgls_ton_null_r2_p <- cbind(pgls_ton_null_r2, pgls_ton_null_p)
  pgls_ton_null_r2_p_df <- as.data.frame(pgls_ton_null_r2_p)
  colnames(pgls_ton_null_r2_p_df) <- c("adj.r.squared ", "p.value")
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  pgls_sp_null_ton_est <- tibble::rownames_to_column(pgls_sp_null_ton_est, "Hypothesis")
  pgls_sp_null_ton_wider <- pgls_sp_null_ton_est %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  pgls_sp_null_ton_wider_df <- as.data.frame(pgls_sp_null_ton_wider[,-1])
  pgls_sp_null_ton_wider_df
  
  final_pgls_df_sp_ton <- cbind(Sim_number, pgls_sp_null_ton_wider_df, pgls_ton_null_r2_p_df)
  
  list_pgls_sp_ton[[i]] <- final_pgls_df_sp_ton
}


final_pgls_null_sp_ton <- as.data.frame(do.call(rbind, list_pgls_sp_ton))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_pgls_null_sp_ton, "final_pgls_null_sp_ton.csv", row.names = FALSE) 
##### --------------------------------------------------------------------------------- #####

########### LEACHE ##############
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_lea <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea", 
                               pattern = ".csv")
list_pgls_sp_lea <- list()

for (i in 1:length(null_list_sp_lea)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea")
  null_data_sp_lea <- read.csv(null_list_sp_lea[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_lea$Long <- scale(null_data_sp_lea$Long, center = T)
  null_data_sp_lea$Lat <- scale(null_data_sp_lea$Lat, center = T)
  null_data_sp_lea$Range_size <- scale(null_data_sp_lea$Range_size, center = T)
  null_data_sp_lea$CCV <- scale(null_data_sp_lea$CCV, center = T)
  null_data_sp_lea$Elevation  <- scale(null_data_sp_lea$Elevation, center = T)
  null_data_sp_lea$CEH_min <- scale(null_data_sp_lea$CEH_min, center = T)
  null_data_sp_lea$CVH <- scale(null_data_sp_lea$CVH, center = T)
  
  comp_null_lea <- comparative.data(tonini, data = null_data_sp_lea, names.col = "Species", vcv = T)
  ##### PERFORM THE PGLS ANALYSIS #####
  options(na.action = "na.fail")
  ###### PGLS #######
  pgls_sp_null_lea <- pgls(Range_size ~ CCV + Elevation + CEH_min + CVH, data = comp_null_lea, lambda = "ML")
  pgls_sp_null_lea_dg <- dredge(pgls_sp_null_lea, rank = "AIC", m.lim = c(1,4))
  pgls_sp_null_lea_dg
  pgls_sp_null_lea_sum <- summary(model.avg(pgls_sp_null_lea_dg, subset = cumsum(weight) <= .95))
  pgls_sp_null_lea_sum <- as.data.frame(pgls_sp_null_lea_sum$coefmat.full)
  pgls_sp_null_lea_est <- pgls_sp_null_lea_sum[1]
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  pgls_sp_null_lea_bm <- summary(get.models(pgls_sp_null_lea_dg, 1)[[1]])
  pgls_sp_null_lea_bm
  
  ##### GET THEIR COEFFICIENTS #######
  pgls_lea_null_r2 <- pgls_sp_null_lea_bm$adj.r.squared
  pgls_lea_null_p <- pgls_sp_null_lea_bm$coefficients[2,4]
  pgls_lea_null_r2_p <- cbind(pgls_lea_null_r2, pgls_lea_null_p)
  pgls_lea_null_r2_p_df <- as.data.frame(pgls_lea_null_r2_p)
  colnames(pgls_lea_null_r2_p_df) <- c("adj.r.squared ", "p.value")
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  pgls_sp_null_lea_est <- tibble::rownames_to_column(pgls_sp_null_lea_est, "Hypothesis")
  pgls_sp_null_lea_wider <- pgls_sp_null_lea_est %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  pgls_sp_null_lea_wider_df <- as.data.frame(pgls_sp_null_lea_wider[,-1])
  pgls_sp_null_lea_wider_df
  
  final_pgls_df_sp_lea <- cbind(Sim_number, pgls_sp_null_lea_wider_df, pgls_lea_null_r2_p_df)
  
  list_pgls_sp_lea[[i]] <- final_pgls_df_sp_lea
}

final_pgls_null_sp_lea <- as.data.frame(do.call(rbind, list_pgls_sp_lea))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_pgls_null_sp_lea, "final_pgls_null_sp_lea.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

###############################################################################################
########################### PERFORM THE OLS AT SITES LEVEL ####################################
###############################################################################################
null_list_sites_all <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All", 
                              pattern = ".csv")
list_ols_sites_all <- list()

for (i in 1:length(null_list_sites_all)){
####### READ THE SIMULATIONS ########
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All")
null_data_sites_all <- read.csv(null_list_sites_all[i])

####### SCALE THE VARIABLES #######
null_data_sites_all$Long <- scale(null_data_sites_all$Long, center = T)
null_data_sites_all$Lat <- scale(null_data_sites_all$Lat, center = T)
null_data_sites_all$Range_size <- scale(null_data_sites_all$Range_size, center = T)
null_data_sites_all$CCV <- scale(null_data_sites_all$CCV, center = T)
null_data_sites_all$Elevation  <- scale(null_data_sites_all$Elevation, center = T)
null_data_sites_all$CEH_min <- scale(null_data_sites_all$CEH_min, center = T)
null_data_sites_all$CVH <- scale(null_data_sites_all$CVH, center = T)

######## PERFORM THE OLS #######
options(na.action = "na.fail")
###### OLS #######
ols_sites_null_all <- lm(Range_size ~ CCV + Elevation + CVH + CEH_min + CCV, null_data_sites_all)
###### DREDGE ######
ols_sites_null_all_dg <- dredge(ols_sites_null_all, rank = "AIC", m.lim = c(1,4))
ols_sites_null_all_dg

##### CUMSUM FOR MODELS WITH MORE THAN ONE ACCUARED MODEL #####
if(length(which(cumsum(ols_sites_null_all_dg$weight) <= .95) == "TRUE") > 1){
###### MODELS WITH <= .95 #######
ols_sites_null_all_sum <- summary(model.avg(ols_sites_null_all_dg, subset = cumsum(weight) <= .95))
###### DATAFRAME WITH THE COEFFICIENTS #######
ols_sites_null_all_sum_df <- as.data.frame(ols_sites_null_all_sum$coefmat.full)
##### SELECT ONLY THE ESTIMATES ######
ols_sites_null_all_est <- ols_sites_null_all_sum_df[1]
}else if(length(which(cumsum(ols_sites_null_all_dg$weight) <= .95) == "TRUE") == 1 & nrow(summary(get.models(ols_sites_null_all_dg, 1)[[1]])$coefficients) <= 5){
  ols_sites_null_all_sum <- summary(model.avg(ols_sites_null_all_dg[c(1,2)]))  
  ols_sites_null_all_sum_df <- as.data.frame(ols_sites_null_all_sum$coefmat.full)
  ols_sites_null_all_est <- ols_sites_null_all_sum_df[1]
  } else{
  ###### DATAFRAME WITH THE ESTIMATES #######
  ols_sites_null_all_sum <- summary(get.models(ols_sites_null_all_dg, 1)[[1]])$coefficients
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sites_null_all_est <- as.data.frame(ols_sites_null_all_sum[,1])}
  
###### R2 AND P VALUE #####
###### GET THE MODEL WITH THE LOWEST AIC #######
ols_sites_null_all_bm <- get.models(ols_sites_null_all_dg, 1)[[1]]
##### GET THEIR COEFFICIENTS #######
ols_sites_null_all_bm_s  <- as.data.frame(broom::glance(ols_sites_null_all_bm))
ols_sites_null_all_bm_s_1 <- ols_sites_null_all_bm_s[,c(2,5,8,9)]  
#### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)

##### GET THE ESTIMATES IN TO COLUMNS #######
ols_sites_null_all_est <- tibble::rownames_to_column(ols_sites_null_all_est, "Hypothesis")
colnames(ols_sites_null_all_est) <- c("Hypothesis", "Estimate")
ols_sites_null_all_wider <- ols_sites_null_all_est %>%
pivot_wider(names_from = Hypothesis, values_from = Estimate)
ols_sites_null_all_wider_df <- as.data.frame(ols_sites_null_all_wider[,-1])
ols_sites_null_all_wider_df

final_ols_df_sites_all <- cbind(Sim_number, ols_sites_null_all_wider_df, ols_sites_null_all_bm_s_1)

list_ols_sites_all[[i]] <- final_ols_df_sites_all
}

final_ols_null_sites_all <- as.data.frame(do.call(rbind, list_ols_sites_all))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sites_all, "final_ols_null_sites_all.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

####### TONINI #####
null_list_sites_ton <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton", 
                                  pattern = ".csv")
list_ols_sites_ton <- list()

for (i in 1:length(null_list_sites_ton)){
  ####### READ THE SIMULATIONS ########
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton")
  null_data_sites_ton <- read.csv(null_list_sites_ton[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sites_ton$Long <- scale(null_data_sites_ton$Long, center = T)
  null_data_sites_ton$Lat <- scale(null_data_sites_ton$Lat, center = T)
  null_data_sites_ton$Range_size <- scale(null_data_sites_ton$Range_size, center = T)
  null_data_sites_ton$CCV <- scale(null_data_sites_ton$CCV, center = T)
  null_data_sites_ton$Elevation  <- scale(null_data_sites_ton$Elevation, center = T)
  null_data_sites_ton$CEH_min <- scale(null_data_sites_ton$CEH_min, center = T)
  null_data_sites_ton$CVH <- scale(null_data_sites_ton$CVH, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sites_null_ton <- lm(Range_size ~ CCV + Elevation + CVH + CEH_min + CCV, null_data_sites_ton)
  ###### DREDGE ######
  ols_sites_null_ton_dg <- dredge(ols_sites_null_ton, rank = "AIC", m.lim = c(1,4))
  ols_sites_null_ton_dg
  
  ##### CUMSUM FOR MODELS WITH MORE THAN ONE ACCUARED MODEL #####
  if(length(which(cumsum(ols_sites_null_ton_dg$weight) <= .95) == "TRUE") > 1){
    ###### MODELS WITH <= .95 #######
    ols_sites_null_ton_sum <- summary(model.avg(ols_sites_null_ton_dg, subset = cumsum(weight) <= .95))
    ###### DATAFRAME WITH THE COEFFICIENTS #######
    ols_sites_null_ton_sum_df <- as.data.frame(ols_sites_null_ton_sum$coefmat.full)
    ##### SELECT ONLY THE ESTIMATES ######
    ols_sites_null_ton_est <- ols_sites_null_ton_sum_df[1]
  }else if(length(which(cumsum(ols_sites_null_ton_dg$weight) <= .95) == "TRUE") == 1 & nrow(summary(get.models(ols_sites_null_ton_dg, 1)[[1]])$coefficients) <= 5){
    ols_sites_null_ton_sum <- summary(model.avg(ols_sites_null_ton_dg[c(1,2)]))  
    ols_sites_null_ton_sum_df <- as.data.frame(ols_sites_null_ton_sum$coefmat.full)
    ols_sites_null_ton_est <- ols_sites_null_ton_sum_df[1]
  } else{
    ###### DATAFRAME WITH THE ESTIMATES #######
    ols_sites_null_ton_sum <- summary(get.models(ols_sites_null_ton_dg, 1)[[1]])$coefficients
    ##### SELECT ONLY THE ESTIMATES ######
    ols_sites_null_ton_est <- as.data.frame(ols_sites_null_ton_sum[,1])}
  
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  ols_sites_null_ton_bm <- get.models(ols_sites_null_ton_dg, 1)[[1]]
  ##### GET THEIR COEFFICIENTS #######
  ols_sites_null_ton_bm_s  <- as.data.frame(broom::glance(ols_sites_null_ton_bm))
  ols_sites_null_ton_bm_s_1 <- ols_sites_null_ton_bm_s[,c(2,5,8,9)]  
  #### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sites_null_ton_est <- tibble::rownames_to_column(ols_sites_null_ton_est, "Hypothesis")
  colnames(ols_sites_null_ton_est) <- c("Hypothesis", "Estimate")
  ols_sites_null_ton_wider <- ols_sites_null_ton_est %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sites_null_ton_wider_df <- as.data.frame(ols_sites_null_ton_wider[,-1])
  ols_sites_null_ton_wider_df
  
  final_ols_df_sites_ton <- cbind(Sim_number, ols_sites_null_ton_wider_df, ols_sites_null_ton_bm_s_1)
  
  list_ols_sites_ton[[i]] <- final_ols_df_sites_ton
}

final_ols_null_sites_ton <- as.data.frame(do.call(rbind, list_ols_sites_ton))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sites_ton, "final_ols_null_sites_ton.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

####### LEACHE ##########
null_list_sites_lea <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea", 
                                  pattern = ".csv")
list_ols_sites_lea <- list()

for (i in 1:length(null_list_sites_lea)){
  ####### READ THE SIMULATIONS ########
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea")
  null_data_sites_lea <- read.csv(null_list_sites_lea[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sites_lea$Long <- scale(null_data_sites_lea$Long, center = T)
  null_data_sites_lea$Lat <- scale(null_data_sites_lea$Lat, center = T)
  null_data_sites_lea$Range_size <- scale(null_data_sites_lea$Range_size, center = T)
  null_data_sites_lea$CCV <- scale(null_data_sites_lea$CCV, center = T)
  null_data_sites_lea$Elevation  <- scale(null_data_sites_lea$Elevation, center = T)
  null_data_sites_lea$CEH_min <- scale(null_data_sites_lea$CEH_min, center = T)
  null_data_sites_lea$CVH <- scale(null_data_sites_lea$CVH, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sites_null_lea <- lm(Range_size ~ CCV + Elevation + CVH + CEH_min + CCV, null_data_sites_lea)
  ###### DREDGE ######
  ols_sites_null_lea_dg <- dredge(ols_sites_null_lea, rank = "AIC", m.lim = c(1,4))
  ols_sites_null_lea_dg
  
  ##### CUMSUM FOR MODELS WITH MORE THAN ONE ACCUARED MODEL #####
  if(length(which(cumsum(ols_sites_null_lea_dg$weight) <= .95) == "TRUE") > 1){
    ###### MODELS WITH <= .95 #######
    ols_sites_null_lea_sum <- summary(model.avg(ols_sites_null_lea_dg, subset = cumsum(weight) <= .95))
    ###### DATAFRAME WITH THE COEFFICIENTS #######
    ols_sites_null_lea_sum_df <- as.data.frame(ols_sites_null_lea_sum$coefmat.full)
    ##### SELECT ONLY THE ESTIMATES ######
    ols_sites_null_lea_est <- ols_sites_null_lea_sum_df[1]
  }else if(length(which(cumsum(ols_sites_null_lea_dg$weight) <= .95) == "TRUE") == 1 & nrow(summary(get.models(ols_sites_null_lea_dg, 1)[[1]])$coefficients) <= 5){
    ols_sites_null_lea_sum <- summary(model.avg(ols_sites_null_lea_dg[c(1,2)]))  
    ols_sites_null_lea_sum_df <- as.data.frame(ols_sites_null_lea_sum$coefmat.full)
    ols_sites_null_lea_est <- ols_sites_null_lea_sum_df[1]
  } else{
    ###### DATAFRAME WITH THE ESTIMATES #######
    ols_sites_null_lea_sum <- summary(get.models(ols_sites_null_lea_dg, 1)[[1]])$coefficients
    ##### SELECT ONLY THE ESTIMATES ######
    ols_sites_null_lea_est <- as.data.frame(ols_sites_null_lea_sum[,1])}
  
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  ols_sites_null_lea_bm <- get.models(ols_sites_null_lea_dg, 1)[[1]]
  ##### GET THEIR COEFFICIENTS #######
  ols_sites_null_lea_bm_s  <- as.data.frame(broom::glance(ols_sites_null_lea_bm))
  ols_sites_null_lea_bm_s_1 <- ols_sites_null_lea_bm_s[,c(2,5,8,9)]  
  #### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sites_null_lea_est <- tibble::rownames_to_column(ols_sites_null_lea_est, "Hypothesis")
  colnames(ols_sites_null_lea_est) <- c("Hypothesis", "Estimate")
  ols_sites_null_lea_wider <- ols_sites_null_lea_est %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sites_null_lea_wider_df <- as.data.frame(ols_sites_null_lea_wider[,-1])
  ols_sites_null_lea_wider_df
  
  final_ols_df_sites_lea <- cbind(Sim_number, ols_sites_null_lea_wider_df, ols_sites_null_lea_bm_s_1)
  
  list_ols_sites_lea[[i]] <- final_ols_df_sites_lea
}

final_ols_null_sites_lea <- as.data.frame(do.call(rbind, list_ols_sites_lea))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sites_lea, "final_ols_null_sites_lea.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

###############################################################################################
#################################### PERFORM THE SARS #########################################
###############################################################################################
library(ncf)
##### ALL SPECIES #####
null_list_sites_all <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All", 
                                  pattern = ".csv")
list_sar_sites_all <- list()

###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")

###### MATRIX OF SPATIAL COORDINATES #####
for(j in 1:length(null_list_sites_all)){
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All")
  
null_data_sites_all <- read.csv(null_list_sites_all[j])

null_data_sites_all$Long <- scale(null_data_sites_all$Long, center = T)
null_data_sites_all$Lat <- scale(null_data_sites_all$Lat, center = T)
null_data_sites_all$Range_size <- scale(null_data_sites_all$Range_size, center = T)
null_data_sites_all$CCV <- scale(null_data_sites_all$CCV, center = T)
null_data_sites_all$Elevation  <- scale(null_data_sites_all$Elevation, center = T)
null_data_sites_all$CEH_min <- scale(null_data_sites_all$CEH_min, center = T)
null_data_sites_all$CVH <- scale(null_data_sites_all$CVH, center = T)
  
sp_all <- SpatialPoints(data.frame(x = null_data_sites_all$Long, y = null_data_sites_all$Lat))
crs(sp_all) <- "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"
coords_all <- coordinates(sp_all)
###### MIN AND MAX NEAREST NEIGHBOR #####
k_all <- knn2nb(knearneigh(coords_all, k = 1))
dist_all <- unlist(nbdists(k_all, coords_all))
max_all <- max(dist_all)
###### CREATE NEIGHBOR MATRIX (DISTANCE IN METERS) #####
d_all <- dnearneigh(sp_all, longlat = FALSE, d1 = 0, d2 = max_all)
###### Create a data frame with path number and equations that will be tested in the pSEM #####
path_equation_all <- data.frame(path_number = paste0("path", 1))
path_equation_all$equation <-c("Range_size ~ CCV + Elevation + CEH_min + CVH")
##### Create models for each path and weighting scheme #####
dir_work_all <- "E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/All"
setwd(dir_work_all)
dir.create(paste0("sim_", j))
dir_work_all_i <- paste0("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/All/sim_", j)
setwd(dir_work_all_i)

for (i in 1:length(scheme)){
  ###### Create a directory by scheme #####
  setwd(dir_work_all_i)
  dir_scheme <- scheme[i]
  if(!dir.exists(dir_scheme)) dir.create(dir_scheme)
  setwd(dir_scheme)
  ##### Create weighting scheme distance matrix #####
  spatial_weights_d_all <- nb2listw(d_all, zero.policy = TRUE, style = scheme[i])
  for (p in 1:nrow(path_equation_all)) {
    schemes <- rep(c(paste(scheme[i])), times = 2)
    paths_all <- rep(c(paste(path_equation_all[p, 1])), times = 2)
    models_all <- c("lm_mod", "error")
    results_path_all <- data.frame(scheme = schemes, path = paths_all, model = models_all)
    ###### OLS model #####
    lm_mod_all <- lm(paste(path_equation_all[p, 2]), data = null_data_sites_all)
    lm_mod_s_all <- summary(lm_mod_all)
    R2_lm_all <- lm_mod_s_all[["adj.r.squared"]]
    coef_lm_all <- lm_mod_all$coefficients
    p_value_lm_all <- glance(lm_mod_all)$p.value
    ##### SAR error models #####
    ##### Minimum distance #####
    error_d_all <- spatialreg::errorsarlm(paste(path_equation_all[p, 2]), data = null_data_sites_all, 
                                          listw = spatial_weights_d_all, tol = 1e-12, zero.policy = T)
    error_d_s_all <-summary(error_d_all, Nagelkerke = TRUE)
    R2_error_d_all <- error_d_s_all$NK
    p_value_d_all <- error_d_s_all$Wald1$p.value
    ceof_d_all <- error_d_s_all$coefficients
    coef_all <- rbind(coef_lm_all, ceof_d_all)
    ##### Save R2 (pseudo R2 for errorsar) and AIC #####
    results_path_all$R2_models <- c(R2_lm_all, R2_error_d_all)
    results_path_all$AIC <- AIC(lm_mod_all, error_d_all)
    results_path_all$p_value <- c(p_value_lm_all, p_value_d_all)
    results_path_all <- cbind(results_path_all, coef_all)
    write.csv(results_path_all, file=paste(path_equation_all[p,1], ".csv"), row.names = F)
    ##### Make correlograms of residual autocorrelation #####
    cor.ols1.res_all <- correlog(as.numeric(null_data_sites_all$Long), as.numeric(null_data_sites_all$Lat), z = as.numeric(residuals(lm_mod_all)), na.rm = TRUE, increment = 1, resamp = 1)
    cor.sar1.res_all <- correlog(as.numeric(null_data_sites_all$Long), as.numeric(null_data_sites_all$Lat), z = as.numeric(residuals(error_d_all)), na.rm = TRUE, increment = 1, resamp = 1)
    ##### Save correlograms in a jpeg file #####
    jpeg(filename = paste(path_equation_all[p, 1], ".jpg"),
         width = 215, height = 279, units = "mm", res = 600)
    par(mfrow = c(2, 1))
    plot(cor.ols1.res_all, xlab = "Distance units", ylab = "Moran's I", ylim = c(-1,1), type = "l",
         lwd= 2, main = paste(models_all[1]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
    abline(h=0, lty=5)
    plot(cor.sar1.res_all, xlab = "Distance units", ylab = "Moran's I", ylim=c(-1,1), type = "l",
         lwd= 2, main =paste(models_all[2]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
    abline(h=0, lty=5)
    dev.off()
    ##### SAVE THE RESULTS #####
    save.image(file = paste(path_equation_all[p, 1], ".RData"))}}
}

###### READ THE SIMULATIONS #######
paths_sim <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/All",
                        pattern ="sim_",
                        full.names = TRUE)
sars_sim_list_all <- list()

##### SUMMARIZE THE RESULTS #####
for (i in 1:100){

  dir_work_all <- paths_sim[i]
  setwd(paste0(dir_work_all, "/S"))
  s_sata_all <- read.csv(paste0(dir_work_all, "/S/path1 .csv"))
  setwd(paste0(dir_work_all, "/C"))
  c_sata_all <- read.csv(paste0(dir_work_all, "/C/path1 .csv"))
  setwd(paste0(dir_work_all, "/W"))
  w_sata_all <- read.csv(paste0(dir_work_all, "/W/path1 .csv"))
  bind_data_all <- rbind(s_sata_all, c_sata_all, w_sata_all)
  bind_data_all$Simulation <- paste0("sim_", i)

results_all <- bind_data_all[which.min(bind_data_all$AIC.AIC),]
results_all1 <- results_all[, c(13, 9, 10, 11, 12, 4, 7, 6, 1, 2, 3, 5, 8)]

sars_sim_list_all[[i]] <- results_all1
}

sars_sim_list_all_1 <- as.data.frame(do.call(rbind, sars_sim_list_all))

##### Back to your working directory #####
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
##### Save the model selection results in a csv #####
write.csv(sars_sim_list_all_1, file="SAR_all.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####

##### TONINI #####
null_list_sites_ton <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton", 
                                  pattern = ".csv")
list_sar_sites_ton <- list()

###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")

###### MATRIX OF SPATIAL COORDINATES #####
for(j in 1:length(null_list_sites_ton)){
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton")
  null_data_sites_ton <- read.csv(null_list_sites_ton[j])
  
  null_data_sites_ton$Long <- scale(null_data_sites_ton$Long, center = T)
  null_data_sites_ton$Lat <- scale(null_data_sites_ton$Lat, center = T)
  null_data_sites_ton$Range_size <- scale(null_data_sites_ton$Range_size, center = T)
  null_data_sites_ton$CCV <- scale(null_data_sites_ton$CCV, center = T)
  null_data_sites_ton$Elevation  <- scale(null_data_sites_ton$Elevation, center = T)
  null_data_sites_ton$CEH_min <- scale(null_data_sites_ton$CEH_min, center = T)
  null_data_sites_ton$CVH <- scale(null_data_sites_ton$CVH, center = T)
  
  
  sp_ton <- SpatialPoints(data.frame(x = null_data_sites_ton$Long, y = null_data_sites_ton$Lat))
  crs(sp_ton) <- "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"
  coords_ton <- coordinates(sp_ton)
  ###### MIN AND MAX NEAREST NEIGHBOR #####
  k_ton <- knn2nb(knearneigh(coords_ton, k = 1))
  dist_ton <- unlist(nbdists(k_ton, coords_ton))
  max_ton <- max(dist_ton)
  ###### CREATE NEIGHBOR MATRIX (DISTANCE IN METERS) #####
  d_ton <- dnearneigh(sp_ton, longlat = FALSE, d1 = 0, d2 = max_ton)
  ###### Create a data frame with path number and equations that will be tested in the pSEM #####
  path_equation_ton <- data.frame(path_number = paste0("path", 1))
  path_equation_ton$equation <-c("Range_size ~ CCV + Elevation + CEH_min + CVH")
  ##### Create models for each path and weighting scheme #####
  dir_work_ton <- "E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Ton"
  setwd(dir_work_ton)
  dir.create(paste0("sim_", j))
  dir_work_ton_i <- paste0("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Ton/sim_", j)
  setwd(dir_work_ton_i)
  
  for (i in 1:length(scheme)){
    ###### Create a directory by scheme #####
    setwd(dir_work_ton_i)
    dir_scheme <- scheme[i]
    if(!dir.exists(dir_scheme)) dir.create(dir_scheme)
    setwd(dir_scheme)
    ##### Create weighting scheme distance matrix #####
    spatial_weights_d_ton <- nb2listw(d_ton, zero.policy = TRUE, style = scheme[i])
    for (p in 1:nrow(path_equation_ton)) {
      schemes <- rep(c(paste(scheme[i])), times = 2)
      paths_ton <- rep(c(paste(path_equation_ton[p, 1])), times = 2)
      models_ton <- c("lm_mod", "error")
      results_path_ton <- data.frame(scheme = schemes, path = paths_ton, model = models_ton)
      ###### OLS model #####
      lm_mod_ton <- lm(paste(path_equation_ton[p, 2]), data = null_data_sites_ton)
      lm_mod_s_ton <- summary(lm_mod_ton)
      R2_lm_ton <- lm_mod_s_ton[["adj.r.squared"]]
      coef_lm_ton <- lm_mod_ton$coefficients
      p_value_lm_ton <- glance(lm_mod_ton)$p.value
      ##### SAR error models #####
      ##### Minimum distance #####
      error_d_ton <- spatialreg::errorsarlm(paste(path_equation_ton[p, 2]), data = null_data_sites_ton, 
                                            listw = spatial_weights_d_ton, tol = 1e-12, zero.policy = T)
      error_d_s_ton <-summary(error_d_ton, Nagelkerke = TRUE)
      R2_error_d_ton <- error_d_s_ton$NK
      p_value_d_ton <- error_d_s_ton$Wald1$p.value
      ceof_d_ton <- error_d_s_ton$coefficients
      coef_ton <- rbind(coef_lm_ton, ceof_d_ton)
      ##### Save R2 (pseudo R2 for errorsar) and AIC #####
      results_path_ton$R2_models <- c(R2_lm_ton, R2_error_d_ton)
      results_path_ton$AIC <- AIC(lm_mod_ton, error_d_ton)
      results_path_ton$p_value <- c(p_value_lm_ton, p_value_d_ton)
      results_path_ton <- cbind(results_path_ton, coef_ton)
      write.csv(results_path_ton, file=paste(path_equation_ton[p,1], ".csv"), row.names = F)
      ##### Make correlograms of residual autocorrelation #####
      cor.ols1.res_ton <- correlog(as.numeric(null_data_sites_ton$Long), as.numeric(null_data_sites_ton$Lat), z = as.numeric(residuals(lm_mod_ton)), na.rm = TRUE, increment = 1, resamp = 1)
      cor.sar1.res_ton <- correlog(as.numeric(null_data_sites_ton$Long), as.numeric(null_data_sites_ton$Lat), z = as.numeric(residuals(error_d_ton)), na.rm = TRUE, increment = 1, resamp = 1)
      ##### Save correlograms in a jpeg file #####
      jpeg(filename = paste(path_equation_ton[p, 1], ".jpg"),
           width = 215, height = 279, units = "mm", res = 600)
      par(mfrow = c(2, 1))
      plot(cor.ols1.res_ton, xlab = "Distance units", ylab = "Moran's I", ylim = c(-1,1), type = "l",
           lwd= 2, main = paste(models_ton[1]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      plot(cor.sar1.res_ton, xlab = "Distance units", ylab = "Moran's I", ylim=c(-1,1), type = "l",
           lwd= 2, main =paste(models_ton[2]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      dev.off()
      ##### SAVE THE RESULTS #####
      save.image(file = paste(path_equation_ton[p, 1], ".RData"))}}
}

###### READ THE SIMULATIONS #######
paths_sim <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Ton",
                        pattern ="sim_",
                        full.names = TRUE)
sars_sim_list_ton <- list() 
##### SUMMARIZE THE RESULTS #####

for (i in 1:100){
  
  dir_work_ton <- paths_sim[i]
  setwd(paste0(dir_work_ton, "/S"))
  s_sata_ton <- read.csv(paste0(dir_work_ton, "/S/path1 .csv"))
  setwd(paste0(dir_work_ton, "/C"))
  c_sata_ton <- read.csv(paste0(dir_work_ton, "/C/path1 .csv"))
  setwd(paste0(dir_work_ton, "/W"))
  w_sata_ton <- read.csv(paste0(dir_work_ton, "/W/path1 .csv"))
  bind_data_ton <- rbind(s_sata_ton, c_sata_ton, w_sata_ton)
  bind_data_ton$Simulation <- paste0("sim_", i)
  
  results_ton <- bind_data_ton[which.min(bind_data_ton$AIC.AIC),]
  results_ton1 <- results_ton[, c(13, 9, 10, 11, 12, 4, 7, 6, 1, 2, 3, 5, 8)]
  
  sars_sim_list_ton[[i]] <- results_ton1
}

sars_sim_list_ton_1 <- as.data.frame(do.call(rbind, sars_sim_list_ton))

##### Back to your working directory #####
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
##### Save the model selection results in a csv #####
write.csv(sars_sim_list_ton_1, file="SAR_ton.csv", row.names = F)
##### --------------------------------------------------------------------------------- ####

##### LEACHE #####
null_list_sites_lea <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea", 
                                  pattern = ".csv")
list_sar_sites_lea <- list()

###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")

###### MATRIX OF SPATIAL COORDINATES #####
for(j in 1:length(null_list_sites_lea)){
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea")
  null_data_sites_lea <- read.csv(null_list_sites_lea[j])
  
  null_data_sites_lea$Long <- scale(null_data_sites_lea$Long, center = T)
  null_data_sites_lea$Lat <- scale(null_data_sites_lea$Lat, center = T) 
  null_data_sites_lea$Range_size <- scale(null_data_sites_lea$Range_size, center = T)
  null_data_sites_lea$CCV <- scale(null_data_sites_lea$CCV, center = T)
  null_data_sites_lea$Elevation  <- scale(null_data_sites_lea$Elevation, center = T)
  null_data_sites_lea$CEH_min <- scale(null_data_sites_lea$CEH_min, center = T)
  null_data_sites_lea$CVH <- scale(null_data_sites_lea$CVH, center = T)
  
  
  sp_lea <- SpatialPoints(data.frame(x = null_data_sites_lea$Long, y = null_data_sites_lea$Lat))
  crs(sp_lea) <- "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"
  coords_lea <- coordinates(sp_lea)
  ###### MIN AND MAX NEAREST NEIGHBOR #####
  k_lea <- knn2nb(knearneigh(coords_lea, k = 1))
  dist_lea <- unlist(nbdists(k_lea, coords_lea))
  max_lea <- max(dist_lea)
  ###### CREATE NEIGHBOR MATRIX (DISTANCE IN METERS) #####
  d_lea <- dnearneigh(sp_lea, longlat = FALSE, d1 = 0, d2 = max_lea)
  ###### Create a data frame with path number and equations that will be tested in the pSEM #####
  path_equation_lea <- data.frame(path_number = paste0("path", 1))
  path_equation_lea$equation <-c("Range_size ~ CCV + Elevation + CEH_min + CVH")
  ##### Create models for each path and weighting scheme #####
  dir_work_lea <- "E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Lea"
  setwd(dir_work_lea)
  dir.create(paste0("sim_", j))
  dir_work_lea_i <- paste0("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Lea/sim_", j)
  setwd(dir_work_lea_i)
  
  for (i in 1:length(scheme)){
    ###### Create a directory by scheme #####
    setwd(dir_work_lea_i)
    dir_scheme <- scheme[i]
    if(!dir.exists(dir_scheme)) dir.create(dir_scheme)
    setwd(dir_scheme)
    ##### Create weighting scheme distance matrix #####
    spatial_weights_d_lea <- nb2listw(d_lea, zero.policy = TRUE, style = scheme[i])
    for (p in 1:nrow(path_equation_lea)) {
      schemes <- rep(c(paste(scheme[i])), times = 2)
      paths_lea <- rep(c(paste(path_equation_lea[p, 1])), times = 2)
      models_lea <- c("lm_mod", "error")
      results_path_lea <- data.frame(scheme = schemes, path = paths_lea, model = models_lea)
      ###### OLS model #####
      lm_mod_lea <- lm(paste(path_equation_lea[p, 2]), data = null_data_sites_lea)
      lm_mod_s_lea <- summary(lm_mod_lea)
      R2_lm_lea <- lm_mod_s_lea[["adj.r.squared"]]
      coef_lm_lea <- lm_mod_lea$coefficients
      p_value_lm_lea <- glance(lm_mod_lea)$p.value
      ##### SAR error models #####
      ##### Minimum distance #####
      error_d_lea <- spatialreg::errorsarlm(paste(path_equation_lea[p, 2]), data = null_data_sites_lea, 
                                            listw = spatial_weights_d_lea, tol = 1e-12, zero.policy = T)
      error_d_s_lea <-summary(error_d_lea, Nagelkerke = TRUE)
      R2_error_d_lea <- error_d_s_lea$NK
      p_value_d_lea <- error_d_s_lea$Wald1$p.value
      ceof_d_lea <- error_d_s_lea$coefficients
      coef_lea <- rbind(coef_lm_lea, ceof_d_lea)
      ##### Save R2 (pseudo R2 for errorsar) and AIC #####
      results_path_lea$R2_models <- c(R2_lm_lea, R2_error_d_lea)
      results_path_lea$AIC <- AIC(lm_mod_lea, error_d_lea)
      results_path_lea$p_value <- c(p_value_lm_lea, p_value_d_lea)
      results_path_lea <- cbind(results_path_lea, coef_lea)
      write.csv(results_path_lea, file=paste(path_equation_lea[p,1], ".csv"), row.names = F)
      ##### Make correlograms of residual autocorrelation #####
      cor.ols1.res_lea <- correlog(as.numeric(null_data_sites_lea$Long), as.numeric(null_data_sites_lea$Lat), z = as.numeric(residuals(lm_mod_lea)), na.rm = TRUE, increment = 1, resamp = 1)
      cor.sar1.res_lea <- correlog(as.numeric(null_data_sites_lea$Long), as.numeric(null_data_sites_lea$Lat), z = as.numeric(residuals(error_d_lea)), na.rm = TRUE, increment = 1, resamp = 1)
      ##### Save correlograms in a jpeg file #####
      jpeg(filename = paste(path_equation_lea[p, 1], ".jpg"),
           width = 215, height = 279, units = "mm", res = 600)
      par(mfrow = c(2, 1))
      plot(cor.ols1.res_lea, xlab = "Distance units", ylab = "Moran's I", ylim = c(-1,1), type = "l",
           lwd= 2, main = paste(models_lea[1]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      plot(cor.sar1.res_lea, xlab = "Distance units", ylab = "Moran's I", ylim=c(-1,1), type = "l",
           lwd= 2, main =paste(models_lea[2]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      dev.off()
      ##### SAVE THE RESULTS #####
      save.image(file = paste(path_equation_lea[p, 1], ".RData"))}}
}

###### READ THE SIMULATIONS #######
paths_sim <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Lea",
                        pattern ="sim_",
                        full.names = TRUE)
sars_sim_list_lea <- list() 
##### SUMMARIZE THE RESULTS #####

for (i in 1:100){
  
  dir_work_lea <- paths_sim[i]
  setwd(paste0(dir_work_lea, "/S"))
  s_sata_lea <- read.csv(paste0(dir_work_lea, "/S/path1 .csv"))
  setwd(paste0(dir_work_lea, "/C"))
  c_sata_lea <- read.csv(paste0(dir_work_lea, "/C/path1 .csv"))
  setwd(paste0(dir_work_lea, "/W"))
  w_sata_lea <- read.csv(paste0(dir_work_lea, "/W/path1 .csv"))
  bind_data_lea <- rbind(s_sata_lea, c_sata_lea, w_sata_lea)
  bind_data_lea$Simulation <- paste0("sim_", i)
  
  results_lea <- bind_data_lea[which.min(bind_data_lea$AIC.AIC),]
  results_lea1 <- results_lea[, c(13, 9, 10, 11, 12, 4, 7, 6, 1, 2, 3, 5, 8)]
  
  sars_sim_list_lea[[i]] <- results_lea1
}

sars_sim_list_lea_1 <- as.data.frame(do.call(rbind, sars_sim_list_lea))

##### Back to your working directory #####
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
##### Save the model selection results in a csv #####
write.csv(sars_sim_list_lea_1, file="SAR_lea.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####

###############################################################################################
########################## PERFORM THE OLS FOR RAPOPORTS RULE SPECIES LEVEL ###################
###############################################################################################
######## ALL ##########
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_all_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/All", 
                               pattern = ".csv")
list_ols_sp_all_rapo <- list()
for (i in 1:length(null_list_sp_all_rapo)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/All")
  null_data_sp_all_rapo <- read.csv(null_list_sp_all_rapo[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_all_rapo$Range_size <- scale(null_data_sp_all_rapo$Range_size, center = T)
  null_data_sp_all_rapo$Lat <- scale(null_data_sp_all_rapo$Lat, center = T)

  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sp_null_all_rapo <- lm(Range_size ~ Lat, null_data_sp_all_rapo)
  ###### DATAFRAME WITH THE ESTIMATES #######
  ols_sp_null_all_sum_rapo <- summary(ols_sp_null_all_rapo)$coefficients
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sp_null_all_est_rapo <- as.data.frame(ols_sp_null_all_sum_rapo[,1])

  
  ###### R2 AND P VALUE #####
  ##### GET THEIR COEFFICIENTS #######
  ols_sp_null_all_bm_s_rapo  <- as.data.frame(broom::glance(ols_sp_null_all_rapo))
  ols_sp_null_all_bm_s_1_rapo <- ols_sp_null_all_bm_s_rapo[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sp_null_all_est_rapo <- tibble::rownames_to_column(ols_sp_null_all_est_rapo, "Hypothesis")
  colnames(ols_sp_null_all_est_rapo) <- c("Hypothesis", "Estimate")
  
  
  ols_sp_null_all_wider_rapo <- ols_sp_null_all_est_rapo %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sp_null_all_wider_df_rapo <- as.data.frame(ols_sp_null_all_wider_rapo[,-1])
  ols_sp_null_all_wider_df_rapo
  
  final_ols_df_sp_all_rapo <- cbind(Sim_number, ols_sp_null_all_wider_df_rapo, ols_sp_null_all_bm_s_1_rapo)
  
  list_ols_sp_all_rapo[[i]] <- final_ols_df_sp_all_rapo
}


final_ols_null_sp_all_rapo <- as.data.frame(do.call(rbind, list_ols_sp_all_rapo))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sp_all_rapo, "final_ols_null_sp_all_rapo.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

######## TONINI ##########
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_ton_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton", 
                                    pattern = ".csv")
list_ols_sp_ton_rapo <- list()
####### OLS - SPECIES LEVEL #####
for (i in 1:length(null_list_sp_ton_rapo)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton")
  null_data_sp_ton_rapo <- read.csv(null_list_sp_ton_rapo[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_ton_rapo$Range_size <- scale(null_data_sp_ton_rapo$Range_size, center = T)
  null_data_sp_ton_rapo$Lat <- scale(null_data_sp_ton_rapo$Lat, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sp_null_ton_rapo <- lm(Range_size ~ Lat, null_data_sp_ton_rapo)
  ###### DATAFRAME WITH THE ESTIMATES #######
  ols_sp_null_ton_sum_rapo <- summary(ols_sp_null_ton_rapo)$coefficients
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sp_null_ton_est_rapo <- as.data.frame(ols_sp_null_ton_sum_rapo[,1])
  
  
  ###### R2 AND P VALUE #####
  ##### GET THEIR COEFFICIENTS #######
  ols_sp_null_ton_bm_s_rapo  <- as.data.frame(broom::glance(ols_sp_null_ton_rapo))
  ols_sp_null_ton_bm_s_1_rapo <- ols_sp_null_ton_bm_s_rapo[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sp_null_ton_est_rapo <- tibble::rownames_to_column(ols_sp_null_ton_est_rapo, "Hypothesis")
  colnames(ols_sp_null_ton_est_rapo) <- c("Hypothesis", "Estimate")
  
  
  ols_sp_null_ton_wider_rapo <- ols_sp_null_ton_est_rapo %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sp_null_ton_wider_df_rapo <- as.data.frame(ols_sp_null_ton_wider_rapo[,-1])
  ols_sp_null_ton_wider_df_rapo
  
  final_ols_df_sp_ton_rapo <- cbind(Sim_number, ols_sp_null_ton_wider_df_rapo, ols_sp_null_ton_bm_s_1_rapo)
  
  list_ols_sp_ton_rapo[[i]] <- final_ols_df_sp_ton_rapo
}


final_ols_null_sp_ton_rapo <- as.data.frame(do.call(rbind, list_ols_sp_ton_rapo))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sp_ton_rapo, "final_ols_null_sp_ton_rapo.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

######## LEACHE ##############
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_lea_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea", 
                                    pattern = ".csv")
list_ols_sp_lea_rapo <- list()
####### OLS - SPECIES LEVEL #####
for (i in 1:length(null_list_sp_lea_rapo)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea")
  null_data_sp_lea_rapo <- read.csv(null_list_sp_lea_rapo[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_lea_rapo$Range_size <- scale(null_data_sp_lea_rapo$Range_size, center = T)
  null_data_sp_lea_rapo$Lat <- scale(null_data_sp_lea_rapo$Lat, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sp_null_lea_rapo <- lm(Range_size ~ Lat, null_data_sp_lea_rapo)
  ###### DATAFRAME WITH THE ESTIMATES #######
  ols_sp_null_lea_sum_rapo <- summary(ols_sp_null_lea_rapo)$coefficients
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sp_null_lea_est_rapo <- as.data.frame(ols_sp_null_lea_sum_rapo[,1])
  
  
  ###### R2 AND P VALUE #####
  ##### GET THEIR COEFFICIENTS #######
  ols_sp_null_lea_bm_s_rapo  <- as.data.frame(broom::glance(ols_sp_null_lea_rapo))
  ols_sp_null_lea_bm_s_1_rapo <- ols_sp_null_lea_bm_s_rapo[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sp_null_lea_est_rapo <- tibble::rownames_to_column(ols_sp_null_lea_est_rapo, "Hypothesis")
  colnames(ols_sp_null_lea_est_rapo) <- c("Hypothesis", "Estimate")
  
  
  ols_sp_null_lea_wider_rapo <- ols_sp_null_lea_est_rapo %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sp_null_lea_wider_df_rapo <- as.data.frame(ols_sp_null_lea_wider_rapo[,-1])
  ols_sp_null_lea_wider_df_rapo
  
  final_ols_df_sp_lea_rapo <- cbind(Sim_number, ols_sp_null_lea_wider_df_rapo, ols_sp_null_lea_bm_s_1_rapo)
  
  list_ols_sp_lea_rapo[[i]] <- final_ols_df_sp_lea_rapo
}


final_ols_null_sp_lea_rapo <- as.data.frame(do.call(rbind, list_ols_sp_lea_rapo))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sp_lea_rapo, "final_ols_null_sp_lea_rapo.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
###############################################################################################
############################### PERFORM THE PGLS FOR RAPOPORTS RULE ###########################
###############################################################################################
######## TONINI ###########
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_ton_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton", 
                                    pattern = ".csv")
list_pgls_sp_ton_rapo <- list()
####### OLS #####
for (i in 1:length(null_list_sp_ton_rapo)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Ton")
  null_data_sp_ton_rapo <- read.csv(null_list_sp_ton_rapo[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_ton_rapo$Range_size <- scale(null_data_sp_ton_rapo$Range_size, center = T)
  null_data_sp_ton_rapo$Lat <- scale(null_data_sp_ton_rapo$Lat, center = T)
  
  
  comp_null_ton_rapo <- comparative.data(tonini, data = null_data_sp_ton_rapo, names.col = "Species", vcv = T)
  ##### PERFORM THE PGLS ANALYSIS #####
  options(na.action = "na.fail")
  ###### PGLS #######
  pgls_sp_null_ton_rapo <- pgls(Range_size ~ Lat, data = comp_null_ton_rapo, lambda = "ML")
  pgls_sp_null_ton_rapo_sum <- summary(pgls_sp_null_ton_rapo)
  pgls_sp_null_ton_rapo_sum <- as.data.frame(pgls_sp_null_ton_rapo_sum$coefficients)
  pgls_sp_null_ton_rapo_est <- pgls_sp_null_ton_rapo_sum[1]
  
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  pgls_sp_null_ton_rapo_bm <- summary(pgls_sp_null_ton_rapo)
  pgls_sp_null_ton_rapo_bm
  
  ##### GET THEIR COEFFICIENTS #######
  pgls_ton_rapo_null_r2 <- pgls_sp_null_ton_rapo_bm$adj.r.squared
  pgls_ton_rapo_null_p <- pgls_sp_null_ton_rapo_bm$coefficients[2,4]
  pgls_ton_rapo_null_r2_p <- cbind(pgls_ton_rapo_null_r2, pgls_ton_rapo_null_p)
  pgls_ton_rapo_null_r2_p_df <- as.data.frame(pgls_ton_rapo_null_r2_p)
  colnames(pgls_ton_rapo_null_r2_p_df) <- c("adj.r.squared ", "p.value")
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  pgls_sp_null_ton_rapo_est <- tibble::rownames_to_column(pgls_sp_null_ton_rapo_est, "Hypothesis")
  pgls_sp_null_ton_rapo_wider <- pgls_sp_null_ton_rapo_est %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  pgls_sp_null_ton_rapo_wider_df <- as.data.frame(pgls_sp_null_ton_rapo_wider[,-1])
  pgls_sp_null_ton_rapo_wider_df
  
  final_pgls_df_sp_ton_rapo <- cbind(Sim_number, pgls_sp_null_ton_rapo_wider_df, pgls_ton_rapo_null_r2_p_df)
  
  list_pgls_sp_ton_rapo[[i]] <- final_pgls_df_sp_ton_rapo
}

final_pgls_null_sp_ton_rapo <- as.data.frame(do.call(rbind, list_pgls_sp_ton_rapo))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_pgls_null_sp_ton_rapo, "final_pgls_null_sp_ton_rapo.csv", row.names = FALSE) 
##### --------------------------------------------------------------------------------- #####
######## LEACHE ###########
###### CREATE A LIST CONTAINING ALL THE SIMULATIONS ########
null_list_sp_lea_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea", 
                                    pattern = ".csv")
list_pgls_sp_lea_rapo <- list()
####### OLS #####
for (i in 1:length(null_list_sp_lea_rapo)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/species_level/Lea")
  null_data_sp_lea_rapo <- read.csv(null_list_sp_lea_rapo[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sp_lea_rapo$Range_size <- scale(null_data_sp_lea_rapo$Range_size, center = T)
  null_data_sp_lea_rapo$Lat <- scale(null_data_sp_lea_rapo$Lat, center = T)
  
  
  comp_null_lea_rapo <- comparative.data(tonini, data = null_data_sp_lea_rapo, names.col = "Species", vcv = T)
  ##### PERFORM THE PGLS ANALYSIS #####
  options(na.action = "na.fail")
  ###### PGLS #######
  pgls_sp_null_lea_rapo <- pgls(Range_size ~ Lat, data = comp_null_lea_rapo, lambda = "ML")
  pgls_sp_null_lea_rapo_sum <- summary(pgls_sp_null_lea_rapo)
  pgls_sp_null_lea_rapo_sum <- as.data.frame(pgls_sp_null_lea_rapo_sum$coefficients)
  pgls_sp_null_lea_rapo_est <- pgls_sp_null_lea_rapo_sum[1]
  
  ###### R2 AND P VALUE #####
  ###### GET THE MODEL WITH THE LOWEST AIC #######
  pgls_sp_null_lea_rapo_bm <- summary(pgls_sp_null_lea_rapo)
  pgls_sp_null_lea_rapo_bm
  
  ##### GET THEIR COEFFICIENTS #######
  pgls_lea_rapo_null_r2 <- pgls_sp_null_lea_rapo_bm$adj.r.squared
  pgls_lea_rapo_null_p <- pgls_sp_null_lea_rapo_bm$coefficients[2,4]
  pgls_lea_rapo_null_r2_p <- cbind(pgls_lea_rapo_null_r2, pgls_lea_rapo_null_p)
  pgls_lea_rapo_null_r2_p_df <- as.data.frame(pgls_lea_rapo_null_r2_p)
  colnames(pgls_lea_rapo_null_r2_p_df) <- c("adj.r.squared ", "p.value")
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  pgls_sp_null_lea_rapo_est <- tibble::rownames_to_column(pgls_sp_null_lea_rapo_est, "Hypothesis")
  pgls_sp_null_lea_rapo_wider <- pgls_sp_null_lea_rapo_est %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  pgls_sp_null_lea_rapo_wider_df <- as.data.frame(pgls_sp_null_lea_rapo_wider[,-1])
  pgls_sp_null_lea_rapo_wider_df
  
  final_pgls_df_sp_lea_rapo <- cbind(Sim_number, pgls_sp_null_lea_rapo_wider_df, pgls_lea_rapo_null_r2_p_df)
  
  list_pgls_sp_lea_rapo[[i]] <- final_pgls_df_sp_lea_rapo
}

final_pgls_null_sp_lea_rapo <- as.data.frame(do.call(rbind, list_pgls_sp_lea_rapo))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_pgls_null_sp_lea_rapo, "final_pgls_null_sp_lea_rapo.csv", row.names = FALSE) 
##### --------------------------------------------------------------------------------- #####

###############################################################################################
########################## PERFORM THE OLS FOR RAPOPORTS RULE SITES LEVEL #####################
###############################################################################################
############ ALL ###############
null_list_sites_all_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All", 
                                       pattern = ".csv")
list_ols_sites_all_rapo <- list()
for (i in 1:length(null_list_sites_all_rapo)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All")
  null_data_sites_all_rapo <- read.csv(null_list_sites_all_rapo[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sites_all_rapo$Range_size <- scale(null_data_sites_all_rapo$Range_size, center = T)
  null_data_sites_all_rapo$Lat <- scale(null_data_sites_all_rapo$Lat, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sites_null_all_rapo <- lm(Range_size ~ Lat, null_data_sites_all_rapo)
  ###### DATAFRAME WITH THE ESTIMATES #######
  ols_sites_null_all_sum_rapo <- summary(ols_sites_null_all_rapo)$coefficients
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sites_null_all_est_rapo <- as.data.frame(ols_sites_null_all_sum_rapo[,1])
  
  
  ###### R2 AND P VALUE #####
  ##### GET THEIR COEFFICIENTS #######
  ols_sites_null_all_bm_s_rapo  <- as.data.frame(broom::glance(ols_sites_null_all_rapo))
  ols_sites_null_all_bm_s_1_rapo <- ols_sites_null_all_bm_s_rapo[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sites_null_all_est_rapo <- tibble::rownames_to_column(ols_sites_null_all_est_rapo, "Hypothesis")
  colnames(ols_sites_null_all_est_rapo) <- c("Hypothesis", "Estimate")
  
  
  ols_sites_null_all_wider_rapo <- ols_sites_null_all_est_rapo %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sites_null_all_wider_df_rapo <- as.data.frame(ols_sites_null_all_wider_rapo[,-1])
  ols_sites_null_all_wider_df_rapo
  
  final_ols_df_sites_all_rapo <- cbind(Sim_number, ols_sites_null_all_wider_df_rapo, ols_sites_null_all_bm_s_1_rapo)
  
  list_ols_sites_all_rapo[[i]] <- final_ols_df_sites_all_rapo
}


final_ols_null_sites_all_rapo <- as.data.frame(do.call(rbind, list_ols_sites_all_rapo))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sites_all_rapo, "final_ols_null_sites_all_rapo.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

############ TONINI ###############
null_list_sites_ton_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton", 
                                       pattern = ".csv")
list_ols_sites_ton_rapo <- list()
for (i in 1:length(null_list_sites_ton_rapo)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton")
  null_data_sites_ton_rapo <- read.csv(null_list_sites_ton_rapo[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sites_ton_rapo$Range_size <- scale(null_data_sites_ton_rapo$Range_size, center = T)
  null_data_sites_ton_rapo$Lat <- scale(null_data_sites_ton_rapo$Lat, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sites_null_ton_rapo <- lm(Range_size ~ Lat, null_data_sites_ton_rapo)
  ###### DATAFRAME WITH THE ESTIMATES #######
  ols_sites_null_ton_sum_rapo <- summary(ols_sites_null_ton_rapo)$coefficients
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sites_null_ton_est_rapo <- as.data.frame(ols_sites_null_ton_sum_rapo[,1])
  
  
  ###### R2 AND P VALUE #####
  ##### GET THEIR COEFFICIENTS #######
  ols_sites_null_ton_bm_s_rapo  <- as.data.frame(broom::glance(ols_sites_null_ton_rapo))
  ols_sites_null_ton_bm_s_1_rapo <- ols_sites_null_ton_bm_s_rapo[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sites_null_ton_est_rapo <- tibble::rownames_to_column(ols_sites_null_ton_est_rapo, "Hypothesis")
  colnames(ols_sites_null_ton_est_rapo) <- c("Hypothesis", "Estimate")
  
  
  ols_sites_null_ton_wider_rapo <- ols_sites_null_ton_est_rapo %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sites_null_ton_wider_df_rapo <- as.data.frame(ols_sites_null_ton_wider_rapo[,-1])
  ols_sites_null_ton_wider_df_rapo
  
  final_ols_df_sites_ton_rapo <- cbind(Sim_number, ols_sites_null_ton_wider_df_rapo, ols_sites_null_ton_bm_s_1_rapo)
  
  list_ols_sites_ton_rapo[[i]] <- final_ols_df_sites_ton_rapo
}


final_ols_null_sites_ton_rapo <- as.data.frame(do.call(rbind, list_ols_sites_ton_rapo))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sites_ton_rapo, "final_ols_null_sites_ton_rapo.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

############ LEACHE ###############
null_list_sites_lea_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea", 
                                       pattern = ".csv")
list_ols_sites_lea_rapo <- list()
for (i in 1:length(null_list_sites_lea_rapo)){
  ####### READ THE SIMULATIONS 
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea")
  null_data_sites_lea_rapo <- read.csv(null_list_sites_lea_rapo[i])
  
  ####### SCALE THE VARIABLES #######
  null_data_sites_lea_rapo$Range_size <- scale(null_data_sites_lea_rapo$Range_size, center = T)
  null_data_sites_lea_rapo$Lat <- scale(null_data_sites_lea_rapo$Lat, center = T)
  
  ######## PERFORM THE OLS #######
  options(na.action = "na.fail")
  ###### OLS #######
  ols_sites_null_lea_rapo <- lm(Range_size ~ Lat, null_data_sites_lea_rapo)
  ###### DATAFRAME WITH THE ESTIMATES #######
  ols_sites_null_lea_sum_rapo <- summary(ols_sites_null_lea_rapo)$coefficients
  ##### SELECT ONLY THE ESTIMATES ######
  ols_sites_null_lea_est_rapo <- as.data.frame(ols_sites_null_lea_sum_rapo[,1])
  
  
  ###### R2 AND P VALUE #####
  ##### GET THEIR COEFFICIENTS #######
  ols_sites_null_lea_bm_s_rapo  <- as.data.frame(broom::glance(ols_sites_null_lea_rapo))
  ols_sites_null_lea_bm_s_1_rapo <- ols_sites_null_lea_bm_s_rapo[,c(2,5,8,9)]  
  ##### SET THE SIMULATION NUMBER ########
  Sim_number <- paste0("Sim_", i)
  
  ##### GET THE ESTIMATES IN TO COLUMNS #######
  ols_sites_null_lea_est_rapo <- tibble::rownames_to_column(ols_sites_null_lea_est_rapo, "Hypothesis")
  colnames(ols_sites_null_lea_est_rapo) <- c("Hypothesis", "Estimate")
  
  
  ols_sites_null_lea_wider_rapo <- ols_sites_null_lea_est_rapo %>%
    pivot_wider(names_from = Hypothesis, values_from = Estimate)
  ols_sites_null_lea_wider_df_rapo <- as.data.frame(ols_sites_null_lea_wider_rapo[,-1])
  ols_sites_null_lea_wider_df_rapo
  
  final_ols_df_sites_lea_rapo <- cbind(Sim_number, ols_sites_null_lea_wider_df_rapo, ols_sites_null_lea_bm_s_1_rapo)
  
  list_ols_sites_lea_rapo[[i]] <- final_ols_df_sites_lea_rapo
}


final_ols_null_sites_lea_rapo <- as.data.frame(do.call(rbind, list_ols_sites_lea_rapo))
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
write.csv(final_ols_null_sites_lea_rapo, "final_ols_null_sites_lea_rapo.csv", row.names = FALSE)

###############################################################################################
################################# PERFORM THE SARS FOR RAPOPORTS RULE #########################
###############################################################################################
############### ALL SPECIES ###############

library(ncf)
##### ALL SPECIES #####
null_list_sites_all_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All", 
                                       pattern = ".csv")
list_sar_sites_all_rapo <- list()

###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")

###### MATRIX OF SPATIAL COORDINATES #####
for(j in 1:length(null_list_sites_all_rapo)){
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/All")
  null_data_sites_all_rapo <- read.csv(null_list_sites_all_rapo[j])
  
  null_data_sites_all_rapo$Range_size <- scale(null_data_sites_all_rapo$Range_size, center = T)
  null_data_sites_all_rapo$Lat <- scale(null_data_sites_all_rapo$Lat, center = T)
  null_data_sites_all_rapo$Long <- scale(null_data_sites_all_rapo$Long, center = T)
  
  sp_all_rapo <- SpatialPoints(data.frame(x = null_data_sites_all_rapo$Long, y = null_data_sites_all_rapo$Lat))
  crs(sp_all_rapo) <- "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"
  coords_all_rapo <- coordinates(sp_all_rapo)
  ###### MIN AND MAX NEAREST NEIGHBOR #####
  k_all_rapo <- knn2nb(knearneigh(coords_all_rapo, k = 1))
  dist_all_rapo <- unlist(nbdists(k_all_rapo, coords_all_rapo))
  max_all_rapo <- max(dist_all_rapo)
  ###### CREATE NEIGHBOR MATRIX (DISTANCE IN METERS) #####
  d_all_rapo <- dnearneigh(sp_all_rapo, longlat = FALSE, d1 = 0, d2 = max_all_rapo)
  ###### Create a data frame with path number and equations that will be tested in the pSEM #####
  path_equation_all_rapo <- data.frame(path_number = paste0("path", 1))
  path_equation_all_rapo$equation <-c("Range_size ~ Lat")
  ##### Create models for each path and weighting scheme #####
  dir_work_all_rapo <- "E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/All_rapo"
  setwd(dir_work_all_rapo)
  dir.create(paste0("sim_", j))
  dir_work_all_rapo_i <- paste0("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/All_rapo/sim_", j)
  setwd(dir_work_all_rapo_i)
  
  for (i in 1:length(scheme)){
    ###### Create a directory by scheme #####
    setwd(dir_work_all_rapo_i)
    dir_scheme <- scheme[i]
    if(!dir.exists(dir_scheme)) dir.create(dir_scheme)
    setwd(dir_scheme)
    ##### Create weighting scheme distance matrix #####
    spatial_weights_d_all_rapo <- nb2listw(d_all_rapo, zero.policy = TRUE, style = scheme[i])
    for (p in 1:nrow(path_equation_all_rapo)) {
      schemes <- rep(c(paste(scheme[i])), times = 2)
      paths_all_rapo <- rep(c(paste(path_equation_all_rapo[p, 1])), times = 2)
      models_all_rapo <- c("lm_mod", "error")
      results_path_all_rapo <- data.frame(scheme = schemes, path = paths_all_rapo, model = models_all_rapo)
      ###### OLS model #####
      lm_mod_all_rapo <- lm(paste(path_equation_all_rapo[p, 2]), data = null_data_sites_all_rapo)
      lm_mod_s_all_rapo <- summary(lm_mod_all_rapo)
      R2_lm_all_rapo <- lm_mod_s_all_rapo[["adj.r.squared"]]
      coef_lm_all_rapo <- lm_mod_all_rapo$coefficients
      p_value_lm_all_rapo <- glance(lm_mod_all_rapo)$p.value
      ##### SAR error models #####
      ##### Minimum distance #####
      error_d_all_rapo <- spatialreg::errorsarlm(paste(path_equation_all_rapo[p, 2]), data = null_data_sites_all_rapo, 
                                                 listw = spatial_weights_d_all_rapo, tol = 1e-12, zero.policy = T)
      error_d_s_all_rapo <-summary(error_d_all_rapo, Nagelkerke = TRUE)
      R2_error_d_all_rapo <- error_d_s_all_rapo$NK
      p_value_d_all_rapo <- error_d_s_all_rapo$Wald1$p.value
      ceof_d_all_rapo <- error_d_s_all_rapo$coefficients
      coef_all_rapo <- rbind(coef_lm_all_rapo, ceof_d_all_rapo)
      ##### Save R2 (pseudo R2 for errorsar) and AIC #####
      results_path_all_rapo$R2_models <- c(R2_lm_all_rapo, R2_error_d_all_rapo)
      results_path_all_rapo$AIC <- AIC(lm_mod_all_rapo, error_d_all_rapo)
      results_path_all_rapo$p_value <- c(p_value_lm_all_rapo, p_value_d_all_rapo)
      results_path_all_rapo <- cbind(results_path_all_rapo, coef_all_rapo)
      write.csv(results_path_all_rapo, file=paste(path_equation_all_rapo[p,1], ".csv"), row.names = F)
      ##### Make correlograms of residual autocorrelation #####
      cor.ols1.res_all_rapo <- correlog(as.numeric(null_data_sites_all_rapo$Long), as.numeric(null_data_sites_all_rapo$Lat), z = as.numeric(residuals(lm_mod_all_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
      cor.sar1.res_all_rapo <- correlog(as.numeric(null_data_sites_all_rapo$Long), as.numeric(null_data_sites_all_rapo$Lat), z = as.numeric(residuals(error_d_all_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
      ##### Save correlograms in a jpeg file #####
      jpeg(filename = paste(path_equation_all_rapo[p, 1], ".jpg"),
           width = 215, height = 279, units = "mm", res = 600)
      par(mfrow = c(2, 1))
      plot(cor.ols1.res_all_rapo, xlab = "Distance units", ylab = "Moran's I", ylim = c(-1,1), type = "l",
           lwd= 2, main = paste(models_all_rapo[1]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      plot(cor.sar1.res_all_rapo, xlab = "Distance units", ylab = "Moran's I", ylim=c(-1,1), type = "l",
           lwd= 2, main =paste(models_all_rapo[2]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      dev.off()
      ##### SAVE THE RESULTS #####
      save.image(file = paste(path_equation_all_rapo[p, 1], ".RData"))}}
}

###### READ THE SIMULATIONS #######
paths_sim <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/All_rapo",
                        pattern ="sim_",
                        full.names = TRUE)
sars_sim_list_All_rapo <- list()

##### SUMMARIZE THE RESULTS #####
for (i in 1:100){
  
  dir_work_All_rapo <- paths_sim[i]
  setwd(paste0(dir_work_All_rapo, "/S"))
  s_sata_All_rapo <- read.csv(paste0(dir_work_All_rapo, "/S/path1 .csv"))
  setwd(paste0(dir_work_All_rapo, "/C"))
  c_sata_All_rapo <- read.csv(paste0(dir_work_All_rapo, "/C/path1 .csv"))
  setwd(paste0(dir_work_All_rapo, "/W"))
  w_sata_All_rapo <- read.csv(paste0(dir_work_All_rapo, "/W/path1 .csv"))
  bind_data_All_rapo <- rbind(s_sata_All_rapo, c_sata_All_rapo, w_sata_All_rapo)
  bind_data_All_rapo$Simulation <- paste0("sim_", i)
  
  results_All_rapo <- bind_data_All_rapo[which.min(bind_data_All_rapo$AIC.AIC),]
  results_All_rapo1 <- results_All_rapo[, c(10, 9, 4, 5, 6, 7, 8, 1, 2, 3)]
  
  sars_sim_list_All_rapo[[i]] <- results_All_rapo1
}

sars_sim_list_All_rapo_1 <- as.data.frame(do.call(rbind, sars_sim_list_All_rapo))

##### Back to your working directory #####
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
##### Save the model selection results in a csv #####
write.csv(sars_sim_list_All_rapo_1, file="SAR_All_rapo.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####
##### TONINI #####
null_list_sites_ton_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton", 
                                       pattern = ".csv")
list_sar_sites_ton_rapo <- list()

###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")

###### MATRIX OF SPATIAL COORDINATES #####
for(j in 1:length(null_list_sites_ton_rapo)){
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Ton")
  null_data_sites_ton_rapo <- read.csv(null_list_sites_ton_rapo[j])
  
  null_data_sites_ton_rapo$Range_size <- scale(null_data_sites_ton_rapo$Range_size, center = T)
  null_data_sites_ton_rapo$Lat <- scale(null_data_sites_ton_rapo$Lat, center = T)
  null_data_sites_ton_rapo$Long <- scale(null_data_sites_ton_rapo$Long, center = T)
  
  sp_ton_rapo <- SpatialPoints(data.frame(x = null_data_sites_ton_rapo$Long, y = null_data_sites_ton_rapo$Lat))
  crs(sp_ton_rapo) <- "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"
  coords_ton_rapo <- coordinates(sp_ton_rapo)
  ###### MIN AND MAX NEAREST NEIGHBOR #####
  k_ton_rapo <- knn2nb(knearneigh(coords_ton_rapo, k = 1))
  dist_ton_rapo <- unlist(nbdists(k_ton_rapo, coords_ton_rapo))
  max_ton_rapo <- max(dist_ton_rapo)
  ###### CREATE NEIGHBOR MATRIX (DISTANCE IN METERS) #####
  d_ton_rapo <- dnearneigh(sp_ton_rapo, longlat = FALSE, d1 = 0, d2 = max_ton_rapo)
  ###### Create a data frame with path number and equations that will be tested in the pSEM #####
  path_equation_ton_rapo <- data.frame(path_number = paste0("path", 1))
  path_equation_ton_rapo$equation <-c("Range_size ~ Lat")
  ##### Create models for each path and weighting scheme #####
  dir_work_ton_rapo <- "E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Ton_rapo"
  setwd(dir_work_ton_rapo)
  dir.create(paste0("sim_", j))
  dir_work_ton_rapo_i <- paste0("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Ton_rapo/sim_", j)
  setwd(dir_work_ton_rapo_i)
  
  for (i in 1:length(scheme)){
    ###### Create a directory by scheme #####
    setwd(dir_work_ton_rapo_i)
    dir_scheme <- scheme[i]
    if(!dir.exists(dir_scheme)) dir.create(dir_scheme)
    setwd(dir_scheme)
    ##### Create weighting scheme distance matrix #####
    spatial_weights_d_ton_rapo <- nb2listw(d_ton_rapo, zero.policy = TRUE, style = scheme[i])
    for (p in 1:nrow(path_equation_ton_rapo)) {
      schemes <- rep(c(paste(scheme[i])), times = 2)
      paths_ton_rapo <- rep(c(paste(path_equation_ton_rapo[p, 1])), times = 2)
      models_ton_rapo <- c("lm_mod", "error")
      results_path_ton_rapo <- data.frame(scheme = schemes, path = paths_ton_rapo, model = models_ton_rapo)
      ###### OLS model #####
      lm_mod_ton_rapo <- lm(paste(path_equation_ton_rapo[p, 2]), data = null_data_sites_ton_rapo)
      lm_mod_s_ton_rapo <- summary(lm_mod_ton_rapo)
      R2_lm_ton_rapo <- lm_mod_s_ton_rapo[["adj.r.squared"]]
      coef_lm_ton_rapo <- lm_mod_ton_rapo$coefficients
      p_value_lm_ton_rapo <- glance(lm_mod_ton_rapo)$p.value
      ##### SAR error models #####
      ##### Minimum distance #####
      error_d_ton_rapo <- spatialreg::errorsarlm(paste(path_equation_ton_rapo[p, 2]), data = null_data_sites_ton_rapo, 
                                                 listw = spatial_weights_d_ton_rapo, tol = 1e-12, zero.policy = T)
      error_d_s_ton_rapo <-summary(error_d_ton_rapo, Nagelkerke = TRUE)
      R2_error_d_ton_rapo <- error_d_s_ton_rapo$NK
      p_value_d_ton_rapo <- error_d_s_ton_rapo$Wald1$p.value
      ceof_d_ton_rapo <- error_d_s_ton_rapo$coefficients
      coef_ton_rapo <- rbind(coef_lm_ton_rapo, ceof_d_ton_rapo)
      ##### Save R2 (pseudo R2 for errorsar) and AIC #####
      results_path_ton_rapo$R2_models <- c(R2_lm_ton_rapo, R2_error_d_ton_rapo)
      results_path_ton_rapo$AIC <- AIC(lm_mod_ton_rapo, error_d_ton_rapo)
      results_path_ton_rapo$p_value <- c(p_value_lm_ton_rapo, p_value_d_ton_rapo)
      results_path_ton_rapo <- cbind(results_path_ton_rapo, coef_ton_rapo)
      write.csv(results_path_ton_rapo, file=paste(path_equation_ton_rapo[p,1], ".csv"), row.names = F)
      ##### Make correlograms of residual autocorrelation #####
      cor.ols1.res_ton_rapo <- correlog(as.numeric(null_data_sites_ton_rapo$Long), as.numeric(null_data_sites_ton_rapo$Lat), z = as.numeric(residuals(lm_mod_ton_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
      cor.sar1.res_ton_rapo <- correlog(as.numeric(null_data_sites_ton_rapo$Long), as.numeric(null_data_sites_ton_rapo$Lat), z = as.numeric(residuals(error_d_ton_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
      ##### Save correlograms in a jpeg file #####
      jpeg(filename = paste(path_equation_ton_rapo[p, 1], ".jpg"),
           width = 215, height = 279, units = "mm", res = 600)
      par(mfrow = c(2, 1))
      plot(cor.ols1.res_ton_rapo, xlab = "Distance units", ylab = "Moran's I", ylim = c(-1,1), type = "l",
           lwd= 2, main = paste(models_ton_rapo[1]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      plot(cor.sar1.res_ton_rapo, xlab = "Distance units", ylab = "Moran's I", ylim=c(-1,1), type = "l",
           lwd= 2, main =paste(models_ton_rapo[2]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      dev.off()
      ##### SAVE THE RESULTS #####
      save.image(file = paste(path_equation_ton_rapo[p, 1], ".RData"))}}
}

###### READ THE SIMULATIONS #######
paths_sim <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/ton_rapo",
                        pattern ="sim_",
                        full.names = TRUE)
sars_sim_list_ton_rapo <- list()
##### SUMMARIZE THE RESULTS #####
for (i in 1:100){
  
  dir_work_ton_rapo <- paths_sim[i]
  setwd(paste0(dir_work_ton_rapo, "/S"))
  s_sata_ton_rapo <- read.csv(paste0(dir_work_ton_rapo, "/S/path1 .csv"))
  setwd(paste0(dir_work_ton_rapo, "/C"))
  c_sata_ton_rapo <- read.csv(paste0(dir_work_ton_rapo, "/C/path1 .csv"))
  setwd(paste0(dir_work_ton_rapo, "/W"))
  w_sata_ton_rapo <- read.csv(paste0(dir_work_ton_rapo, "/W/path1 .csv"))
  bind_data_ton_rapo <- rbind(s_sata_ton_rapo, c_sata_ton_rapo, w_sata_ton_rapo)
  bind_data_ton_rapo$Simulation <- paste0("sim_", i)
  
  results_ton_rapo <- bind_data_ton_rapo[which.min(bind_data_ton_rapo$AIC.AIC),]
  results_ton_rapo1 <- results_ton_rapo[, c(10, 9, 4, 5, 6, 7, 8, 1, 2, 3)]
  
  sars_sim_list_ton_rapo[[i]] <- results_ton_rapo1
}

sars_sim_list_ton_rapo_1 <- as.data.frame(do.call(rbind, sars_sim_list_ton_rapo))

##### Back to your working directory #####
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")
##### Save the model selection results in a csv #####
write.csv(sars_sim_list_ton_rapo_1, file="SAR_ton_rapo.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####

##### LEACHE #####
null_list_sites_lea_rapo <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea", 
                                       pattern = ".csv")
list_sar_sites_lea_rapo <- list()

###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")

###### MATRIX OF SPATIAL COORDINATES #####
for(j in 1:length(null_list_sites_lea_rapo)){
  setwd("E:/7_Doctorado/Cap_1/4_data/Null_models/sites_level/Lea")
  null_data_sites_lea_rapo <- read.csv(null_list_sites_lea_rapo[j])
  
  null_data_sites_lea_rapo$Range_size <- scale(null_data_sites_lea_rapo$Range_size, center = T)
  null_data_sites_lea_rapo$Lat <- scale(null_data_sites_lea_rapo$Lat, center = T)
  null_data_sites_lea_rapo$Long <- scale(null_data_sites_lea_rapo$Long, center = T)
  
  sp_lea_rapo <- SpatialPoints(data.frame(x = null_data_sites_lea_rapo$Long, y = null_data_sites_lea_rapo$Lat))
  crs(sp_lea_rapo) <- "+proj=aea +lat_1=5 +lat_2=60 +lon_0=-100 +units=m"
  coords_lea_rapo <- coordinates(sp_lea_rapo)
  ###### MIN AND MAX NEAREST NEIGHBOR #####
  k_lea_rapo <- knn2nb(knearneigh(coords_lea_rapo, k = 1))
  dist_lea_rapo <- unlist(nbdists(k_lea_rapo, coords_lea_rapo))
  max_lea_rapo <- max(dist_lea_rapo)
  ###### CREATE NEIGHBOR MATRIX (DISTANCE IN METERS) #####
  d_lea_rapo <- dnearneigh(sp_lea_rapo, longlat = FALSE, d1 = 0, d2 = max_lea_rapo)
  ###### Create a data frame with path number and equations that will be tested in the pSEM #####
  path_equation_lea_rapo <- data.frame(path_number = paste0("path", 1))
  path_equation_lea_rapo$equation <-c("Range_size ~ Lat")
  ##### Create models for each path and weighting scheme #####
  dir_work_lea_rapo <- "E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Lea_rapo"
  setwd(dir_work_lea_rapo)
  dir.create(paste0("sim_", j))
  dir_work_lea_rapo_i <- paste0("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Lea_rapo/sim_", j)
  setwd(dir_work_lea_rapo_i)
  
  for (i in 1:length(scheme)){
    ###### Create a directory by scheme #####
    setwd(dir_work_lea_rapo_i)
    dir_scheme <- scheme[i]
    if(!dir.exists(dir_scheme)) dir.create(dir_scheme)
    setwd(dir_scheme)
    ##### Create weighting scheme distance matrix #####
    spatial_weights_d_lea_rapo <- nb2listw(d_lea_rapo, zero.policy = TRUE, style = scheme[i])
    for (p in 1:nrow(path_equation_lea_rapo)) {
      schemes <- rep(c(paste(scheme[i])), times = 2)
      paths_lea_rapo <- rep(c(paste(path_equation_lea_rapo[p, 1])), times = 2)
      models_lea_rapo <- c("lm_mod", "error")
      results_path_lea_rapo <- data.frame(scheme = schemes, path = paths_lea_rapo, model = models_lea_rapo)
      ###### OLS model #####
      lm_mod_lea_rapo <- lm(paste(path_equation_lea_rapo[p, 2]), data = null_data_sites_lea_rapo)
      lm_mod_s_lea_rapo <- summary(lm_mod_lea_rapo)
      R2_lm_lea_rapo <- lm_mod_s_lea_rapo[["adj.r.squared"]]
      coef_lm_lea_rapo <- lm_mod_lea_rapo$coefficients
      p_value_lm_lea_rapo <- glance(lm_mod_lea_rapo)$p.value
      ##### SAR error models #####
      ##### Minimum distance #####
      error_d_lea_rapo <- spatialreg::errorsarlm(paste(path_equation_lea_rapo[p, 2]), data = null_data_sites_lea_rapo, 
                                                 listw = spatial_weights_d_lea_rapo, tol = 1e-12, zero.policy = T)
      error_d_s_lea_rapo <-summary(error_d_lea_rapo, Nagelkerke = TRUE)
      R2_error_d_lea_rapo <- error_d_s_lea_rapo$NK
      p_value_d_lea_rapo <- error_d_s_lea_rapo$Wald1$p.value
      ceof_d_lea_rapo <- error_d_s_lea_rapo$coefficients
      coef_lea_rapo <- rbind(coef_lm_lea_rapo, ceof_d_lea_rapo)
      ##### Save R2 (pseudo R2 for errorsar) and AIC #####
      results_path_lea_rapo$R2_models <- c(R2_lm_lea_rapo, R2_error_d_lea_rapo)
      results_path_lea_rapo$AIC <- AIC(lm_mod_lea_rapo, error_d_lea_rapo)
      results_path_lea_rapo$p_value <- c(p_value_lm_lea_rapo, p_value_d_lea_rapo)
      results_path_lea_rapo <- cbind(results_path_lea_rapo, coef_lea_rapo)
      write.csv(results_path_lea_rapo, file=paste(path_equation_lea_rapo[p,1], ".csv"), row.names = F)
      ##### Make correlograms of residual autocorrelation #####
      cor.ols1.res_lea_rapo <- correlog(as.numeric(null_data_sites_lea_rapo$Long), as.numeric(null_data_sites_lea_rapo$Lat), z = as.numeric(residuals(lm_mod_lea_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
      cor.sar1.res_lea_rapo <- correlog(as.numeric(null_data_sites_lea_rapo$Long), as.numeric(null_data_sites_lea_rapo$Lat), z = as.numeric(residuals(error_d_lea_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
      ##### Save correlograms in a jpeg file #####
      jpeg(filename = paste(path_equation_lea_rapo[p, 1], ".jpg"),
           width = 215, height = 279, units = "mm", res = 600)
      par(mfrow = c(2, 1))
      plot(cor.ols1.res_lea_rapo, xlab = "Distance units", ylab = "Moran's I", ylim = c(-1,1), type = "l",
           lwd= 2, main = paste(models_lea_rapo[1]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      plot(cor.sar1.res_lea_rapo, xlab = "Distance units", ylab = "Moran's I", ylim=c(-1,1), type = "l",
           lwd= 2, main =paste(models_lea_rapo[2]), cex.main = 2, cex.lab=1.8, cex.axis=1.5)
      abline(h=0, lty=5)
      dev.off()
      ##### SAVE THE RESULTS #####
      save.image(file = paste(path_equation_lea_rapo[p, 1], ".RData"))}}
}

###### READ THE SIMULATIONS #######
paths_sim <- list.files("E:/7_Doctorado/Cap_1/4_data/Null_models/SAR/Lea_rapo",
                        pattern ="sim_",
                        full.names = TRUE)
sars_sim_list_lea_rapo <- list()

##### SUMMARIZE THE RESULTS #####
for (i in 1:100){
  
  dir_work_lea_rapo <- paths_sim[i]
  setwd(paste0(dir_work_lea_rapo, "/S"))
  s_sata_lea_rapo <- read.csv(paste0(dir_work_lea_rapo, "/S/path1 .csv"))
  setwd(paste0(dir_work_lea_rapo, "/C"))
  c_sata_lea_rapo <- read.csv(paste0(dir_work_lea_rapo, "/C/path1 .csv"))
  setwd(paste0(dir_work_lea_rapo, "/W"))
  w_sata_lea_rapo <- read.csv(paste0(dir_work_lea_rapo, "/W/path1 .csv"))
  bind_data_lea_rapo <- rbind(s_sata_lea_rapo, c_sata_lea_rapo, w_sata_lea_rapo)
  bind_data_lea_rapo$Simulation <- paste0("sim_", i)
  
  results_lea_rapo <- bind_data_lea_rapo[which.min(bind_data_lea_rapo$AIC.AIC),]
  results_lea_rapo1 <- results_lea_rapo[, c(10, 9, 4, 5, 6, 7, 8, 1, 2, 3)]
  
  sars_sim_list_lea_rapo[[i]] <- results_lea_rapo1
}

sars_sim_list_lea_rapo_1 <- as.data.frame(do.call(rbind, sars_sim_list_lea_rapo))

##### Back to your working directory #####
setwd("E:/7_Doctorado/Cap_1/4_data/Null_models")

##### Save the model selection results in a csv #####
write.csv(sars_sim_list_lea_rapo_1, file="SAR_lea_rapo.csv", row.names = F)
save.image("E:/7_Doctorado/Cap_1/7_Codes/Rdata/6_Null_models_regresions.RData")
