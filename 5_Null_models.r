
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
r_mask <- crop(rast_full, dist_ahull) ###### "dist_ahull" es un multipoligono que contiene las areas de distribución de todas las 
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

##### --------------------------------------------------------------------------------- #####
#############################################################################################
####################################### END OF THE CODE #####################################
#############################################################################################
##### --------------------------------------------------------------------------------- #####


