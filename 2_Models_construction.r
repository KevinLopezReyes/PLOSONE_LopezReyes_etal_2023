#############################################################################################
################################### LOAD THE PACKAGES #######################################
#############################################################################################
library(rangeBuilder)
library(raster)
library(maps)
library(maptools)
library(rgdal)
library(rangemap)
library(kuenm)
library(sf)
library(nngeo)
library(smoothr)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################# CREATE THE ALPHA HULLS ####################################
#############################################################################################
##### LOAD THE PRESENCE RECORDS
setwd("E:/7_Doctorado/Cap_1/1_Records/Joint")
data <- list.files (pattern = ".csv")
##### RUN THE FUNCTION TO CREATE THE ALPHA HULLS #####
for(i in 1:length(data)){
spp <- read.csv(data[[i]], header = T, sep = ",")
name <- substr(data[i],1,nchar(data[i])-4)
range <- getDynamicAlphaHull(spp, fraction = 1, partCount = 1,
buff = 50000, initialAlpha = 0, coordHeaders = c("Long", "Lat"),
clipToCoast = "terrestrial", alphaIncrement = 0.5, verbose = FALSE,
alphaCap = 10)
##### MAKE THE POLYGON SPATIAL #####
range_df <- as(range[[1]], "Spatial")
##### CHANGE THE COLUMNS NAMES #####
range_df$dummy = name
range_df@data
names(range_df@data)[names(range_df@data)=="dummy"] <- "BINOMIAL"
##### REMOVE THE SMALL POLYGONS #####
range_df <- disaggregate(range_df)
range_df$areakm2 <- area(range_df)
max <- max(range_df$areakm2) / 2
range_df <- drop_crumbs(range_df, threshold = max, drop_empty = TRUE)
##### SAVE THE RESULTS #####
setwd("E:/7_Doctorado/Cap_1/3_Models/ahull")
writeOGR(range_df, dsn = getwd(), layer = name, driver = "ESRI Shapefile")
setwd("E:/7_Doctorado/Cap_1/1_Records/Joint")
next}

##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################ CREATE THE CONVEX HULLS ####################################
#############################################################################################
##### LOAD THE PRESENCE RECORDS #####
setwd("E:/7_Doctorado/Cap_1/1_Records/Joint")
data <- list.files (pattern = ".csv")
##### RUN THE CODE FOR CREATE THE CONVEX HULLS #####
for(i in 1:length(data)){
spp <- read.csv(data[[i]], header=T, sep=",")
spp_name <- substr(data[i],1,nchar(data[i])-4)
setwd("E:/7_Doctorado/Cap_1/3_Models/chull1")
chull1 <- rangemap_hull(occurrences = spp, hull_type = "convex", buffer_distance = 50000, 
						extent_of_occurrence = FALSE, area_of_occupancy = FALSE,
						save_shp = TRUE, name = spp_name,
						overwrite = TRUE)
setwd("E:/7_Doctorado/Cap_1/1_Records/Joint")
next}
##### --------------------------------------------------------------------------------- #####
#############################################################################################
######################### CREATE THE ECOLOGICAL NICHE MODELING ##############################
#############################################################################################
##### SET THW WORKING DIRECTORY AND THE NAMES OF THE OBJECTS AND PATHS THAT I NEED #####
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species")
sps <- dir()
occ_joint <- paste0(sps, "/occ_joint.csv")
occ_tra <- paste0(sps, "/occ_train.csv")
occ_test <- paste0(sps, "/occ_test.csv")
back_folder <- paste0(sps, "/Background")
batch_cal <- paste0(sps, "/batch_calibration")
out_candidates <- paste0(sps, "/Candidate_models")
out_calibration <- paste0(sps, "/Calibration_results")
regm <- seq(1,4,1)
fclas <- "basic"
mxpath <- "E:/7_Doctorado/Cap_1/7_anexos/kuenm"
##### WITH THE NEXT LOOP, I RUN THE CALIBRATION PROCESS #####
for (i in 1:length(sps)){
kuenm_cal_swd(occ.joint = occ_joint[i], occ.tra = occ_tra[i], occ.test = occ_test[i], 
				back.dir = back_folder[i], batch = batch_cal[i],
              out.dir.models = out_candidates[i], reg.mult = regm, f.clas = fclas,
              max.memory = 1000, args = NULL, maxent.path = mxpath,
              selection = "OR_AICc", threshold = 10,
              rand.percent = 50, iterations = 500,
              kept = TRUE, out.dir.eval = out_calibration[i])
next}
##### AFTER THE CALIBRATION PROCESS FINISH, I NEED TO RUN THE EVALUATION AND THE CONSTRUCTION OF THE FINAL MODEL #####
batch_fin <- paste0(sps, "/batch_models")
g_variables <- paste0(sps, "/G_var")
out_models <- paste0(sps, "/Final_models")
for (i in 1:length(sps)){
kuenm_mod_swd(occ.joint = occ_joint[i], back.dir = back_folder[i], out.eval = out_calibration[i], 
			  batch = batch_fin[i], rep.n = 10, rep.type = "Bootstrap",
			  jackknife = FALSE, max.memory = 1000, out.format = "cloglog", project = TRUE, G.var.dir = g_variables[i],
              ext.type = "no_ext", write.mess = FALSE, write.clamp = FALSE, maxent.path = mxpath,
			  args = NULL, out.dir = out_models[i], wait = TRUE, run = TRUE)}
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### END OF THE CODE ######################################
#############################################################################################
##### --------------------------------------------------------------------------------- #####

















