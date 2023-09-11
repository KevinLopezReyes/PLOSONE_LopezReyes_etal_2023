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
##### CREATE A LIST TO SAVE THE ALPHA VALUES #####
list_alpha <- list()
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
#### SAVE THE ALPHA HULLS VALUES IN THE LIST #####
alpha_value <- range$alpha
alpha_value <- substr(alpha_value, 6, nchar(alpha_value))
list_alpha[i] <- alpha_value
setwd("E:/7_Doctorado/Cap_1/1_Records/Joint")
next}
##### SAVE THE ALPHA VALUES IN A DATAFRAME #####
list_alpha <- unlist(list_alpha)
Species <- substr(data,1,nchar(data)-4)
list_alpha <- cbind(Species, list_alpha)
setwd("E:/7_Doctorado/Cap_1/5_tables")
write.csv(list_alpha, "list_alpha.csv", row.names = FALSE)
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
##### I NEED TO CHANGE THE COLUMNS NAMES TO MATCH WITH THE ALPHA HULLS ####
##### BUT BEFORE YOU CAN DO THAT YOU NEED TO REMOVE UNIQUE RECORDS FROM THE CHULL1 FOLDER #####
mydir <- "E:/7_Doctorado/Cap_1/3_Models/chull1/"
delfiles <- dir(path = mydir, pattern="unique_records")
file.remove(file.path(mydir, delfiles))
##### NOW YOU CAN CHANGE THE COLUMNS NAMES #####
setwd("E:/7_Doctorado/Cap_1/3_Models/chull1")
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
setwd("E:/7_Doctorado/Cap_1/3_Models/chull")
writeOGR(mod_1, dsn=getwd(), layer = name, driver="ESRI Shapefile")
setwd("E:/7_Doctorado/Cap_1/3_Models/chull1")}
##### AFTER THIS, YOU CAN DELETE THE CHULL1 FOLDER #####
##### --------------------------------------------------------------------------------- #####
#############################################################################################
######################### CREATE THE ECOLOGICAL NICHE MODELING ##############################
#############################################################################################
##### WE ARE GONNA USE THE DATA CREATED IN THE "SCRIPT_1" #####
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
############################## SELECT THE BEST MODELS FOR MAXENT ############################
#############################################################################################
##### I HAVE SO MANY REPLICATES AND FINAL MODELS, SO I NEED TO SUMMARIZE ALL OF THEM #####
library(raster)
library(ntbox)
library(rangemap)
##### FIRST I GET ALL THE DIRECTORIES PER SPECIES WHERE ARE ALLOCATE THE FINAL MODELS #####
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species")
sps <- dir()
##### THE NEXT LOOP READ THE FOLDER NAMES, THEN EXTRACT THE NAME OF THE PATHS FOR EVERY FINAL MODEL #####
##### AND SAVE THE MEDIANS OF EACH FINAL MODEL IN A NEW FOLDER #####
for (i in 1:length(sps)){
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species")
names <- substr(sps[i],1,nchar(sps[i]))
path_cal <- paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", names, "/Calibration_results")
setwd(path_cal)
cal_res <- read.csv("selected_models.csv")
row_names <- cal_res[,1]
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species")
sps <- dir()
sps1 <- sps[i]
for (i in 1:length(row_names)){
path <- paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species/", sps1, "/Final_Models")
path_final <- paste0(path,"/", row_names[1],"_NE")
file <- sps1
models <- list.files(path_final,
                     pattern = paste0(file,"_current_median.asc"),
                     full.names = TRUE)
mod <- raster(models)
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/F_models")
dir.create(names)
p_m <- paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/F_models/", names)
setwd(p_m)
writeRaster(mod, filename=paste0(names, "_", row_names[i]), format = "ascii")}}
##### THE NEXT LOOP READ THE MODELS SAVED IN THE LAST LOOP AND CALCULATES THEIR MEDIANS ####
##### IF THERE ARE MORE THAN ONE FINAL MODEL (MEDIAN OF THE MEDIANS), IF ONLY ONE FINAL MODEL WERE SELECTED ####
#### IT ONLY TAKES THE MODELS AND SAVE IN THE FINAL FOLDER FOR THE BINARIZATION #####
for (i in 1:length(sps)){
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species")
sps <- dir()
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/species")
names <- substr(sps[i], 1, nchar(sps[i]))
p_m <- paste0("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/F_models/", names)
models <- list.files(p_m,
                     pattern = ".asc",
                     full.names = TRUE)
if(length(models) > 1){
  stack1 <- raster::stack(models)
  stack1 <- calc(stack1, median)
  } else {stack1 <- raster(models)}
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/Medians")
writeRaster(stack1, filename = paste0(names,"_median"), format="ascii")}
##### AFTER I HAVE THE FINAL MODELS (ONE MODEL PER SPECIES) I NEED TO MAKE IT BINARY AND SAVE IT AS SHAPEFILE #####
library(ntbox)
library(sf)
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
data <-list.files(pattern=".csv")
mods <- list.files("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/Medians",
                     pattern = ".asc",
                     full.names = TRUE)
					 
for(i in 1:length(data)){
setwd("E:/7_Doctorado/Cap_1/7_anexos/kuenm/bases/1_Joint")
spp <- read.csv(data[[i]], header=T, sep=",")
name <- substr(data[i],1,nchar(data[i])-4)
mod <- raster(mods[i])
setwd("E:/7_Doctorado/Cap_1/3_Models/kuenm")
bin_m <- bin_model(mod, spp[,c(2,3)], percent = 10)
pol <- rasterToPolygons(bin_m, dissolve = TRUE)
pol@data
pol$species = name
names(pol@data)[names(pol@data)=="species"] <- "BINOMIAL"
pol$layer
pol1 <- subset(pol, layer == "1")
writeOGR(pol1, dsn=  getwd(), layer = name, driver="ESRI Shapefile")}
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### END OF THE CODE ######################################
#############################################################################################


















