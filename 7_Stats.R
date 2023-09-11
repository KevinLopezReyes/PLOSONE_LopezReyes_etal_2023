#############################################################################################
########################################## LOAD AND PACKAGES ################################
#############################################################################################
library(MuMIn)
library(ggplot2)
library(tibble)
library(ggpubr)
library(usdm)
library(tidyr)
library(spdep)
library(ncf)
library(broom)
#############################################################################################
################################## LOAD AND TRANSFORM THE DATA ##############################
#############################################################################################
##### LOAD THE DATA #####
##### SPECIES LEVEL #####
sp_vars_all <- read.csv("G:/7_Doctorado/Cap_1/4_data/all_data_specieslevel.csv")
sp_vars_ton <- read.csv("G:/7_Doctorado/Cap_1/4_data/ton_data_specieslevel.csv")
sp_vars_lea <- read.csv("G:/7_Doctorado/Cap_1/4_data/lea_data_specieslevel.csv")
##### SITES LEVEL #####
sites_vars_all <- read.csv("G:/7_Doctorado/Cap_1/4_data/all_data_siteslevel.csv")
sites_vars_ton <- read.csv("G:/7_Doctorado/Cap_1/4_data/ton_data_siteslevel.csv")
sites_vars_lea <- read.csv("G:/7_Doctorado/Cap_1/4_data/lea_data_siteslevel.csv")
##### TRANSFORM THE DATA TO MAKE THE COMPARABLE COEFFICIENTS #####
##### SPECIES #####
##### ALL #####
sp_vars_all$Long <- scale(sp_vars_all$Long, center = T)
sp_vars_all$Lat <- scale(sp_vars_all$Lat, center = T)
sp_vars_all$Range_size <- scale(sp_vars_all$Range_size, center = T)
sp_vars_all$CCV <- scale(sp_vars_all$CCV, center = T)
sp_vars_all$Elevation  <- scale(sp_vars_all$Elevation, center = T)
sp_vars_all$CEH_min <- scale(sp_vars_all$CEH_min, center = T)
sp_vars_all$CEH_max <- scale(sp_vars_all$CEH_max, center = T)
sp_vars_all$CVH <- scale(sp_vars_all$CVH, center = T)
##### TONINI ######
sp_vars_ton$Long <- scale(sp_vars_ton$Long, center = T)
sp_vars_ton$Lat <- scale(sp_vars_ton$Lat, center = T)
sp_vars_ton$Range_size <- scale(sp_vars_ton$Range_size, center = T)
sp_vars_ton$CCV <- scale(sp_vars_ton$CCV, center = T)
sp_vars_ton$Elevation  <- scale(sp_vars_ton$Elevation, center = T)
sp_vars_ton$CEH_min <- scale(sp_vars_ton$CEH_min, center = T)
sp_vars_ton$CEH_max <- scale(sp_vars_ton$CEH_max, center = T)
sp_vars_ton$CVH <- scale(sp_vars_ton$CVH, center = T)
##### LEACHE #####
sp_vars_lea$Long <- scale(sp_vars_lea$Long, center = T)
sp_vars_lea$Lat <- scale(sp_vars_lea$Lat, center = T)
sp_vars_lea$Range_size <- scale(sp_vars_lea$Range_size, center = T)
sp_vars_lea$CCV <- scale(sp_vars_lea$CCV, center = T)
sp_vars_lea$Elevation  <- scale(sp_vars_lea$Elevation, center = T)
sp_vars_lea$CEH_min <- scale(sp_vars_lea$CEH_min, center = T)
sp_vars_lea$CEH_max <- scale(sp_vars_lea$CEH_max, center = T)
sp_vars_lea$CVH <- scale(sp_vars_lea$CVH, center = T)
##### SITES #####
##### ALL #####
sites_vars_all$Long <- scale(sites_vars_all$Long, center = T)
sites_vars_all$Lat <- scale(sites_vars_all$Lat, center = T)
sites_vars_all$Range_size <- scale(sites_vars_all$Range_size, center = T)
sites_vars_all$CCV <- scale(sites_vars_all$CCV, center = T)
sites_vars_all$Elevation  <- scale(sites_vars_all$Elevation, center = T)
sites_vars_all$CEH_min <- scale(sites_vars_all$CEH_min, center = T)
sites_vars_all$CEH_max <- scale(sites_vars_all$CEH_max, center = T)
sites_vars_all$CVH <- scale(sites_vars_all$CVH, center = T)
##### TONINI ######
sites_vars_ton$Long <- scale(sites_vars_ton$Long, center = T)
sites_vars_ton$Lat <- scale(sites_vars_ton$Lat, center = T)
sites_vars_ton$Range_size <- scale(sites_vars_ton$Range_size, center = T)
sites_vars_ton$CCV <- scale(sites_vars_ton$CCV, center = T)
sites_vars_ton$Elevation  <- scale(sites_vars_ton$Elevation, center = T)
sites_vars_ton$CEH_min <- scale(sites_vars_ton$CEH_min, center = T)
sites_vars_ton$CEH_max <- scale(sites_vars_ton$CEH_max, center = T)
sites_vars_ton$CVH <- scale(sites_vars_ton$CVH, center = T)
##### LEACHE #####
sites_vars_lea$Long <- scale(sites_vars_lea$Long, center = T)
sites_vars_lea$Lat <- scale(sites_vars_lea$Lat, center = T)
sites_vars_lea$Range_size <- scale(sites_vars_lea$Range_size, center = T)
sites_vars_lea$CCV <- scale(sites_vars_lea$CCV, center = T)
sites_vars_lea$Elevation  <- scale(sites_vars_lea$Elevation, center = T)
sites_vars_lea$CEH_min <- scale(sites_vars_lea$CEH_min, center = T)
sites_vars_lea$CEH_max <- scale(sites_vars_lea$CEH_max, center = T)
sites_vars_lea$CVH <- scale(sites_vars_lea$CVH, center = T)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
######################################### PHYLOGENIES #######################################
#############################################################################################
library(caper)
leache <- read.nexus("C:/Users/Kevin/Downloads/Filogenias/leache.nex")
tonini <- read.nexus("C:/Users/Kevin/Downloads/Filogenias/tonini.nex")
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################ CHECK COLLINEARITY AMONG PREDICTORS ########################
#############################################################################################
options(scipen=999)
##### SPECIES LEVEL - VIF #####
vif(sp_vars_all[,c(5:9)])
vif(sp_vars_ton[,c(5:9)])
vif(sp_vars_ton[,c(5:9)])
###### SITES LEVEL - VIF #####
vif(sites_vars_all[4:8])
vif(sites_vars_ton[4:8])
vif(sites_vars_ton[4:8])
##### --------------------------------------------------------------------------------- #####
#############################################################################################
############################# PERFORM THE OLS AND PGLS (SPECIES) ############################
#############################################################################################
##### OLS - SPECIES LEVEL #####
##### ALL SPECIES #####
options(na.action = "na.fail")
ols_sp_full <- lm(Range_size ~ CCV + Elevation + CVH + CEH_min, sp_vars_all)
ols_sp_full_dg <- dredge(ols_sp_full, rank = "AIC")
ols_sp_full_dg
ols_sp_full_sum <- summary(model.avg(ols_sp_full_dg, subset = cumsum(weight) <= .95))
ols_sp_full_sum <- as.data.frame(ols_sp_full_sum$coefmat.full)
ols_sp_full_bm <- get.models(ols_sp_full_dg, 1)[[1]]
ols_sp_full_bm_s  <- as.data.frame(broom::glance(ols_sp_full_bm))
ols_sp_full_bm_s$Set <- "All"
##### --------------------------------------------------------------------------------- #####
##### TONINI #####
ols_sp_ton <- lm(Range_size ~ CCV + Elevation + CEH_min + CVH, sp_vars_ton)
ols_sp_ton_dg <- dredge(ols_sp_ton, rank = "AIC")
ols_sp_ton_dg
ols_sp_ton_sum <- summary(model.avg(ols_sp_ton_dg, subset = cumsum(weight) <= .95))
ols_sp_ton_sum <- as.data.frame(ols_sp_ton_sum$coefmat.full)
ols_sp_ton_bm <- get.models(ols_sp_ton_dg, 1)[[1]]
ols_sp_ton_bm_s  <- as.data.frame(broom::glance(ols_sp_ton_bm))
ols_sp_ton_bm_s$Set <- "Tonini"
##### --------------------------------------------------------------------------------- #####
##### LEACHE #####
ols_sp_lea <- lm(Range_size ~ CCV + Elevation + CEH_min + CVH, sp_vars_lea)
ols_sp_lea_dg <- dredge(ols_sp_lea, rank = "AIC")
ols_sp_lea_dg
ols_sp_lea_sum <- summary(model.avg(ols_sp_lea_dg, subset = cumsum(weight) <= .95))
ols_sp_lea_sum <- as.data.frame(ols_sp_lea_sum$coefmat.full)
ols_sp_lea_bm <- get.models(ols_sp_lea_dg, 1)[[1]]
ols_sp_lea_bm_s  <- as.data.frame(broom::glance(ols_sp_lea_bm))
ols_sp_lea_bm_s$Set <- "Leache"
##### SAVE ALL #####
ols_sp_sum_all <- cbind(ols_sp_full_sum, ols_sp_ton_sum, ols_sp_lea_sum)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(ols_sp_sum_all, "ols_sp_sum_all.csv")
##### SAVE DREDGE #####
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(ols_sp_full_dg, "ols_sp_full_dg.csv", row.names = FALSE)
write.csv(ols_sp_ton_dg, "ols_sp_ton_dg.csv", row.names = FALSE)
write.csv(ols_sp_lea_dg, "ols_sp_lea_dg.csv", row.names = FALSE)
##### BEST MODEL PARAMETERS #####
ols_bm_sp <- rbind(ols_sp_full_bm_s, ols_sp_ton_bm_s, ols_sp_lea_bm_s)
write.csv(ols_bm_sp, "ols_bm_sp.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
##### PGLS #####
##### COMPARATIVE DATA #####
comp_ton <- comparative.data(tonini, data = sp_vars_ton, names.col = "Species", vcv = T)
comp_lea <- comparative.data(leache, data = sp_vars_lea, names.col = "Species", vcv = T)
##### PERFORM THE PGLS ANALYSIS #####
##### TONINI #####
pgls_sp_ton <- pgls(Range_size ~ CCV + Elevation + CEH_min + CVH, data = comp_ton, lambda="ML")
pgls_sp_ton_dg <- dredge(pgls_sp_ton, rank = "AIC")
pgls_sp_ton_dg
pgls_sp_ton_sum <- summary(model.avg(pgls_sp_ton_dg, subset = cumsum(weight) <= .95))
pgls_sp_ton_sum <- as.data.frame(pgls_sp_ton_sum$coefmat.full)
pgls_sp_ton_bm <- summary(get.models(pgls_sp_ton_dg, 1)[[1]])
pgls_sp_ton_bm
##### LEACHE #####
pgls_sp_lea <- pgls(Range_size ~ CCV + Elevation + CEH_min + CVH, data = comp_lea, lambda="ML")
pgls_sp_lea_dg <- dredge(pgls_sp_lea, rank = "AIC")
pgls_sp_lea_dg
pgls_sp_lea_sum <- summary(model.avg(pgls_sp_lea_dg, subset = cumsum(weight) <= .95))
pgls_sp_lea_sum <- as.data.frame(pgls_sp_lea_sum$coefmat.full)
pgls_sp_lea_bm <- summary(get.models(pgls_sp_lea_dg, 1)[[1]])
pgls_sp_lea_bm
##### SAVE ALL #####
pgls_sp_sum_all <- cbind(pgls_sp_ton_sum, pgls_sp_lea_sum)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(pgls_sp_sum_all, "pgls_sum_all.csv")
##### SAVE ALL DREDGE #####
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(pgls_sp_ton_dg, "pgls_sp_ton_dg.csv", row.names = FALSE)
write.csv(pgls_sp_lea_dg, "pgls_sp_lea_dg.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
############################### PERFORM THE OLS AND SARS (SITES) ############################
#############################################################################################
##### OLS - SITES LEVEL #####
##### ALL SPECIES - SITES #####
ols_sites_full <- lm(Range_size ~ CCV + Elevation + CEH_min + CVH, sites_vars_all)
ols_sites_full_dg <- dredge(ols_sites_full, rank = "AIC")
ols_sites_full_dg
ols_sites_full_fmod <- summary(get.models(ols_sites_full_dg, 1)[[1]])$coefficients
ols_sites_full_bm <- get.models(ols_sites_full_dg, 1)[[1]]
ols_sites_full_bm_s  <- as.data.frame(broom::glance(ols_sites_full_bm))
ols_sites_full_bm_s$Set <- "All"
##### --------------------------------------------------------------------------------- #####
##### TONINI #####
ols_sites_ton <- lm(Range_size ~ CCV + Elevation + CEH_min + CVH, sites_vars_ton)
ols_sites_ton_dg <- dredge(ols_sites_ton, rank = "AIC")
ols_sites_ton_dg
ols_sites_ton_fmod <- summary(get.models(ols_sites_ton_dg, 1)[[1]])$coefficients
ols_sites_ton_bm <- get.models(ols_sites_ton_dg, 1)[[1]]
ols_sites_ton_bm_s  <- as.data.frame(broom::glance(ols_sites_ton_bm))
ols_sites_ton_bm_s$Set <- "Tonini"
##### LEACHE #####
ols_sites_lea <- lm(Range_size ~ CCV + Elevation + CEH_min + CVH, sites_vars_lea)
ols_sites_lea_dg <- dredge(ols_sites_lea, rank = "AIC")
ols_sites_lea_dg
ols_sites_lea_fmod <- summary(get.models(ols_sites_lea_dg, 1)[[1]])$coefficients
ols_sites_lea_bm <- get.models(ols_sites_lea_dg, 1)[[1]]
ols_sites_lea_bm_s  <- as.data.frame(broom::glance(ols_sites_lea_bm))
ols_sites_lea_bm_s$Set <- "Leache"
##### SAVE ALL #####
ols_sites_sum_all <- cbind(ols_sites_full_fmod, ols_sites_ton_fmod, ols_sites_lea_fmod)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(ols_sites_sum_all, "ols_sites_sum_all.csv")
##### SAVE ALL DREDGE #####
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(ols_sites_full_dg, "ols_sites_full_dg.csv", row.names = FALSE)
write.csv(ols_sites_ton_dg, "ols_sites_ton_dg.csv", row.names = FALSE)
write.csv(ols_sites_lea_dg, "ols_sites_lea_dg.csv", row.names = FALSE)
##### BEST MODEL PARAMETERS #####
ols_bm_sites <- rbind(ols_sites_full_bm_s, ols_sites_ton_bm_s, ols_sites_lea_bm_s)
write.csv(ols_bm_sites, "ols_bm_sites.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################## PERFORM THE SARS (SITES) #################################
#############################################################################################
library(spdep)
library(ncf)
library(broom)
###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")
##### --------------------------------------------------------------------------------- #####
##### ALL SPECIES #####
###### MATRIX OF SPATIAL COORDINATES #####
sp_all <- SpatialPoints(data.frame(x = sites_vars_all$Long, y = sites_vars_all$Lat))
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
dir_work_all <- "G:/7_Doctorado/Cap_1/4_data/SAR/All"

for (i in 1:length(scheme)){
  ###### Create a directory by scheme #####
  setwd(dir_work_all)
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
    lm_mod_all <- lm(paste(path_equation_all[p, 2]), data = sites_vars_all)
    lm_mod_s_all <- summary(lm_mod_all)
    R2_lm_all <- lm_mod_s_all[["adj.r.squared"]]
    coef_lm_all <- lm_mod_all$coefficients
    p_value_lm_all <- glance(lm_mod_all)$p.value
    ##### SAR error models #####
    ##### Minimum distance #####
    error_d_all <- spatialreg::errorsarlm(paste(path_equation_all[p, 2]), data = sites_vars_all, 
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
    cor.ols1.res_all <- correlog(as.numeric(sites_vars_all$Long), as.numeric(sites_vars_all$Lat), z = as.numeric(residuals(lm_mod_all)), na.rm = TRUE, increment = 1, resamp = 1)
    cor.sar1.res_all <- correlog(as.numeric(sites_vars_all$Long), as.numeric(sites_vars_all$Lat), z = as.numeric(residuals(error_d_all)), na.rm = TRUE, increment = 1, resamp = 1)
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
##### SUMMARIZE THE RESULTS #####
results_all <- data.frame()
for (i in 1:length(scheme)){
  setwd(dir_work_all)
  ##### Extract results from each scheme directory #####
  dir_scheme <- scheme[i]
  setwd(dir_scheme)
  scheme_results_list_all <- list.files(path = ".", pattern= '.csv$')
  ##### Join all results #####
  scheme_results_all <- lapply(scheme_results_list_all, read.csv)
  scheme_results2_all <- do.call("rbind", scheme_results_all)
  results_all <- rbind(results_all, scheme_results2_all)}
##### Back to your working directory #####
setwd(dir_work_all)
##### Save the model selection results in a csv #####
write.csv(results_all, file="results_model_selection_all.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####

##### TONINI #####
###### MATRIX OF SPATIAL COORDINATES #####
sp_ton <- SpatialPoints(data.frame(x = sites_vars_ton$Long, y = sites_vars_ton$Lat))
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
dir_work_ton <- "G:/7_Doctorado/Cap_1/4_data/SAR/Ton"

for (i in 1:length(scheme)){
  ###### Create a directory by scheme #####
  setwd(dir_work_ton)
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
    lm_mod_ton <- lm(paste(path_equation_ton[p, 2]), data = sites_vars_ton)
    lm_mod_s_ton <- summary(lm_mod_ton)
    R2_lm_ton <- lm_mod_s_ton[["adj.r.squared"]]
    coef_lm_ton <- lm_mod_ton$coefficients
    p_value_lm_ton <- glance(lm_mod_ton)$p.value
    ##### SAR error models #####
    ##### Minimum distance #####
    error_d_ton <- spatialreg::errorsarlm(paste(path_equation_ton[p, 2]), data = sites_vars_ton, 
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
    cor.ols1.res_ton <- correlog(as.numeric(sites_vars_ton$Long), as.numeric(sites_vars_ton$Lat), z = as.numeric(residuals(lm_mod_ton)), na.rm = TRUE, increment = 1, resamp = 1)
    cor.sar1.res_ton <- correlog(as.numeric(sites_vars_ton$Long), as.numeric(sites_vars_ton$Lat), z = as.numeric(residuals(error_d_ton)), na.rm = TRUE, increment = 1, resamp = 1)
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
##### SUMMARIZE THE RESULTS #####
results_ton <- data.frame()
for (i in 1:length(scheme)){
  setwd(dir_work_ton)
  ##### Extract results from each scheme directory #####
  dir_scheme <- scheme[i]
  setwd(dir_scheme)
  scheme_results_list_ton <- list.files(path = ".", pattern= '.csv$')
  ##### Join all results #####
  scheme_results_ton <- lapply(scheme_results_list_ton, read.csv)
  scheme_results2_ton <- do.call("rbind", scheme_results_ton)
  results_ton <- rbind(results_ton, scheme_results2_ton)}
##### Back to your working directory #####
setwd(dir_work_ton)
##### Save the model selection results in a csv #####
write.csv(results_ton, file="results_model_selection_ton.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####

##### LEACHE #####
###### MATRIX OF SPATIAL COORDINATES #####
sp_lea <- SpatialPoints(data.frame(x = sites_vars_lea$Long, y = sites_vars_lea$Lat))
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
dir_work_lea <- "G:/7_Doctorado/Cap_1/4_data/SAR/Lea"

for (i in 1:length(scheme)){
  ###### Create a directory by scheme #####
  setwd(dir_work_lea)
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
    lm_mod_lea <- lm(paste(path_equation_lea[p, 2]), data = sites_vars_lea)
    lm_mod_s_lea <- summary(lm_mod_lea)
    R2_lm_lea <- lm_mod_s_lea[["adj.r.squared"]]
    coef_lm_lea <- lm_mod_lea$coefficients
    p_value_lm_lea <- glance(lm_mod_lea)$p.value
    ##### SAR error models #####
    ##### Minimum distance #####
    error_d_lea <- spatialreg::errorsarlm(paste(path_equation_lea[p, 2]), data = sites_vars_lea, 
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
    cor.ols1.res_lea <- correlog(as.numeric(sites_vars_lea$Long), as.numeric(sites_vars_lea$Lat), z = as.numeric(residuals(lm_mod_lea)), na.rm = TRUE, increment = 1, resamp = 1)
    cor.sar1.res_lea <- correlog(as.numeric(sites_vars_lea$Long), as.numeric(sites_vars_lea$Lat), z = as.numeric(residuals(error_d_lea)), na.rm = TRUE, increment = 1, resamp = 1)
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
##### SUMMARIZE THE RESULTS #####
results_lea <- data.frame()
for (i in 1:length(scheme)){
  setwd(dir_work_lea)
  ##### Extract results from each scheme directory #####
  dir_scheme <- scheme[i]
  setwd(dir_scheme)
  scheme_results_list_lea <- list.files(path = ".", pattern= '.csv$')
  ##### Join all results #####
  scheme_results_lea <- lapply(scheme_results_list_lea, read.csv)
  scheme_results2_lea <- do.call("rbind", scheme_results_lea)
  results_lea <- rbind(results_lea, scheme_results2_lea)}
##### Back to your working directory #####
setwd(dir_work_lea)
##### Save the model selection results in a csv #####
write.csv(results_lea, file="results_model_selection_lea.csv", row.names = F)

save.image("G:/7_Doctorado/Cap_1/7_Codes/Rdata/5_stats_plots.RData")
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### REGRESSIONS (RAPOPORT) ###############################
#############################################################################################
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################## OLS SPECIES LEVEL (RAPOPORT) #############################
#############################################################################################
###### ALL SPECIES ######
ols_sp_rapo_all <- lm(Range_size ~ Lat, sp_vars_all)
ols_sp_rapo_sum_all <- summary(ols_sp_rapo_all)$coefficients
ols_sp_rapo_est_all <- as.data.frame(ols_sp_rapo_sum_all[,1])
##### GET THE ESTIMATES IN TO COLUMNS #######
ols_sp_rapo_est_all <- tibble::rownames_to_column(ols_sp_rapo_est_all, "Hypothesis")
colnames(ols_sp_rapo_est_all) <- c("Hypothesis", "Estimate")
###### R2 AND P VALUE #####
##### GET THEIR COEFFICIENTS #######
ols_sp_rapo_bm_all  <- as.data.frame(broom::glance(ols_sp_rapo_all))
ols_sp_rapo_bm_s_all <- ols_sp_rapo_bm_all[,c(2,5,8,9)]  

ols_sp_rapo_all_wider <- ols_sp_rapo_est_all %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
ols_sp_rapo_all_wider_df <- as.data.frame(ols_sp_rapo_all_wider[,-1])
ols_sp_rapo_all_wider_df

final_ols_df_sp_all_rapo <- cbind(ols_sp_rapo_all_wider_df, ols_sp_rapo_bm_s_all)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(final_ols_df_sp_all_rapo, "OLS_rapoport_All_specieslevel.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
######## TONINI #########
sp_vars_tonols_sp_rapo_ton <- lm(Range_size ~ Lat, sp_vars_ton)
ols_sp_rapo_sum_ton <- summary(ols_sp_rapo_ton)$coefficients
ols_sp_rapo_est_ton <- as.data.frame(ols_sp_rapo_sum_ton[,1])
##### GET THE ESTIMATES IN TO COLUMNS #######
ols_sp_rapo_est_ton <- tibble::rownames_to_column(ols_sp_rapo_est_ton, "Hypothesis")
colnames(ols_sp_rapo_est_ton) <- c("Hypothesis", "Estimate")
###### R2 AND P VALUE #####
##### GET THEIR COEFFICIENTS #######
ols_sp_rapo_bm_ton  <- as.data.frame(broom::glance(ols_sp_rapo_ton))
ols_sp_rapo_bm_s_ton <- ols_sp_rapo_bm_ton[,c(2,5,8,9)]  

ols_sp_rapo_ton_wider <- ols_sp_rapo_est_ton %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
ols_sp_rapo_ton_wider_df <- as.data.frame(ols_sp_rapo_ton_wider[,-1])
ols_sp_rapo_ton_wider_df

final_ols_df_sp_ton_rapo <- cbind(ols_sp_rapo_ton_wider_df, ols_sp_rapo_bm_s_ton)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(final_ols_df_sp_ton_rapo, "OLS_rapoport_ton_specieslevel.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
######## LEACHE #########
ols_sp_rapo_lea <- lm(Range_size ~ Lat, sp_vars_lea)
ols_sp_rapo_sum_lea <- summary(ols_sp_rapo_lea)$coefficients
ols_sp_rapo_est_lea <- as.data.frame(ols_sp_rapo_sum_lea[,1])
##### GET THE ESTIMATES IN TO COLUMNS #######
ols_sp_rapo_est_lea <- tibble::rownames_to_column(ols_sp_rapo_est_lea, "Hypothesis")
colnames(ols_sp_rapo_est_lea) <- c("Hypothesis", "Estimate")
###### R2 AND P VALUE #####
##### GET THEIR COEFFICIENTS #######
ols_sp_rapo_bm_lea  <- as.data.frame(broom::glance(ols_sp_rapo_lea))
ols_sp_rapo_bm_s_lea <- ols_sp_rapo_bm_lea[,c(2,5,8,9)]  

ols_sp_rapo_lea_wider <- ols_sp_rapo_est_lea %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
ols_sp_rapo_lea_wider_df <- as.data.frame(ols_sp_rapo_lea_wider[,-1])
ols_sp_rapo_lea_wider_df

final_ols_df_sp_lea_rapo <- cbind(ols_sp_rapo_lea_wider_df, ols_sp_rapo_bm_s_lea)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(final_ols_df_sp_lea_rapo, "OLS_rapoport_lea_specieslevel.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################## PGLS SPECIES LEVEL (RAPOPORT) ############################
#############################################################################################
######### TONINI #########
pgls_sp_ton_rapo <- pgls(Range_size ~ Lat, data = comp_ton, lambda = "ML")
pgls_sp_ton_rapo_sum <- summary(pgls_sp_ton_rapo)
pgls_sp_ton_rapo_sum <- as.data.frame(pgls_sp_ton_rapo_sum$coefficients)
pgls_sp_ton_rapo_est <- pgls_sp_ton_rapo_sum[1]

###### R2 AND P VALUE #####
###### GET THE MODEL WITH THE LOWEST AIC #######
pgls_sp_ton_rapo_bm <- summary(pgls_sp_ton_rapo)
pgls_sp_ton_rapo_bm

##### GET THEIR COEFFICIENTS #######
pgls_ton_rapo_r2 <- pgls_sp_ton_rapo_bm$adj.r.squared
pgls_ton_rapo_p <- pgls_sp_ton_rapo_bm$coefficients[2,4]
pgls_ton_rapo_r2_p <- cbind(pgls_ton_rapo_r2, pgls_ton_rapo_p)
pgls_ton_rapo_r2_p_df <- as.data.frame(pgls_ton_rapo_r2_p)
colnames(pgls_ton_rapo_r2_p_df) <- c("adj.r.squared ", "p.value")

##### GET THE ESTIMATES IN TO COLUMNS #######
pgls_sp_ton_rapo_est <- tibble::rownames_to_column(pgls_sp_ton_rapo_est, "Hypothesis")
pgls_sp_ton_rapo_wider <- pgls_sp_ton_rapo_est %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
pgls_sp_ton_rapo_wider_df <- as.data.frame(pgls_sp_ton_rapo_wider[,-1])
pgls_sp_ton_rapo_wider_df

final_pgls_df_sp_ton_rapo <- cbind(pgls_sp_ton_rapo_wider_df, pgls_ton_rapo_r2_p_df)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(final_pgls_df_sp_ton_rapo, "PGLS_rapoport_ton.csv", row.names = FALSE)
######### LEACHE #########
pgls_sp_lea_rapo <- pgls(Range_size ~ Lat, data = comp_lea, lambda = "ML")
pgls_sp_lea_rapo_sum <- summary(pgls_sp_lea_rapo)
pgls_sp_lea_rapo_sum <- as.data.frame(pgls_sp_lea_rapo_sum$coefficients)
pgls_sp_lea_rapo_est <- pgls_sp_lea_rapo_sum[1]

###### R2 AND P VALUE #####
###### GET THE MODEL WITH THE LOWEST AIC #######
pgls_sp_lea_rapo_bm <- summary(pgls_sp_lea_rapo)
pgls_sp_lea_rapo_bm

##### GET THEIR COEFFICIENTS #######
pgls_lea_rapo_r2 <- pgls_sp_lea_rapo_bm$adj.r.squared
pgls_lea_rapo_p <- pgls_sp_lea_rapo_bm$coefficients[2,4]
pgls_lea_rapo_r2_p <- cbind(pgls_lea_rapo_r2, pgls_lea_rapo_p)
pgls_lea_rapo_r2_p_df <- as.data.frame(pgls_lea_rapo_r2_p)
colnames(pgls_lea_rapo_r2_p_df) <- c("adj.r.squared ", "p.value")

##### GET THE ESTIMATES IN TO COLUMNS #######
pgls_sp_lea_rapo_est <- tibble::rownames_to_column(pgls_sp_lea_rapo_est, "Hypothesis")
pgls_sp_lea_rapo_wider <- pgls_sp_lea_rapo_est %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
pgls_sp_lea_rapo_wider_df <- as.data.frame(pgls_sp_lea_rapo_wider[,-1])
pgls_sp_lea_rapo_wider_df

final_pgls_df_sp_lea_rapo <- cbind(pgls_sp_lea_rapo_wider_df, pgls_lea_rapo_r2_p_df)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(final_pgls_df_sp_lea_rapo, "PGLS_rapoport_lea.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################## OLS SITES LEVEL (RAPOPORT) ###############################
#############################################################################################
###### ALL SPECIES ######
ols_sites_rapo_all <- lm(Range_size ~ Lat, sites_vars_all)
ols_sites_rapo_sum_all <- summary(ols_sites_rapo_all)$coefficients
ols_sites_rapo_est_all <- as.data.frame(ols_sites_rapo_sum_all[,1])
##### GET THE ESTIMATES IN TO COLUMNS #######
ols_sites_rapo_est_all <- tibble::rownames_to_column(ols_sites_rapo_est_all, "Hypothesis")
colnames(ols_sites_rapo_est_all) <- c("Hypothesis", "Estimate")
###### R2 AND P VALUE #####
##### GET THEIR COEFFICIENTS #######
ols_sites_rapo_bm_all  <- as.data.frame(broom::glance(ols_sites_rapo_all))
ols_sites_rapo_bm_s_all <- ols_sites_rapo_bm_all[,c(2,5,8,9)]  

ols_sites_rapo_all_wider <- ols_sites_rapo_est_all %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
ols_sites_rapo_all_wider_df <- as.data.frame(ols_sites_rapo_all_wider[,-1])
ols_sites_rapo_all_wider_df

final_ols_df_sites_all_rapo <- cbind(ols_sites_rapo_all_wider_df, ols_sites_rapo_bm_s_all)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(final_ols_df_sites_all_rapo, "OLS_rapoport_All_siteslevel.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
######## TONINI #########
ols_sites_rapo_ton <- lm(Range_size ~ Lat, sites_vars_ton)
ols_sites_rapo_sum_ton <- summary(ols_sites_rapo_ton)$coefficients
ols_sites_rapo_est_ton <- as.data.frame(ols_sites_rapo_sum_ton[,1])
##### GET THE ESTIMATES IN TO COLUMNS #######
ols_sites_rapo_est_ton <- tibble::rownames_to_column(ols_sites_rapo_est_ton, "Hypothesis")
colnames(ols_sites_rapo_est_ton) <- c("Hypothesis", "Estimate")
###### R2 AND P VALUE #####
##### GET THEIR COEFFICIENTS #######
ols_sites_rapo_bm_ton  <- as.data.frame(broom::glance(ols_sites_rapo_ton))
ols_sites_rapo_bm_s_ton <- ols_sites_rapo_bm_ton[,c(2,5,8,9)]  

ols_sites_rapo_ton_wider <- ols_sites_rapo_est_ton %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
ols_sites_rapo_ton_wider_df <- as.data.frame(ols_sites_rapo_ton_wider[,-1])
ols_sites_rapo_ton_wider_df

final_ols_df_sites_ton_rapo <- cbind(ols_sites_rapo_ton_wider_df, ols_sites_rapo_bm_s_ton)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(final_ols_df_sites_ton_rapo, "OLS_rapoport_ton_siteslevel.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
######## LEACHE #########
ols_sites_rapo_lea <- lm(Range_size ~ Lat, sites_vars_lea)
ols_sites_rapo_sum_lea <- summary(ols_sites_rapo_lea)$coefficients
ols_sites_rapo_est_lea <- as.data.frame(ols_sites_rapo_sum_lea[,1])
##### GET THE ESTIMATES IN TO COLUMNS #######
ols_sites_rapo_est_lea <- tibble::rownames_to_column(ols_sites_rapo_est_lea, "Hypothesis")
colnames(ols_sites_rapo_est_lea) <- c("Hypothesis", "Estimate")
###### R2 AND P VALUE #####
##### GET THEIR COEFFICIENTS #######
ols_sites_rapo_bm_lea  <- as.data.frame(broom::glance(ols_sites_rapo_lea))
ols_sites_rapo_bm_s_lea <- ols_sites_rapo_bm_lea[,c(2,5,8,9)]  

ols_sites_rapo_lea_wider <- ols_sites_rapo_est_lea %>%
  pivot_wider(names_from = Hypothesis, values_from = Estimate)
ols_sites_rapo_lea_wider_df <- as.data.frame(ols_sites_rapo_lea_wider[,-1])
ols_sites_rapo_lea_wider_df

final_ols_df_sites_lea_rapo <- cbind(ols_sites_rapo_lea_wider_df, ols_sites_rapo_bm_s_lea)
setwd("G:/7_Doctorado/Cap_1/4_data")
write.csv(final_ols_df_sites_lea_rapo, "OLS_rapoport_lea_siteslevel.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
################################## SARS SITES LEVEL (RAPOPORT) ##############################
#############################################################################################
###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")
##### --------------------------------------------------------------------------------- #####
##### ALL SPECIES #####
###### MATRIX OF SPATIAL COORDINATES #####
sp_all_rapo <- SpatialPoints(data.frame(x = sites_vars_all$Long, y = sites_vars_all$Lat))
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
dir_work_all_rapo <- "G:/7_Doctorado/Cap_1/4_data/SAR/All_rapo"

for (i in 1:length(scheme)){
  ###### Create a directory by scheme #####
  setwd(dir_work_all_rapo)
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
    lm_mod_all_rapo <- lm(paste(path_equation_all_rapo[p, 2]), data = sites_vars_all)
    lm_mod_s_all_rapo <- summary(lm_mod_all_rapo)
    R2_lm_all_rapo <- lm_mod_s_all_rapo[["adj.r.squared"]]
    coef_lm_all_rapo <- lm_mod_all_rapo$coefficients
    p_value_lm_all_rapo <- glance(lm_mod_all_rapo)$p.value
    ##### SAR error models #####
    ##### Minimum distance #####
    error_d_all_rapo <- spatialreg::errorsarlm(paste(path_equation_all_rapo[p, 2]), data = sites_vars_all, 
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
    cor.ols1.res_all_rapo <- correlog(as.numeric(sites_vars_all$Long), as.numeric(sites_vars_all$Lat), z = as.numeric(residuals(lm_mod_all_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
    cor.sar1.res_all_rapo <- correlog(as.numeric(sites_vars_all$Long), as.numeric(sites_vars_all$Lat), z = as.numeric(residuals(error_d_all_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
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
##### SUMMARIZE THE RESULTS #####
results_all_rapo <- data.frame()
for (i in 1:length(scheme)){
  setwd(dir_work_all_rapo)
  ##### Extract results from each scheme directory #####
  dir_scheme <- scheme[i]
  setwd(dir_scheme)
  scheme_results_list_all_rapo <- list.files(path = ".", pattern= '.csv$')
  ##### Join all results #####
  scheme_results_all_rapo <- lapply(scheme_results_list_all_rapo, read.csv)
  scheme_results2_all_rapo <- do.call("rbind", scheme_results_all_rapo)
  results_all_rapo <- rbind(results_all_rapo, scheme_results2_all_rapo)}
##### Back to your working directory #####
setwd(dir_work_all_rapo)
##### Save the model selection results in a csv #####
write.csv(results_all_rapo, file="results_model_selection_all_rapo.csv", row.names = F)

####### TONINI #########
###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")
##### --------------------------------------------------------------------------------- #####
##### TONINI #####
###### MATRIX OF SPATIAL COORDINATES #####
sp_ton_rapo <- SpatialPoints(data.frame(x = sites_vars_ton$Long, y = sites_vars_ton$Lat))
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
dir_work_ton_rapo <- "G:/7_Doctorado/Cap_1/4_data/SAR/Ton_rapo"

for (i in 1:length(scheme)){
  ###### Create a directory by scheme #####
  setwd(dir_work_ton_rapo)
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
    lm_mod_ton_rapo <- lm(paste(path_equation_ton_rapo[p, 2]), data = sites_vars_ton)
    lm_mod_s_ton_rapo <- summary(lm_mod_ton_rapo)
    R2_lm_ton_rapo <- lm_mod_s_ton_rapo[["adj.r.squared"]]
    coef_lm_ton_rapo <- lm_mod_ton_rapo$coefficients
    p_value_lm_ton_rapo <- glance(lm_mod_ton_rapo)$p.value
    ##### SAR error models #####
    ##### Minimum distance #####
    error_d_ton_rapo <- spatialreg::errorsarlm(paste(path_equation_ton_rapo[p, 2]), data = sites_vars_ton, 
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
    cor.ols1.res_ton_rapo <- correlog(as.numeric(sites_vars_ton$Long), as.numeric(sites_vars_ton$Lat), z = as.numeric(residuals(lm_mod_ton_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
    cor.sar1.res_ton_rapo <- correlog(as.numeric(sites_vars_ton$Long), as.numeric(sites_vars_ton$Lat), z = as.numeric(residuals(error_d_ton_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
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
##### SUMMARIZE THE RESULTS #####
results_ton_rapo <- data.frame()
for (i in 1:length(scheme)){
  setwd(dir_work_ton_rapo)
  ##### Extract results from each scheme directory #####
  dir_scheme <- scheme[i]
  setwd(dir_scheme)
  scheme_results_list_ton_rapo <- list.files(path = ".", pattern= '.csv$')
  ##### Join all results #####
  scheme_results_ton_rapo <- lapply(scheme_results_list_ton_rapo, read.csv)
  scheme_results2_ton_rapo <- do.call("rbind", scheme_results_ton_rapo)
  results_ton_rapo <- rbind(results_ton_rapo, scheme_results2_ton_rapo)}
##### Back to your working directory #####
setwd(dir_work_ton_rapo)
##### Save the model selection results in a csv #####
write.csv(results_ton_rapo, file="results_model_selection_ton_rapo.csv", row.names = F)
######## LEACHE ##########
###### Matrix weighting schemes to prove #####
scheme <- c("W", "C", "S")
##### --------------------------------------------------------------------------------- #####
##### LEACHE #####
###### MATRIX OF SPATIAL COORDINATES #####
sp_lea_rapo <- SpatialPoints(data.frame(x = sites_vars_lea$Long, y = sites_vars_lea$Lat))
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
dir_work_lea_rapo <- "G:/7_Doctorado/Cap_1/4_data/SAR/Lea_rapo"

for (i in 1:length(scheme)){
  ###### Create a directory by scheme #####
  setwd(dir_work_lea_rapo)
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
    lm_mod_lea_rapo <- lm(paste(path_equation_lea_rapo[p, 2]), data = sites_vars_lea)
    lm_mod_s_lea_rapo <- summary(lm_mod_lea_rapo)
    R2_lm_lea_rapo <- lm_mod_s_lea_rapo[["adj.r.squared"]]
    coef_lm_lea_rapo <- lm_mod_lea_rapo$coefficients
    p_value_lm_lea_rapo <- glance(lm_mod_lea_rapo)$p.value
    ##### SAR error models #####
    ##### Minimum distance #####
    error_d_lea_rapo <- spatialreg::errorsarlm(paste(path_equation_lea_rapo[p, 2]), data = sites_vars_lea, 
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
    cor.ols1.res_lea_rapo <- correlog(as.numeric(sites_vars_lea$Long), as.numeric(sites_vars_lea$Lat), z = as.numeric(residuals(lm_mod_lea_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
    cor.sar1.res_lea_rapo <- correlog(as.numeric(sites_vars_lea$Long), as.numeric(sites_vars_lea$Lat), z = as.numeric(residuals(error_d_lea_rapo)), na.rm = TRUE, increment = 1, resamp = 1)
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
##### SUMMARIZE THE RESULTS #####
results_lea_rapo <- data.frame()
for (i in 1:length(scheme)){
  setwd(dir_work_lea_rapo)
  ##### Extract results from each scheme directory #####
  dir_scheme <- scheme[i]
  setwd(dir_scheme)
  scheme_results_list_lea_rapo <- list.files(path = ".", pattern= '.csv$')
  ##### Join all results #####
  scheme_results_lea_rapo <- lapply(scheme_results_list_lea_rapo, read.csv)
  scheme_results2_lea_rapo <- do.call("rbind", scheme_results_lea_rapo)
  results_lea_rapo <- rbind(results_lea_rapo, scheme_results2_lea_rapo)}
##### Back to your working directory #####
setwd(dir_work_lea_rapo)
##### Save the model selection results in a csv #####
write.csv(results_lea_rapo, file="results_model_selection_lea_rapo.csv", row.names = F)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
####################################### END OF THE CODE #####################################
#############################################################################################

