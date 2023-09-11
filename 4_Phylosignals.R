#############################################################################################
################################## PHYLOGENETIC SIGNALS #####################################
#############################################################################################
library(phytools)
#############################################################################################
library(caper)
leache <- read.nexus("C:/Users/Kevin/Downloads/Filogenias/leache.nex")
tonini <- read.nexus("C:/Users/Kevin/Downloads/Filogenias/tonini.nex")
##### TONININ #####
sp_vars_all <- read.csv("G:/7_Doctorado/Cap_1/4_data/all_data_specieslevel.csv")
comp_ton <- comparative.data(tonini, data = sp_vars_all, names.col = "Species", vcv = T)
##### K #####
ton_k_lat <- phylosig(comp_ton$phy, comp_ton$data$Lat, method = "K", test=TRUE, nsim = 1000)
ton_k_rs <- phylosig(comp_ton$phy, comp_ton$data$Range_size, method = "K", test=TRUE, nsim = 1000)
ton_k_CCV <- phylosig(comp_ton$phy, comp_ton$data$CCV, method="K", test = TRUE, nsim = 1000)
ton_k_Elev <- phylosig(comp_ton$phy, comp_ton$data$Elevation, method ="K", test = TRUE, nsim=1000)
ton_k_CEH_min <- phylosig(comp_ton$phy, comp_ton$data$CEH_min, method ="K", test = TRUE, nsim=1000)
ton_k_CEH_max <- phylosig(comp_ton$phy, comp_ton$data$CEH_max, method ="K", test = TRUE, nsim=1000)
ton_k_CVH <- phylosig(comp_ton$phy, comp_ton$data$CVH, method = "K", test = TRUE, nsim=1000)
##### LAMBDA #####
ton_lambda_lat <- phylosig(comp_ton$phy, comp_ton$data$Lat, method = "lambda", test=TRUE, nsim = 1000)
ton_lambda_rs <- phylosig(comp_ton$phy, comp_ton$data$Range_size, method = "lambda", test=TRUE, nsim = 1000)
ton_lambda_CCV <- phylosig(comp_ton$phy, comp_ton$data$CCV, method="lambda", test = TRUE, nsim = 1000)
ton_lambda_Elev <- phylosig(comp_ton$phy, comp_ton$data$Elevation, method ="lambda", test = TRUE, nsim=1000)
ton_lambda_CEH_min <- phylosig(comp_ton$phy, comp_ton$data$CEH_min, method ="lambda", test = TRUE, nsim=1000)
ton_lambda_CEH_max <- phylosig(comp_ton$phy, comp_ton$data$CEH_max, method ="lambda", test = TRUE, nsim=1000)
ton_lambda_CVH <- phylosig(comp_ton$phy, comp_ton$data$CVH, method = "lambda", test = TRUE, nsim=1000)
##### --------------------------------------------------------------------------------- #####
##### LEACHE #####
comp_lea <- comparative.data(leache, data = sp_vars_all, names.col = "Species", vcv = T)
##### K #####
lea_k_lat <- phylosig(comp_lea$phy, comp_lea$data$Lat, method = "K", test=TRUE, nsim = 1000)
lea_k_rs <- phylosig(comp_lea$phy, comp_lea$data$Range_size, method = "K", test=TRUE, nsim = 1000)
lea_k_CCV <- phylosig(comp_lea$phy, comp_lea$data$CCV, method="K", test = TRUE, nsim = 1000)
lea_k_Elev <- phylosig(comp_lea$phy, comp_lea$data$Elevation, method ="K", test = TRUE, nsim=1000)
lea_k_CEH_min <- phylosig(comp_lea$phy, comp_lea$data$CEH_min, method ="K", test = TRUE, nsim=1000)
lea_k_CEH_max <- phylosig(comp_lea$phy, comp_lea$data$CEH_max, method ="K", test = TRUE, nsim=1000)
lea_k_CVH <- phylosig(comp_lea$phy, comp_lea$data$CVH, method = "K", test = TRUE, nsim=1000)
##### LAMBDA #####
lea_lambda_lat <- phylosig(comp_lea$phy, comp_lea$data$Lat, method = "lambda", test=TRUE, nsim = 1000)
lea_lambda_rs <- phylosig(comp_lea$phy, comp_lea$data$Range_size, method = "lambda", test=TRUE, nsim = 1000)
lea_lambda_CCV <- phylosig(comp_lea$phy, comp_lea$data$CCV, method="lambda", test = TRUE, nsim = 1000)
lea_lambda_Elev <- phylosig(comp_lea$phy, comp_lea$data$Elevation, method ="lambda", test = TRUE, nsim=1000)
lea_lambda_CEH_min <- phylosig(comp_lea$phy, comp_lea$data$CEH_min, method ="lambda", test = TRUE, nsim=1000)
lea_lambda_CEH_max <- phylosig(comp_lea$phy, comp_lea$data$CEH_max, method ="lambda", test = TRUE, nsim=1000)
lea_lambda_CVH <- phylosig(comp_lea$phy, comp_lea$data$CVH, method = "lambda", test = TRUE, nsim=1000)
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### END OF THE CODE ######################################
#############################################################################################
##### --------------------------------------------------------------------------------- #####



