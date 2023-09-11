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
##### GET THE VALUES TO A DATAFRAMES #####
##### TONINI - K #####
ton_k_lat_df <- data.frame(Trait = "Latitude", K = ton_k_lat$K, pvalue = ton_k_lat$P)
ton_k_rs_df <- data.frame(Trait = "Range size", K = ton_k_rs$K, pvalue = ton_k_rs$P)
ton_k_CCV_df <- data.frame(Trait = "Climate change velocity", K = ton_k_CCV$K, pvalue = ton_k_CCV$P)
ton_k_Elev_df <- data.frame(Trait = "Elevation", K = ton_k_Elev$K, pvalue = ton_k_Elev$P)
ton_k_CEH_min_df <- data.frame(Trait = "Temperature limit (lower)", K = ton_k_CEH_min$K, pvalue = ton_k_CEH_min$P)
ton_k_CEH_max_df <- data.frame(Trait = "Temperature limit (upper)", K = ton_k_CEH_max$K, pvalue = ton_k_CEH_max$P)
ton_k_CVH_df <- data.frame(Trait = "Temperature range", K = ton_k_CVH$K, pvalue = ton_k_CVH$P)
all_ton_k_df <- rbind(ton_k_lat_df, ton_k_rs_df, ton_k_CCV_df, ton_k_Elev_df, ton_k_CEH_min_df, ton_k_CEH_max_df, ton_k_CVH_df)
##### TONINI - LAMBDA #####
ton_lambda_lat_df <- data.frame(Trait = "Latitude", lambda = ton_lambda_lat$lambda, pvalue = ton_lambda_lat$P)
ton_lambda_rs_df <- data.frame(Trait = "Range size", lambda = ton_lambda_rs$lambda, pvalue = ton_lambda_rs$P)
ton_lambda_CCV_df <- data.frame(Trait = "Climate change velocity", lambda = ton_lambda_CCV$lambda, pvalue = ton_lambda_CCV$P)
ton_lambda_Elev_df <- data.frame(Trait = "Elevation", lambda = ton_lambda_Elev$lambda, pvalue = ton_lambda_Elev$P)
ton_lambda_CEH_min_df <- data.frame(Trait = "Temperature limit (lower)", lambda = ton_lambda_CEH_min$lambda, pvalue = ton_lambda_CEH_min$P)
ton_lambda_CEH_max_df <- data.frame(Trait = "Temperature limit (upper)", lambda = ton_lambda_CEH_max$lambda, pvalue = ton_lambda_CEH_max$P)
ton_lambda_CVH_df <- data.frame(Trait = "Temperature range", lambda = ton_lambda_CVH$lambda, pvalue = ton_lambda_CVH$P)
all_ton_lambda_df <- rbind(ton_lambda_lat_df, ton_lambda_rs_df, ton_lambda_CCV_df, ton_lambda_Elev_df, ton_lambda_CEH_min_df, ton_lambda_CEH_max_df, ton_lambda_CVH_df)
##### LEACHE - K #####
lea_k_lat_df <- data.frame(Trait = "Latitude", K = lea_k_lat$K, pvalue = lea_k_lat$P)
lea_k_rs_df <- data.frame(Trait = "Range size", K = lea_k_rs$K, pvalue = lea_k_rs$P)
lea_k_CCV_df <- data.frame(Trait = "Climate change velocity", K = lea_k_CCV$K, pvalue = lea_k_CCV$P)
lea_k_Elev_df <- data.frame(Trait = "Elevation", K = lea_k_Elev$K, pvalue = lea_k_Elev$P)
lea_k_CEH_min_df <- data.frame(Trait = "Temperature limit (lower)", K = lea_k_CEH_min$K, pvalue = lea_k_CEH_min$P)
lea_k_CEH_max_df <- data.frame(Trait = "Temperature limit (upper)", K = lea_k_CEH_max$K, pvalue = lea_k_CEH_max$P)
lea_k_CVH_df <- data.frame(Trait = "Temperature range", K = lea_k_CVH$K, pvalue = lea_k_CVH$P)
all_lea_k_df <- rbind(lea_k_lat_df, lea_k_rs_df, lea_k_CCV_df, lea_k_Elev_df, lea_k_CEH_min_df, lea_k_CEH_max_df, lea_k_CVH_df)
##### LEACHE - LAMBDA #####
lea_lambda_lat_df <- data.frame(Trait = "Latitude", lambda = lea_lambda_lat$lambda, pvalue = lea_lambda_lat$P)
lea_lambda_rs_df <- data.frame(Trait = "Range size", lambda = lea_lambda_rs$lambda, pvalue = lea_lambda_rs$P)
lea_lambda_CCV_df <- data.frame(Trait = "Climate change velocity", lambda = lea_lambda_CCV$lambda, pvalue = lea_lambda_CCV$P)
lea_lambda_Elev_df <- data.frame(Trait = "Elevation", lambda = lea_lambda_Elev$lambda, pvalue = lea_lambda_Elev$P)
lea_lambda_CEH_min_df <- data.frame(Trait = "Temperature limit (lower)", lambda = lea_lambda_CEH_min$lambda, pvalue = lea_lambda_CEH_min$P)
lea_lambda_CEH_max_df <- data.frame(Trait = "Temperature limit (upper)", lambda = lea_lambda_CEH_max$lambda, pvalue = lea_lambda_CEH_max$P)
lea_lambda_CVH_df <- data.frame(Trait = "Temperature range", lambda = lea_lambda_CVH$lambda, pvalue = lea_lambda_CVH$P)
all_lea_lambda_df <- rbind(lea_lambda_lat_df, lea_lambda_rs_df, lea_lambda_CCV_df, lea_lambda_Elev_df, lea_lambda_CEH_min_df, lea_lambda_CEH_max_df, lea_lambda_CVH_df)
##### MERGE ALL OF THEM IN A SINGLE DATAFRAME #####
all_ton_k_lambda <- cbind(all_ton_k_df, all_ton_lambda_df[c(2,3)])
all_ton_k_lambda$Phylogeny <- "Tonini"
all_lea_k_lambda <- cbind(all_lea_k_df, all_lea_lambda_df[c(2,3)])
all_lea_k_lambda$Phylogeny <- "Leache"
##### ALL ####
all_phylosig <- rbind(all_ton_k_lambda, all_lea_k_lambda)
setwd("G:/7_Doctorado/Cap_1/6_tables")
options(scipen=999)
write.csv(all_phylosig, "phylosig_traits.csv", row.names = FALSE)
##### --------------------------------------------------------------------------------- #####

save.image("G:/7_Doctorado/Cap_1/7_Codes/Rdata/4_Phylosignals.RData")


