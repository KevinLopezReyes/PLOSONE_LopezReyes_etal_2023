#############################################################################################
###################################### PLOTS (RAPOPORT) #####################################
#############################################################################################
library(ggplot2)
library(ggpubr)
library(caper)
library(dplyr)
library(patchwork)
##### ALL SPECIES #####
##### LOAD THE DATA #####
##### SPECIES LEVEL #####
sp_vars_all <- read.csv("G:/7_Doctorado/Cap_1/4_data/all_data_specieslevel.csv")
sp_vars_ton <- read.csv("G:/7_Doctorado/Cap_1/4_data/ton_data_specieslevel.csv")
sp_vars_lea <- read.csv("G:/7_Doctorado/Cap_1/4_data/lea_data_specieslevel.csv")
##### SITES LEVEL #####
sites_vars_all <- read.csv("G:/7_Doctorado/Cap_1/4_data/all_data_siteslevel.csv")
sites_vars_ton <- read.csv("G:/7_Doctorado/Cap_1/4_data/ton_data_siteslevel.csv")
sites_vars_lea <- read.csv("G:/7_Doctorado/Cap_1/4_data/lea_data_siteslevel.csv")

##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### FIGURE 1 HISTOGRAMS ##################################
#############################################################################################
RS_ahull <- read.csv("G:/7_Doctorado/Cap_1/4_data/RS/rs_ahull.csv")
RS_ahull$Species <- gsub(" ", "_", RS_ahull$Species)
RS_ahull_log <- read.csv("G:/7_Doctorado/Cap_1/4_data/RS/rs_ahull_log.csv")
leache <- read.nexus("C:/Users/Kevin/Downloads/Filogenias/leache.nex")
tonini <- read.nexus("C:/Users/Kevin/Downloads/Filogenias/tonini.nex")

####### GET THE DATASETS #######
comp_ton <- comparative.data(tonini, data = RS_ahull, names.col = "Species", vcv = T)
comp_lea <- comparative.data(leache, data = RS_ahull, names.col = "Species", vcv = T)
ton_names <- comp_ton$phy$tip.label
lea_names <- comp_lea$phy$tip.label
RS_ahull_ton <- RS_ahull %>% filter(Species %in% ton_names)
RS_ahull_log_ton <- RS_ahull_log %>% filter(Species %in% ton_names)
RS_ahull_lea <- RS_ahull %>% filter(Species %in% lea_names)
RS_ahull_log_lea <- RS_ahull_log %>% filter(Species %in% lea_names)

options(scipen = 0)
###### PLOT THE HISTOGRAMS ########
all_RS_plot <- ggplot() +  
  geom_histogram(aes(x = RS_ahull$Range_size), 
                 color = "black", fill = "#7bd5ee", size = 0.15) + 
  labs(x = "Range size", y = "Frequency") + ggtitle("All species") + theme_bw() 

ton_RS_plot <- ggplot() +  
  geom_histogram(aes(x = RS_ahull_ton$Range_size), 
                 color = "black", fill = "#7bd5ee", size = 0.15) + 
  labs(x = "Range size", y = "Frequency") + ggtitle("Tonini et al. 2016") + theme_bw() 

lea_RS_plot <- ggplot() +  
  geom_histogram(aes(x = RS_ahull_lea$Range_size), 
                 color = "black", fill = "#7bd5ee", size = 0.15) + 
  labs(x = "Range size", y = "Frequency") + ggtitle("Leache et al. 2016") + theme_bw() 

rs_hist_plot <- ggarrange(all_RS_plot, ton_RS_plot, lea_RS_plot, nrow = 3)

all_RS_log_plot <- ggplot() +  
  geom_histogram(aes(x = RS_ahull_log$Range_size), 
                 color = "black", fill = "#f6e6ac", size = 0.15) + ggtitle("") + 
  labs(x = "log(Range size)", y = "") + theme_bw() 

ton_RS_log_plot <- ggplot() +  
  geom_histogram(aes(x = RS_ahull_log_ton$Range_size), 
                 color = "black", fill = "#f6e6ac", size = 0.15) + ggtitle("") + 
  labs(x = "log(Range size)", y = "") + theme_bw() 

lea_RS_log_plot <- ggplot() +  
  geom_histogram(aes(x = RS_ahull_log_lea$Range_size), 
                 color = "black", fill = "#f6e6ac", size = 0.15) + ggtitle("") + 
  labs(x = "log(Range size)", y = "") + theme_bw() 

log_hist_plot <- ggarrange(all_RS_log_plot, ton_RS_log_plot, lea_RS_log_plot, nrow = 3)

rs_hist_plot + log_hist_plot

setwd("G:/7_Doctorado/Cap_1/5_figures")
ggsave("Fig_1.tiff", width = 6, height = 8, dpi = 1000, bg = "white")
#ggsave("Fig_1.pdf", width = 6, height = 8, dpi = 1000, bg = "white")

##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### FIGURE 2 (RAPOPORT - SPECIES) ########################
#############################################################################################
Rap_sps <- ggplot(sp_vars_all, aes(x = Lat/111000, y = Range_size)) + 
  geom_smooth(method = "lm", se = TRUE, colour = "#297383", size = 0.5, fill = "#8bc5cd") + 
  geom_point(color = "black", fill = "#b4e6ee", stroke = 0.15, shape = 21, size = 0.8, alpha = 1) + 
  labs(x = "Latitudinal Midpoints", y = "log(Range size)")  + 
  ggtitle("All species") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) 
##### TONINI #####
Rap_sps_ton <- ggplot(sp_vars_ton, aes(x = Lat/111000, y = Range_size)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#297383", size = 0.5, fill = "#8bc5cd") + 
  geom_point(color = "black", fill = "#b4e6ee", stroke = 0.15, shape = 21, size = 0.8, alpha = 1) +
  labs(x = "Latitudinal Midpoints", y = "") + 
  ggtitle("Tonini et al. 2016") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6))
##### LEACHE ######
Rap_sps_lea <- ggplot(sp_vars_lea, aes(x = Lat/111000, y = Range_size)) + 
  geom_smooth(method = "lm", se = TRUE, colour = "#297383", size = 0.5, fill = "#8bc5cd") + 
  geom_point(color = "black", fill = "#b4e6ee", stroke = 0.15, shape = 21, size = 0.8, alpha = 1) +
  labs(x = "Latitudinal Midpoints", y = "") + 
  ggtitle("Leache et al. 2016") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) 

ggarrange(Rap_sps, Rap_sps_ton, Rap_sps_lea, ncol = 3, nrow = 1)

setwd("G:/7_Doctorado/Cap_1/5_figures")
ggsave("Fig_2.tiff", width = 6, height = 2, dpi = 1000, bg = "white")
#ggsave("Fig_2.pdf", width = 6, height = 2, dpi = 1000, bg = "white")
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### FIGURE 3 (RAPOPORT - SITES) ##########################
#############################################################################################
##### ALL SPECIES - SITES LEVEL #####
Rap_sites <- ggplot(sites_vars_all, aes(x = Lat/111000, y = Range_size)) + 
  geom_smooth(method = "lm", se = TRUE, colour = "#297383", size = 0.8, fill = "#8bc5cd") + 
  geom_point(color = "black", fill = "#b4e6ee", shape = 21, stroke = 0.06, size = 0.6, alpha = 1) + 
  labs(x = "Absolute latitude", y = "log(Range size)") + 
  ggtitle("All species") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) 
##### TONINI #####
Rap_sites_ton <- ggplot(sites_vars_ton, aes(x = Lat/111000, y = Range_size)) + 
  geom_smooth(method = "lm", se = TRUE, colour = "#297383", size = 0.8, fill = "#8bc5cd") + 
  geom_point(color = "black", fill = "#b4e6ee", shape = 21, stroke = 0.06, size = 0.6, alpha = 1) + 
  labs(x = "Absolute latitude ", y = "") + 
  ggtitle("Tonini et al. 2016") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) 
##### LEACHE #####
Rap_sites_lea <- ggplot(sites_vars_lea, aes(x = Lat/111000, y = Range_size)) + 
  geom_smooth(method = "lm", se = TRUE, colour = "#297383", size = 0.8, fill = "#8bc5cd")+ 
  geom_point(color = "black", fill = "#b4e6ee", shape = 21, stroke = 0.06, size = 0.6, alpha = 1) + 
  labs(x = "Absolute latitude", y = "") + 
  ggtitle("Leache et al. 2016") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6))

ggarrange(Rap_sites, Rap_sites_ton, Rap_sites_lea, ncol = 3, nrow = 1)
setwd("G:/7_Doctorado/Cap_1/5_figures")
ggsave("Fig_3.tiff", width = 6, height = 2, dpi = 1000, bg = "white")
#ggsave("Fig_3.pdf", width = 6, height = 2, dpi = 1000, bg = "white")
##### --------------------------------------------------------------------------------- #####
#############################################################################################
######################## FIGURE 4 (NULL VS OBSERVED COEFFICIENTES - PGLS) ###################
#############################################################################################
####### LOAD THE NULL COEFICIENTES #######
pgls_null_sp_ton <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_pgls_null_sp_ton.csv")
pgls_null_sp_ton_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_pgls_null_sp_ton_rapo.csv")

pgls_null_sp_lea <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_pgls_null_sp_lea.csv")
pgls_null_sp_lea_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_pgls_null_sp_lea_rapo.csv")

##### GET THE INTERVAL CONFIDENCE ########
###### RAPOPORT #######
pgls_quant_rapoport_ton_sp <- quantile(pgls_null_sp_ton_rapo$Lat, c(0, 0.95))
pgls_quant_rapoport_lea_sp <- quantile(pgls_null_sp_lea_rapo$Lat, c(0, 0.95))

####### CEH_min #######
pgls_quant_ceh_ton_sp <- quantile(pgls_null_sp_ton$CEH_min, c(0.05, 1))
pgls_quant_ceh_lea_sp <- quantile(pgls_null_sp_lea$CEH_min, c(0.05, 1))

####### Elevation #######
pgls_quant_elevation_ton_sp <- quantile(pgls_null_sp_ton$Elevation, c(0.05, 1))
pgls_quant_elevation_lea_sp <- quantile(pgls_null_sp_lea$Elevation, c(0.05, 1))

####### CCV #######
pgls_quant_ccv_ton_sp <- quantile(pgls_null_sp_ton$CCV, c(0, 0.95))
pgls_quant_ccv_lea_sp <- quantile(pgls_null_sp_lea$CCV, c(0, 0.95))

####### CVH #######
pgls_quant_cvh_ton_sp <- quantile(pgls_null_sp_ton$CVH, c(0, 0.95))
pgls_quant_cvh_lea_sp <- quantile(pgls_null_sp_lea$CVH, c(0, 0.95))


########### TONINI SPECIES ##############
pgls_null_sp_ton_rapoport <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_ton_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_rapoport_ton_sp[1]) + 
  labs(x = "", y = "Frequency") +
  ggtitle("Tonini et al. 2016 - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.52, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_rapoport_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


pgls_null_sp_ton_ceh <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_ton$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_ceh_ton_sp[2]) + 
  labs(x = "", y = "Frequency") +
  ggtitle("Tonini et al. 2016 - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.33, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_ceh_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


pgls_null_sp_ton_elev <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_ton$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_elevation_ton_sp[2]) + 
  labs(x = "", y = "Frequency") +
  ggtitle("Tonini et al. 2016 - EH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.16, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_elevation_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


pgls_null_sp_ton_ccv <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_ton$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_ccv_ton_sp[1]) + 
  labs(x = "", y = "Frequency") +
  ggtitle("Tonini et al. 2016 - HCSH") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.13, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_ccv_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


pgls_null_sp_ton_cvh <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_ton$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_cvh_ton_sp[1]) + 
  labs(x = "Slope", y = "Frequency") +
  ggtitle("Tonini et al. 2016 - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.17, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_cvh_ton_sp, linetype = "dashed", 
             color = "#3162ac", size = 0.6)

########## Leache et al. 2016 #########
pgls_null_sp_lea_rapoport <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_lea_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15,  boundary = pgls_quant_rapoport_lea_sp[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache et al. 2016 - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_rapoport_lea_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


pgls_null_sp_lea_ceh <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_lea$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_ceh_lea_sp[2]) + 
  labs(x = "", y = "") +
  ggtitle("Leache et al. 2016 - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.36, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_ceh_lea_sp, linetype = "dashed", 
             color = "#3162ac", size = 0.6)


pgls_null_sp_lea_elev <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_lea$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_elevation_lea_sp[2]) + 
  labs(x = "", y = "") +
  ggtitle("Leache et al. 2016 - EH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.14, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_elevation_lea_sp, linetype = "dashed", 
             color = "#3162ac", size = 0.6)


pgls_null_sp_lea_ccv <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_lea$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_ccv_lea_sp[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache et al. 2016 - HCSH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.16, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_ccv_lea_sp, linetype = "dashed", 
             color = "#3162ac", size = 0.6)


pgls_null_sp_lea_cvh <- ggplot() + 
  geom_histogram(aes(x = pgls_null_sp_lea$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = pgls_quant_cvh_lea_sp[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("Leache et al. 2016 - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.08, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = pgls_quant_cvh_lea_sp, linetype = "dashed", 
             color = "#3162ac", size = 0.6)


pgls_sp_null <- ggarrange(pgls_null_sp_ton_rapoport, pgls_null_sp_lea_rapoport, pgls_null_sp_ton_ceh, 
                          pgls_null_sp_lea_ceh, pgls_null_sp_ton_elev, pgls_null_sp_lea_elev, 
                          pgls_null_sp_ton_ccv, pgls_null_sp_lea_ccv, pgls_null_sp_ton_cvh, pgls_null_sp_lea_cvh,
                          ncol = 2, nrow = 5)

pgls_sp_null <- annotate_figure(pgls_sp_null, top = text_grob("PGLS", 
                                                              color = "black", face = "bold", size = 15))

pgls_sp_null
##### --------------------------------------------------------------------------------- #####
#############################################################################################
########################## FIGURE 4 (NULL VS OBSERVED COEFFICIENTES ) SARS ##################
#############################################################################################
####### LOAD THE NULL COEFICIENTES #######
SAR_null_sites_all <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/SAR_All.csv")
SAR_null_sites_all_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/SAR_All_rapo.csv")

SAR_null_sites_ton <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/SAR_ton.csv")
SAR_null_sites_ton_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/SAR_ton_rapo.csv")

SAR_null_sites_lea <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/SAR_lea.csv")
SAR_null_sites_lea_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/SAR_lea_rapo.csv")

##### GET THE INTERVAL CONFIDENCE ########
###### RAPOPORT #######
SAR_quant_rapoport_all_sites <- quantile(SAR_null_sites_all_rapo$Lat, c(0, 0.95))
SAR_quant_rapoport_ton_sites <- quantile(SAR_null_sites_ton_rapo$Lat, c(0, 0.95))
SAR_quant_rapoport_lea_sites <- quantile(SAR_null_sites_lea_rapo$Lat, c(0, 0.95))

####### CEH_min #######
SAR_quant_ceh_all_sites <- quantile(SAR_null_sites_all$CEH_min, c(0.05, 1))
SAR_quant_ceh_ton_sites <- quantile(SAR_null_sites_ton$CEH_min, c(0.05, 1))
SAR_quant_ceh_lea_sites <- quantile(SAR_null_sites_lea$CEH_min, c(0.05, 1))

####### Elevation #######
SAR_quant_elevation_all_sites <- quantile(SAR_null_sites_all$Elevation, c(0.05, 1))
SAR_quant_elevation_ton_sites <- quantile(SAR_null_sites_ton$Elevation, c(0.05, 1))
SAR_quant_elevation_lea_sites <- quantile(SAR_null_sites_lea$Elevation, c(0.05, 1))

####### CCV #######
SAR_quant_ccv_all_sites <- quantile(SAR_null_sites_all$CCV, c(0, 0.95))
SAR_quant_ccv_ton_sites <- quantile(SAR_null_sites_ton$CCV, c(0, 0.95))
SAR_quant_ccv_lea_sites <- quantile(SAR_null_sites_lea$CCV, c(0, 0.95))

####### CVH #######
SAR_quant_cvh_all_sites <- quantile(SAR_null_sites_all$CVH, c(0, 0.95))
SAR_quant_cvh_ton_sites <- quantile(SAR_null_sites_ton$CVH, c(0, 0.95))
SAR_quant_cvh_lea_sites <- quantile(SAR_null_sites_lea$CVH, c(0, 0.95))

######### ALL SITES #############
SAR_null_sites_all_rapoport <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_all_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_rapoport_all_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("All species - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.68, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_rapoport_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_all_ceh <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_all$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_ceh_all_sites[2]) + 
  labs(x = "", y = "") +
  ggtitle("All species - CEH") + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.31, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_ceh_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_all_elev <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_all$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_elevation_all_sites[2]) + 
  labs(x = "", y = "") +
  ggtitle("All species - EH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.14, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_elevation_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_all_ccv <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_all$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_ccv_all_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("All species - HCSH") + theme_bw() + theme(text = element_text(size = 6)) + 
  geom_vline(xintercept = 0.02, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_ccv_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_all_cvh <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_all$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_cvh_all_sites[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("All species - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.08, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_cvh_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

########### TONINI SITES ##############
SAR_null_sites_ton_rapoport <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_ton_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_rapoport_ton_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini et al. 2016 - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.68, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_rapoport_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_ton_ceh <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_ton$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_ceh_ton_sites[2]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini et al. 2016 - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.3, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_ceh_ton_sites, linetype = "dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_ton_elev <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_ton$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_elevation_ton_sites[2]) +
  labs(x = "", y = "") +
  ggtitle("Tonini et al. 2016 - EH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.12, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_elevation_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_ton_ccv <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_ton$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_ccv_ton_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini et al. 2016 - HCSH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.02, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_ccv_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_ton_cvh <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_ton$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_cvh_ton_sites[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("Tonini et al. 2016 - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.07, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_cvh_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

########## LEACHE SITES #########
SAR_null_sites_lea_rapoport <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_lea_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_rapoport_lea_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache et al. 2016 - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.7, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_rapoport_lea_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_lea_ceh <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_lea$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_ceh_lea_sites[2]) + 
  labs(x = "", y = "") +
  ggtitle("Leache et al. 2016 - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.35, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_ceh_lea_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)


SAR_null_sites_lea_elev <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_lea$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_elevation_lea_sites[2]) + 
  labs(x = "", y = "") +
  ggtitle("Leache et al. 2016 - EH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.13, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_elevation_lea_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

SAR_null_sites_lea_ccv <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_lea$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_ccv_lea_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache et al. 2016 - HCSH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.02, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_ccv_lea_sites, linetype = "dashed", 
             color = "#3162ac", size = 0.6)

SAR_null_sites_lea_cvh <- ggplot() + 
  geom_histogram(aes(x = SAR_null_sites_lea$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = SAR_quant_cvh_lea_sites[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("Leache et al. 2016 - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.05, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = SAR_quant_cvh_lea_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

SAR_sites_null <- ggarrange(SAR_null_sites_all_rapoport, SAR_null_sites_ton_rapoport, SAR_null_sites_lea_rapoport, SAR_null_sites_all_ceh, SAR_null_sites_ton_ceh, 
                            SAR_null_sites_lea_ceh, SAR_null_sites_all_elev, SAR_null_sites_ton_elev, SAR_null_sites_lea_elev, SAR_null_sites_all_ccv, 
                            SAR_null_sites_ton_ccv, SAR_null_sites_lea_ccv, SAR_null_sites_all_cvh, SAR_null_sites_ton_cvh, SAR_null_sites_lea_cvh,
                            ncol = 3, nrow = 5)

SAR_sites_null <- annotate_figure(SAR_sites_null, top = text_grob("SARs", 
                                                                  color = "black", face = "bold", size = 15))

setwd("G:/7_Doctorado/Cap_1/5_figures")
pgls_sp_null + SAR_sites_null + plot_layout(widths = c(0.5, 0.75))

ggsave("Fig_4.tiff", width = 10, height = 8, dpi = 1000, bg = "white")
save.image("G:/7_Doctorado/Cap_1/7_Codes/Rdata/7_Plots.RData")
##### --------------------------------------------------------------------------------- #####

#############################################################################################
##### -------------------------------- SUPPLEMENTARY FIGURES -------------------------- #####
#############################################################################################
##### --------------------------------------------------------------------------------- #####
#############################################################################################
##################################### FIGURE S1 (COMPARATION) ###############################
#############################################################################################
load("G:/7_Doctorado/Cap_1/7_Codes/Rdata/5_PAM.RData")
##### RANGE SIZE VS RANGE ZISE FOR ALL METHODS #####
a_c_hull_rs <- ggplot(all_met_RS, aes(RS_ahull, RS_chull)) +
  geom_point(color = "gray30", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Range Size - Alpha Hull", y = "Range Size - Convex Hull") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("Number of species = 103") +  geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) + 
  theme(text = element_text(size = 11), plot.title = element_text(size = 13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) 

a_sdm_rs <- ggplot(all_met_RS, aes(RS_ahull, RS_sdm)) +
  geom_point(color = "gray30", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Range Size - Alpha Hull", y = "Range Size - SDM") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("Number of species = 95") +  geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) + 
  theme(text = element_text(size = 11), plot.title = element_text(size = 13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) 

c_sdm_rs <- ggplot(all_met_RS, aes(RS_chull, RS_sdm)) +
  geom_point(color = "gray30", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Range Size - Convex Hull", y = "Range Size - SDM") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("Number of species = 95") +  geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) + 
  theme(text = element_text(size = 11), plot.title = element_text(size = 13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) 
##### MIDPOINTS VS MIDPOINTS FOR ALL METHODS #####
a_c_hull_mp <- ggplot(all_met_RS, aes(Lat_ahull, Lat_chull)) +
  geom_point(color = "gray30", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Midpoints - Alpha Hull", y = "Midpoints - Convex Hull") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("") +  geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) + 
  theme(text = element_text(size = 11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) 

a_sdm_mp <- ggplot(all_met_RS, aes(Lat_ahull, Lat_sdm)) +
  geom_point(color = "gray30", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Midpoints - Alpha Hull", y = "Midpoints - SDM") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("") +  geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) + 
  theme(text = element_text(size = 11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) 

c_sdm_mp <- ggplot(all_met_RS, aes(Lat_chull, Lat_sdm)) +
  geom_point(color = "gray30", shape = 16, size = 1.5, alpha = 0.6) + 
  stat_cor(method="spearman", cor.coef.name ="rho", label.x.npc="left", label.y.npc="top") +
  labs(x = "Midpoints - Convex Hull", y = "Midpoints - SDM") + 
  geom_smooth(color="red", method ="lm", se = FALSE, size = 0.3) + 
  ggtitle("") +  geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text = element_text(size = 6)) + theme(text = element_text(size = 6)) + 
  theme(text = element_text(size = 11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) 

ggarrange(a_c_hull_rs, a_sdm_rs, c_sdm_rs, a_c_hull_mp, a_sdm_mp, c_sdm_mp, ncol = 3, nrow = 2)
setwd("G:/7_Doctorado/Cap_1/5_figures")
ggsave("Fig_S1.tiff", width = 9, height = 6, dpi = 1000, bg = "white") 

#############################################################################################
############################ FIGURE S2 (NULL VS OBSERVED COEFFICIENTES ) ####################
#############################################################################################
#############################################################################################
#######################################  OLS - SPECIES ######################################  
#############################################################################################
####### LOAD THE NULL COEFICIENTES #######
ols_null_sp_all <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sp_all.csv")
ols_null_sp_all_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sp_all_rapo.csv")

ols_null_sp_ton <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sp_ton.csv")
ols_null_sp_ton_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sp_ton_rapo.csv")

ols_null_sp_lea <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sp_lea.csv")
ols_null_sp_lea_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sp_lea_rapo.csv")

##### GET THE INTERVAL CONFIDENCE ########
###### RAPOPORT #######
ols_quant_rapoport_all_sp <- quantile(ols_null_sp_all_rapo$Lat, c(0, 0.95))
ols_quant_rapoport_ton_sp <- quantile(ols_null_sp_ton_rapo$Lat, c(0, 0.95))
ols_quant_rapoport_lea_sp <- quantile(ols_null_sp_lea_rapo$Lat, c(0, 0.95))

####### CEH_min #######
ols_quant_ceh_all_sp <- quantile(ols_null_sp_all$CEH_min, c(0.05, 1))
ols_quant_ceh_ton_sp <- quantile(ols_null_sp_ton$CEH_min, c(0.05, 1))
ols_quant_ceh_lea_sp <- quantile(ols_null_sp_lea$CEH_min, c(0.05, 1))

####### Elevation #######
ols_quant_elevation_all_sp <- quantile(ols_null_sp_all$Elevation, c(0.05, 1))
ols_quant_elevation_ton_sp <- quantile(ols_null_sp_ton$Elevation, c(0.05, 1))
ols_quant_elevation_lea_sp <- quantile(ols_null_sp_lea$Elevation, c(0.05, 1))

####### CCV #######
ols_quant_ccv_all_sp <- quantile(ols_null_sp_all$CCV, c(0, 0.95))
ols_quant_ccv_ton_sp <- quantile(ols_null_sp_ton$CCV, c(0, 0.95))
ols_quant_ccv_lea_sp <- quantile(ols_null_sp_lea$CCV, c(0, 0.95))

####### CVH #######
ols_quant_cvh_all_sp <- quantile(ols_null_sp_all$CVH, c(0, 0.95))
ols_quant_cvh_ton_sp <- quantile(ols_null_sp_ton$CVH, c(0, 0.95))
ols_quant_cvh_lea_sp <- quantile(ols_null_sp_lea$CVH, c(0, 0.95))

######### ALL SPECIES #############
ols_null_sp_all_rapoport <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_all_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_rapoport_all_sp[1]) + 
  labs(x = "", y = "Frequency") +
  ggtitle("All species - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.53, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_rapoport_all_sp, linetype = "dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_all_ceh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_all$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ceh_all_sp[2]) + 
  labs(x = "", y = "Frequency") +
  ggtitle("All species - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.45, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ceh_all_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_all_elev <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_all$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_elevation_all_sp[2]) + 
  labs(x = "", y = "Frequency") +
  ggtitle("All species - Elevation") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.18, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_elevation_all_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_all_ccv <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_all$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ccv_all_sp[1]) + 
  labs(x = "", y = "Frequency") +
  ggtitle("All species - CCV") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.15, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ccv_all_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_all_cvh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_all$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_cvh_all_sp[1]) + 
  labs(x = "Slope", y = "Frequency") +
  ggtitle("All species - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.04, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_cvh_all_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)

########### TONINI ##############
ols_null_sp_ton_rapoport <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_ton_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_rapoport_ton_sp[1]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini species - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.52, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_rapoport_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_ton_ceh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_ton$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ceh_ton_sp[2]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini species - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.33, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ceh_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_ton_elev <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_ton$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_elevation_ton_sp[2]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini species - Elevation") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.16, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_elevation_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_ton_ccv <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_ton$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ccv_ton_sp[1]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini species - CCV") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.13, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ccv_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_ton_cvh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_ton$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_cvh_ton_sp[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("Tonini species - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.17, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_cvh_ton_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)

########## LEACHE SPECIES #########
ols_null_sp_lea_rapoport <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_lea_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_rapoport_lea_sp[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache species - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.50, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_rapoport_lea_sp, linetype = "dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_lea_ceh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_lea$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ceh_lea_sp[2]) + 
  labs(x = "", y = "") +
  ggtitle("Leache species - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.36, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ceh_lea_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_lea_elev <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_lea$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_elevation_lea_sp[2]) + 
  labs(x = "", y = "") +
  ggtitle("Leache species - Elevation") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.14, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_elevation_lea_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_lea_ccv <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_lea$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ccv_lea_sp[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache species - CCV") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.16, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ccv_lea_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)


ols_null_sp_lea_cvh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sp_lea$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_cvh_lea_sp[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("Leache species - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.08, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_cvh_lea_sp, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_sp_null <- ggarrange(ols_null_sp_all_rapoport, ols_null_sp_ton_rapoport, ols_null_sp_lea_rapoport, ols_null_sp_all_ceh, ols_null_sp_ton_ceh, 
                         ols_null_sp_lea_ceh, ols_null_sp_all_elev, ols_null_sp_ton_elev, ols_null_sp_lea_elev, ols_null_sp_all_ccv, 
                         ols_null_sp_ton_ccv, ols_null_sp_lea_ccv, ols_null_sp_all_cvh, ols_null_sp_ton_cvh, ols_null_sp_lea_cvh,
                         ncol = 3, nrow = 5)

ols_sp_null <- annotate_figure(ols_sp_null, top = text_grob("OLS species level", 
                                                            color = "black", face = "bold", size = 15))
ols_sp_null 
##### --------------------------------------------------------------------------------- #####
#############################################################################################
#######################################  OLS - SITES ######################################  
#############################################################################################
####### LOAD THE NULL COEFICIENTES #######
ols_null_sites_all <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sites_all.csv")
ols_null_sites_all_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sites_all_rapo.csv")

ols_null_sites_ton <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sites_ton.csv")
ols_null_sites_ton_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sites_ton_rapo.csv")

ols_null_sites_lea <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sites_lea.csv")
ols_null_sites_lea_rapo <- read.csv("G:/7_Doctorado/Cap_1/4_data/Null_models/final_ols_null_sites_lea_rapo.csv")

##### GET THE INTERVAL CONFIDENCE ########
###### RAPOPORT #######
ols_quant_rapoport_all_sites <- quantile(ols_null_sites_all_rapo$Lat, c(0, 0.95))
ols_quant_rapoport_ton_sites <- quantile(ols_null_sites_ton_rapo$Lat, c(0, 0.95))
ols_quant_rapoport_lea_sites <- quantile(ols_null_sites_lea_rapo$Lat, c(0, 0.95))

####### CEH_min #######
ols_quant_ceh_all_sites <- quantile(ols_null_sites_all$CEH_min, c(0.05, 1))
ols_quant_ceh_ton_sites <- quantile(ols_null_sites_ton$CEH_min, c(0.05, 1))
ols_quant_ceh_lea_sites <- quantile(ols_null_sites_lea$CEH_min, c(0.05, 1))

####### Elevation #######
ols_quant_elevation_all_sites <- quantile(ols_null_sites_all$Elevation, c(0.05, 1))
ols_quant_elevation_ton_sites <- quantile(ols_null_sites_ton$Elevation, c(0.05, 1))
ols_quant_elevation_lea_sites <- quantile(ols_null_sites_lea$Elevation, c(0.05, 1))

####### CCV #######
ols_quant_ccv_all_sites <- quantile(ols_null_sites_all$CCV, c(0, 0.95))
ols_quant_ccv_ton_sites <- quantile(ols_null_sites_ton$CCV, c(0, 0.95))
ols_quant_ccv_lea_sites <- quantile(ols_null_sites_lea$CCV, c(0, 0.95))

####### CVH #######
ols_quant_cvh_all_sites <- quantile(ols_null_sites_all$CVH, c(0, 0.95))
ols_quant_cvh_ton_sites <- quantile(ols_null_sites_ton$CVH, c(0, 0.95))
ols_quant_cvh_lea_sites <- quantile(ols_null_sites_lea$CVH, c(0, 0.95))

######### ALL SITES #############
ols_null_sites_all_rapoport <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_all_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_rapoport_all_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("All species - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.72, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_rapoport_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_all_ceh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_all$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ceh_all_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("All species - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.84, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ceh_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_all_elev <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_all$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_elevation_all_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("All species - Elevation") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.43, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_elevation_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_all_ccv <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_all$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ccv_all_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("All species - CCV") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.21, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ccv_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_all_cvh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_all$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_cvh_all_sites[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("All species - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.08, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_cvh_all_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

########### TONINI SITES ##############
ols_null_sites_ton_rapoport <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_ton_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_rapoport_ton_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini species - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.68, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_rapoport_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_ton_ceh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_ton$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ceh_ton_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini species - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.85, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ceh_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_ton_elev <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_ton$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_elevation_ton_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini species - Elevation") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.42, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_elevation_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_ton_ccv <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_ton$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ccv_ton_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Tonini species - CCV") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.22, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ccv_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_ton_cvh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_ton$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_cvh_ton_sites[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("Tonini species - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.089, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_cvh_ton_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

########## LEACHE SITES #########
ols_null_sites_lea_rapoport <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_lea_rapo$Lat), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_rapoport_lea_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache species - Rapoport's rule") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.71, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_rapoport_lea_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_lea_ceh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_lea$CEH_min), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ceh_lea_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache species - CEH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.80, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ceh_lea_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_lea_elev <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_lea$Elevation), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_elevation_lea_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache species - Elevation") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.42, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_elevation_lea_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_lea_ccv <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_lea$CCV), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_ccv_lea_sites[1]) + 
  labs(x = "", y = "") +
  ggtitle("Leache species - CCV") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 0.22, linetype="dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_ccv_lea_sites, linetype="dashed", 
             color = "#3162ac", size = 0.6)

ols_null_sites_lea_cvh <- ggplot() + 
  geom_histogram(aes(x = ols_null_sites_lea$CVH), 
                 color = "black", fill = "#d5dede", size = 0.15, boundary = ols_quant_cvh_lea_sites[1]) + 
  labs(x = "Slope", y = "") +
  ggtitle("Leache species - CVH") + theme_bw() + theme(text = element_text(size = 6)) +
  geom_vline(xintercept = -0.05, linetype = "dashed", 
             color = "#ac2020", size = 0.6) + 
  geom_vline(xintercept = ols_quant_cvh_lea_sites, linetype = "dashed", 
             color = "#3162ac", size = 0.6)

ols_sites_null <- ggarrange(ols_null_sites_all_rapoport, ols_null_sites_ton_rapoport, ols_null_sites_lea_rapoport, ols_null_sites_all_ceh, ols_null_sites_ton_ceh, 
                            ols_null_sites_lea_ceh, ols_null_sites_all_elev, ols_null_sites_ton_elev, ols_null_sites_lea_elev, ols_null_sites_all_ccv, 
                            ols_null_sites_ton_ccv, ols_null_sites_lea_ccv, ols_null_sites_all_cvh, ols_null_sites_ton_cvh, ols_null_sites_lea_cvh,
                            ncol = 3, nrow = 5)

ols_sites_null <- annotate_figure(ols_sites_null, top = text_grob("OLS sites level", 
                                                                  color = "black", face = "bold", size = 15))

library(patchwork)
setwd("G:/7_Doctorado/Cap_1/5_figures")
ols_sp_null + ols_sites_null
ggsave("Fig_S2.tiff", width = 12, height = 8, dpi = 1000, bg = "white")
##### --------------------------------------------------------------------------------- #####
#############################################################################################
###################################### END OF THE CODE ######################################
#############################################################################################
##### --------------------------------------------------------------------------------- #####

