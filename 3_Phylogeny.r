
library(ape)
library(caper)
library(stringr)
library(gdata)
setwd("C:/Users/Kevin/Downloads/Filogenias")

### load the phylogenetic trees ###
tonini <- read.nexus("tonini.nex")
leache <- read.nexus("leache.nex") 

### Prune Leache et al. 2016 tree ###
## create a vector containing the names that i will prune ###
names <- leache$tip.label
### sub species names ###
names <- gsub("_", " ", names)

words_list <- list()
for(i in 1:length(names)){
name <- stringi::stri_count_words(names[i])
if(name > 2){
words_list[i] <- names[i]}}
words <- unlist(words_list)

w_list <- list()
for(i in 1:length(names)){
name <- grepl("Sceloporus", names[i], fixed = TRUE)
if(name == "FALSE"){
w_list[i] <- names[i]}}
w <- unlist(w_list)

fin_list <- combine(w, words)
fin_list <- fin_list[,1]
fin_list <- gsub(" ", "_", fin_list)

#Export
leache_prune <- drop.tip(leache, c(fin_list,"Sceloporus_lineolateralis", 
"Sceloporus_vandenburgianus"))
leache_prune <- ladderize(leache_prune, right = F)

leache_prune$tip.label

setwd("C:/Users/Kevin/Downloads/Filogenias")
ape::write.nexus(leache_prune, file = "Leache_pruned.nex")

### Obtain times from each species per phylogeny ###
ages_leache <- setNames(leache_prune$edge.length[sapply(1:length(leache_prune$tip.label), 
function(x,y) which (y==x), 
y = leache_prune$edge[,2])], leache_prune$tip.label)

ages_tonini <- setNames(tonini$edge.length[sapply(1:length(tonini$tip.label), 
function(x,y) which (y==x), 
y = tonini$edge[,2])], tonini$tip.label)

as.data.frame(ages_leache)
as.data.frame(ages_tonini)
write.csv(ages_leache, "ages_leache.csv")
write.csv(ages_tonini, "ages_tonini.csv")

ages_leache <- read.csv("C:/Users/Kevin/Downloads/Filogenias/ages_leache.csv")
ages_tonini <- read.csv("C:/Users/Kevin/Downloads/Filogenias/ages_tonini.csv")

colnames(ages_leache) <- c("species", "ages")
colnames(ages_tonini) <- c("species", "ages")

ages_final <- merge(ages_leache, ages_tonini, by="species", all = TRUE)
colnames(ages_final) <- c("species", "Leache", "Tonini")

write.csv(ages_final, "ages_final.csv", row.names = FALSE)
ages_final <- read.csv("C:/Users/Kevin/Downloads/Filogenias/ages_final.csv")

library(dplyr)
ages_final <- ages_final %>% filter_at(vars(Tonini,Leache), all_vars(!is.na(.)))
#ages_final$Leache <- (ages_final$Leache - min(ages_final$Leache)) / (max(ages_final$Leache) - min(ages_final$Leache))

#ages_final$Tonini <- (ages_final$Tonini - min(ages_final$Tonini)) / (max(ages_final$Tonini) - min(ages_final$Tonini))

plot(ages_final$Leache, ages_final$Tonini, main="Species diversification time",
xlab="Leache et al. 2016", 
ylab="Tonini et al. 2016")
abline(lm(ages_final$Tonini ~ ages_final$Leache))
#abline(coef = c(0, 1), col = "red")
shapiro.test(ages_final$Leache)
shapiro.test(ages_final$Tonini)
cor.test(ages_final$Leache, ages_final$Tonini, method = "spearman")
#
library(reshape2)
ages_f2 <- melt(ages_final)
#install.packages("lawstat")
library(lawstat)
levene.test(ages_f2$value, ages_f2$variable, location = "median")
wilcox.test(value ~ variable, data = ages_f2, paired = FALSE)

##install.packages("MKinfer")
library(MKinfer)
boot.t.test(value ~ variable, data = ages_f2, reps = 10000)
boxplot(ages_final$Leache, ages_final$Tonini)
hist(ages_final$Leache)
hist(ages_final$Tonini)

library(ggplot2)
line <- ggplot(ages_final, aes(x = Leache, y = Tonini)) + geom_point(shape=1, alpha = 1) + 
geom_smooth(method = "lm") + theme_bw()
boxplot <- ggplot(ages_f2, aes(x = variable, y = value)) + geom_boxplot() + theme_bw()
library(ggpubr)
arg <- ggarrange(line, boxplot, ncol = 2, nrow = 1)

save.image("G:/7_Doctorado/Cap_1/7_Codes/Rdata/3_Phylogeny.RData")


