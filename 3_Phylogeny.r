
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


save.image("G:/7_Doctorado/Cap_1/7_Codes/Rdata/3_Phylogeny.RData")


