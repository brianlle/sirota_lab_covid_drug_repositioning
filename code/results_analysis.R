#####################################################################
# Comp Immuno Mini-course Workshop: Computational Drug Repositioning
# Case study: COVID-19
#
# Component #3: Inspect the results! What drugs come up? How do the signatures look in comparison
# to the drug profiles?
# Adapted from code written by Brian Le, Bin Chen, and Marina Sirota
# For more details, see: Le et al.: https://insight.jci.org/articles/view/133761
# For original pipelien code, see: Chen et al.: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5447464/
#####################################################################

library(pheatmap)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(AnnotationDbi)
library(reshape2)
library(pheatmap)

dataset <- "ALV" # one of "ALV", "EXP", "BALF"

load(paste0("results/cmap_predictions_", dataset, ".RData"))  #pipeline results
load('data/cmap_signatures.RData')   #cmap_signatures
cmap_experiments <- read.csv("data/cmap_drug_experiments_new.csv", stringsAsFactors =  F) #cmap profiles metadata
valid_instances <- read.csv("data/cmap_valid_instances.csv", stringsAsFactors = F)

drug_preds <- results[[1]]
dz_sig <- results[[2]]

#keep valid (concordant) profiles; keep drugs listed in DrugBank
cmap_experiments_valid <- merge(cmap_experiments, valid_instances, by="id")
cmap_experiments_valid <- subset(cmap_experiments_valid, valid == 1 & DrugBank.ID != "NULL")

drug_instances_all <- merge(drug_preds, cmap_experiments_valid, by.x="exp_id", by.y="id")

#We can inspect the distribution of reversal scores
hist(drug_instances_all$cmap_score, breaks = 10)

#Apply thresholds for significant hits: here, we apply FDR < 0.05 and keep only the reversed profiles (cmap_score < 0)
drug_instances <- subset(drug_instances_all, q < 0.05 & cmap_score < 0)

#Since drugs in cmap have been tested multiple times, we keep the most negative score to be as inclusive as possible
#Alternatively, you could aggregate the scores in a different manner (e.g. averages, or based on metadata)
drug_instances <- drug_instances %>% 
  group_by(name) %>% 
  dplyr::slice(which.min(cmap_score))

drug_instances <- drug_instances[order(drug_instances$cmap_score), ]
drug_instances_id <- c(drug_instances$exp_id) + 1 #the first column is the gene id
#get candidate drugs
drug_signatures <- cmap_signatures[,c(1, drug_instances_id)] #the first column is the gene id

write.csv(drug_instances, file = paste0("results/", dataset, "_hits.csv"))

drug_dz_signature <- merge(dz_sig[, c("GeneID", "log2FoldChange")], drug_signatures, by.x = "GeneID", by.y="V1")
colnames(drug_dz_signature)[2] <- "value"
drug_dz_signature <- drug_dz_signature[order(drug_dz_signature$value),]

#Convert disease and drug values to ranks from 1:numgenes
#Higher rank corresponds to more overexpressed, so we need to reverse order of disease sig
drug_dz_signature[,2] <- -drug_dz_signature[,2] 
for (i in 2:ncol(drug_dz_signature)){
  drug_dz_signature[,i] <- rank(drug_dz_signature[,i] )
}
drug_dz_signature <- drug_dz_signature[order(drug_dz_signature[,2]),] #order by disease expression

gene_ids <- drug_dz_signature[,1]
drug_dz_signature <- drug_dz_signature[, -1]

drug_names <- sapply(2:ncol(drug_dz_signature), function(id){
  #need to subtract 1 as in cmap_signatures, V1 is gene id.
  new_id <- strtoi(paste(unlist(strsplit(as.character(colnames(drug_dz_signature)[id]),""))[-1], collapse="")) - 1 
  cmap_experiments_valid$name[cmap_experiments_valid$id == new_id]
})
colnames(drug_dz_signature)[-1] <- drug_names

write.csv(drug_dz_signature, paste0("results/", dataset, "drug_dz_signature_all_hits.csv"))

#FIGURE: all hits from cmap, using a red/blue color scheme
pdf(paste0("results/", dataset, "_heatmap_cmap_hits.pdf"), width = 12, height = 15)
layout(matrix(1))
par(mar=c(6, 4, 1, 0.5))
colPal <- redblue(100)
image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels= F)
axis(2,  at=seq(0,1,length.out= length (gene_ids) ), labels= F)
text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.015),
     labels = c(dataset, drug_names), srt = 45, pos=2, offset=-0.2, xpd = TRUE, cex=0.7, 
     col = "black")
dev.off()


