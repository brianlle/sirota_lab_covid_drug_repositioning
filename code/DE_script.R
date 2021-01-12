# ==================================================================================== #
# This for differential expression analyses using COVID data from 
# https://www.tandfonline.com/doi/full/10.1080/22221751.2020.1747363
# ==================================================================================== #
#!/usr/bin/r

## get wording directory
getwd()

## create a directory for Mac and PC
dir.create("~/Google Drive/COVID") #  working directory (mac)
dir.create("C:/Users/tim/Google Drive/COVID")  # working directory (PC)

## move inside the created directory
setwd("../COVID")

## set working directory for Mac and PC
setwd("~/Google Drive/COVID")      #  working directory (mac)
setwd("C:/Users/tim/Google Drive/COVID")    #  working directory (PC)


## installlibraries
install.packages("DESeq2") ## for diff exp analysis 
install.packages("readr") ## fast and friendly way to read rectangular data (like csv, tsv, and fwf). 
## install.packages("ReportingTools")
##install.packages("AnnotationDbi")
##install.packages("ensembldb")
## install.packages("EnsDb.Hsapiens.v79")
install.packages("devtools") # devtools simplifies many common tasks
devtools::install_github("stephenturner/annotables") # gene and trascript annotation package
#BiocManager::install("DEGreport")
BiocManager::install('EnhancedVolcano') ## volcano plots
install.packages("pheatmap") ## A package for drawing pretty heatmaps in R
install.packages("fgsea") ## The package implements an algorithm for fast gene set enrichment analysis.

## load the libraries
library("DESeq2")
library("readr")
library("AnnotationDbi")
library(ensembldb)
library("annotables")
library(EnhancedVolcano)
library(fgsea)


##  Diff exoresiion on PBMC 
PBMC_counts <- read.csv("PBMC_counts.csv", head = TRUE)
rownames(PBMC_counts) <- PBMC_counts$X
PBMC_counts$X <- NULL
## write.csv(PBMC_counts,"PBMC_counts.csv")

dim(PBMC_counts)       #Shows the dimensions of the table, i.e. rows and columns.
colnames(PBMC_counts)    #Shows the names of columns, if any.
summary(PBMC_counts)    #Outputs summary statistics of data. Ex- minimum, max, median, etc.
head(PBMC_counts) #Displays first 6 rows of data.
tail(PBMC_counts)      #Displays last 6 rows of data. 

## load experimental design
PBMC_design <- read.table("PBMC_design.txt", header = T)
rownames(PBMC_design)<- PBMC_design$sampleID

dim(PBMC_design)       #Shows the dimensions of the table, i.e. rows and columns.
colnames(PBMC_design)    #Shows the names of columns, if any.
summary(PBMC_design)    #Outputs summary statistics of data. Ex- minimum, max, median, etc.
head(PBMC_design) #Displays first 6 rows of data.
tail(PBMC_design)      #Displays last 6 rows of data. 


##  run DE analysis
dds <- DESeqDataSetFromMatrix(countData = round(PBMC_counts),
 		colData = PBMC_design,
 	  	design = ~disease)

dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, fitType='local')
summary(dds)
resultsNames(dds) 

## extract results 
res <- results(dds,contrast=c("disease","COVID", "HD"))
head(res,2)
dim(res)

## how many genes are significant?
table(res$padj<0.05) #
res <- res[order(res$padj), ]
head(res,2)
res_data <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res_data)[1] <- "Genes"
head(res_data,2)
mcols(res,use.names = TRUE)

## how many are significant and with a log2 FC > 2 or < -2
## log2FC = log2(B) - log2(A)
## FC = 2 ^ log2FC

table(res$pvalue<0.05 & res$log2FoldChange > 2) #1951 (792) more expressed in unstimualted siCtrl
table(res$pvalue<0.05 & res$log2FoldChange < -2) #1713 (880) more expressed in TGFb siCtrl

## add annotation from ensamble gene ID to HUGO gene symbol
res_tidy <- tidy.DESeqResults(res)
res_tidy_temp2 <- res_tidy %>% inner_join(grch37, by=c("gene"="ensgene")) 
#write.csv(res_tidy_temp2, file="Covid_Cntrols_Res_annotated.csv")

##if you want to know how many protein coding genes are significant
res_tidy_temp_protein_coding = res_tidy_temp2 %>% dplyr::filter(biotype=="protein_coding")
table(res_tidy_temp_protein_coding$p.value<0.05 & res_tidy_temp_protein_coding$estimate > 2) ##518
table(res_tidy_temp_protein_coding$p.value<0.05 & res_tidy_temp_protein_coding$estimate < -2) #421

## for heatmap and batch correction you need vst or rlog
vsd <- vst(dds, blind=FALSE)

## PCA ##
DESeq2::plotPCA(vsd, intgroup="disease")+ ggtitle("disease") 

## volcano plot
res_tidy_temp_protein_coding <- as.data.frame(res_tidy_temp_protein_coding)
EnhancedVolcano(res_tidy_temp_protein_coding, lab = res_tidy_temp_protein_coding$symbol, 
                x = 'estimate',      y = 'p.adjusted', 
                xlim = c(-8, 8),      pCutoff = 0.05, 
                FCcutoff = 2,      transcriptPointSize = 1.5,      transcriptLabSize = 3.0) 


## heatmap on genes with adj pvalue < 0.05 and logFC 2
padj.cutoff <- 0.05 
lfc.cutoff <- 2
threshold <- res$padj < padj.cutoff & abs(res$log2FoldChange)> lfc.cutoff 
length(which(threshold))
res$threshold <- threshold     
sigOE <- data.frame(subset(res, threshold==TRUE))
normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[sigOE$gene,] 
select <- rownames(subset(res, threshold==TRUE))
df <- as.data.frame(colData(dds)[,c("disease")])
rownames(df)<- PBMC_design$sampleID
data <- assay(vsd)[select,]
dim(data)
head(data,2)
data <- as.data.frame(data)
counts_filtered_df <- data[apply(data, MARGIN = 1, FUN = function(x) sd(x) != 0),]
pheatmap(counts_filtered_df, show_rownames=F, annotation_col=df, fontsize = 4, fontsize_row  = 6, scale = "row",width = 2, height = 1)

### Pathway enrichment analysis - GSEA
ranks <- res_tidy_temp2$estimate
names(ranks) <- res_tidy_temp2$symbol
head(ranks, 20)

### Let’s use the Hallmark gene set from MSigDB. Hallmark gene sets summarize and represent specific well-defined biological states or processes and display coherent expression. These gene sets were generated by a computational methodology based on identifying overlaps between gene sets in other MSigDB collections and retaining genes that display coordinate expression. The gmtPathways() function will take a GMT file you downloaded from MSigDB and turn it into a list. Each element in the list is a character vector of genes in the pathway.

# Load the pathways into a named list
pathways.hallmark <- gmtPathways("../h.all.v7.0.symbols.gmt")

# Look at them all if you want (uncomment)
# pathways.hallmark

# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

## Now, run the fgsea algorithm with 1000 permutations:
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

## Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

  # Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

 ## Plot the normalized enrichment scores. Color the bar indicating whether or not the pathway was significant:
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="up/down reg pathways covid / healthy from GSEA") + 
  theme_minimal()

topUp <- fgseaRes %>% 
  dplyr::filter(ES > 0) %>% 
  top_n(10, wt=-padj)

topDown <- fgseaRes %>% 
  dplyr::filter(ES < 0) %>% 
  top_n(10, wt=-padj)

topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)

fwrite(topPathways, file ="COVID_pathways_.csv")