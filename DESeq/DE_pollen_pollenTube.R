# pollen vs leaf, polletube vs pollen DE analysis
library(dplyr)
library("DESeq2")

# Part1: pollen-biased genes
# prep the dataset
readcnt<-read.table("exon_cnts_TX.txt",header=T)
head(readcnt)
readcnt_matrix<-as.matrix(readcnt[ , -1])
rownames(readcnt_matrix) <- readcnt[ , 1]
head(readcnt_matrix)
str(readcnt_matrix)
meta_data <- data.frame(colnames(readcnt_matrix))
write.csv(meta_data, file="meta_data.csv", row.names = FALSE) 
meta_data<-read.csv("meta_data.csv") # labels added manually
rownames(meta_data) <- colnames(readcnt_matrix)
head(meta_data)


# run DE analysis
dds_all <- DESeqDataSetFromMatrix(countData = readcnt_matrix, 
                                  colData = meta_data,
                                  design = ~ tissue)
dds_all <- DESeq(dds_all)
resultsNames(dds_all)
res_all <- results(dds_all)
head(res_all)
summary(res_all)
result_all <- res_all[complete.cases(res_all),]
nrow(res_all) 
nrow(result_all) 


# extract pollen-biased genes
DE_tissue_TX<-as.data.frame(subset(result_all, padj<.1 & abs(log2FoldChange)>1))
DE_tissue_TX <- DE_tissue_TX %>% mutate(gene = rownames(DE_tissue_TX))
head(DE_tissue_TX) 

DE_tissue_TX <- DE_tissue_TX[,c(7,1:6)]

DE_tissue_TX_pollen <-  DE_tissue_TX %>% filter(log2FoldChange>0) 
DE_tissue_TX_leaf <-  DE_tissue_TX %>% filter(log2FoldChange<0) 
head(DE_tissue_TX_pollen)
write.csv(DE_tissue_TX_pollen, file = "pollen_biased_TX.csv", row.names = FALSE) 


# Part2: pollenTube-biased genes
# prep the dataset
readcnt<-read.table("pollenTube_cnts_TX.txt",header=T)
head(readcnt)
readcnt_matrix<-as.matrix(readcnt[ , -1])
rownames(readcnt_matrix) <- readcnt[ , 1]
head(readcnt_matrix)
str(readcnt_matrix)
meta_data <- data.frame(colnames(readcnt_matrix))
write.csv(meta_data, file="meta_data.csv", row.names = FALSE)
meta_data<-read.csv("meta_data1.csv") # labels added manually
rownames(meta_data) <- colnames(readcnt_matrix)
meta_data <- meta_data %>% select(2)
head(meta_data)


# run DE analysis
dds_all <- DESeqDataSetFromMatrix(countData = readcnt_matrix, 
                                  colData = meta_data,
                                  design = ~ tissue)
dds_all <- DESeq(dds_all)
resultsNames(dds_all)
res_all <- results(dds_all)
head(res_all)
summary(res_all)
result_all <- res_all[complete.cases(res_all),]
nrow(res_all) 
nrow(result_all) 


# extract pollenTube-biased genes
DE_tissue_TX <-
    as.data.frame(subset(result_all, padj < .1 & abs(log2FoldChange) > 1))
DE_tissue_TX <- DE_tissue_TX %>% mutate(gene = rownames(DE_tissue_TX))
head(DE_tissue_TX) 
DE_tissue_TX <- DE_tissue_TX[,c(7,1:6)]
DE_tissue_TX_pollenTube <-  DE_tissue_TX %>% filter(log2FoldChange>0) 
DE_tissue_TX_pollen <-  DE_tissue_TX %>% filter(log2FoldChange<0) 
head(DE_tissue_TX_pollen)
write.csv(DE_tissue_TX_pollenTube, file = "pollenTube_biased_TX.csv", row.names = FALSE) 

