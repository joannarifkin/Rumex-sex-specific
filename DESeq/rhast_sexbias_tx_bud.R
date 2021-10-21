## redoing previous pipeline from july using new annotation created by yunchen
## lifted over annotation - lots of lost info
library(tidyverse)
library(BiocGenerics)
library(DESeq2)
#setwd("~/Dropbox/rhast_DE_july_2021")

output_all <-read.table("flowerBud_ref_LA_exon.txt", row.names=1, header=T, skip =1, quote="")
colnames(output_all)
colnames_new <- str_replace_all(colnames(output_all), 
                            "X.ohta2.meng.yuan.rumex.pollen_expr.rhast_remap.georgeFlowerBud.AnalysisReady_",
                            "sample") %>% 
  str_replace(., ".bam","")
colnames_new
colnames(output_all) <- colnames_new 
head(output_all)
#output_fixed <- output_all %>% mutate(.,Chr = gsub(";/*","",Chr)) #filter(grepl("LG", Chr)) %>% 
#output_fixed$Chr <- gsub("LG4/*", "X", output_fixed$Chr)
#output_fixed$Chr <- gsub("LG1/*", "A1", output_fixed$Chr)
#output_fixed$Chr <- gsub("LG3/*", "A2", output_fixed$Chr)
#output_fixed$Chr <- gsub("LG5/*", "A3", output_fixed$Chr)
#output_fixed$Chr <- gsub("LG2/*", "A4", output_fixed$Chr)
#output_fixed$Chr = factor(output_fixed$Chr, levels=c('A1','A2','A4','X','A3'))
#head(output_fixed)
output_filtered <- output_all %>%
  dplyr::filter(., rowSums(across(sample17:sample24))>=20) %>% 
  select(Chr,Start,End,sample17:sample24) %>%
  rownames_to_column("GeneID")
head(output_filtered)
write.csv(output_filtered$GeneID, "expressed_genes_flower.csv")
coldata <- read.csv("FlowerBud_info.csv") 
coldata <- coldata %>% mutate(SC = paste(Sex,Cytotype, sep='')) # needed for DESEQ
##row info

rowdata <- output_filtered %>% 
  #rownames_to_column("GeneID") %>% 
  select(GeneID,Chr,Start,End) 
head(rowdata)
gene_num_chr <- rowdata %>% dplyr::count(Chr)
gene_num_chr
#write.csv(gene_num_chr,"expressed_gene_num_chr.csv", quote =F)

counts_tx <- output_filtered %>% column_to_rownames("GeneID") %>% select(sample17:sample24) 
coldata_tx <- coldata %>% slice_tail(n=8) 
df_tx <- DESeqDataSetFromMatrix(counts_tx, coldata_tx, design = ~Sex)
dds_tx <-DESeq(df_tx)
resultsNames(dds_tx) # sex_M vs F --> LFC >0 is male biased, <0 is female-biased.
res_tx <- results(dds_tx)
summary(res_tx) #padj<.1 & abs(log2FoldChange)>1)
resSig_tx <- subset(res_tx, padj<0.1 & abs(log2FoldChange)>1)
dim(resSig_tx) # 850 genes?
head(resSig_tx)
resSig_tx_ord <- as.data.frame(resSig_tx) %>% rownames_to_column("GeneID") %>%
  arrange(desc(log2FoldChange))
head(resSig_tx_ord)
sexDEGS_tx <- inner_join(resSig_tx_ord, rowdata, "GeneID", ) %>%## list of DE genes
  mutate(Bias=(ifelse(log2FoldChange > 0, 'MaleBias', 'FemaleBias')))  #%>% filter(grep
##
write.csv(sexDEGS_tx,"sex_biased_flowerbud_TX.csv", quote =F)
