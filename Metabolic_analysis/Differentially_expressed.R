library(DESeq2)
# Load library for RColorBrewer
library(RColorBrewer)
# Load library for pheatmap
library(pheatmap)
# Load library for tidyverse
library(tidyverse)

#setwd("..\\Data")

gene_expr = read.csv('..\\Data\\gene_count.csv', row.names = 'ï..gene_id')
gene_expr = gene_expr %>% select(c("APC_WT1","APC_WT2", "APC_WT3", "APC_m4_1", "APC_m4_2", "APC_m4_3" ))

metadata = data.frame(c('WT', 'WT', 'WT', 'M4', 'M4', 'M4'))
rownames(metadata) = colnames(gene_expr)
colnames(metadata) = c("condition")
metadata$condition <- as.factor(metadata$condition)
dds <- DESeqDataSetFromMatrix(countData= gene_expr,
                              colData = metadata,
                              design = ~condition)
dds <- DESeq(dds)

dds <- dds[rowSums(counts(dds)) >= 10,]
dds$condition <- relevel(dds$condition, ref = "M4")

res <-results(dds)
res <- res[order(res$padj),]
a <- data.frame(res)

a <- a[!is.na(a$padj),]


write.csv(a, "DEG.csv")

padj.cutoff = 0.05

upregulated <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% filter(padj < padj.cutoff & log2FoldChange > 0)
downregulated <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% filter(padj < padj.cutoff & log2FoldChange < 0)

write.csv(upregulated, "upregulated_genes.csv")
write.csv(downregulated, "downregulated_genes.csv")
