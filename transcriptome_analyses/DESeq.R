#!/bin/R
#### Author: KeChen		Emali:chenke57@163.com

library(DESeq2)
library(tidyverse)
setwd("./")
mycounts<-read.table("counts.tsv", header=TRUE, sep="\t")
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)
condition <- factor(c(rep("Control",3),rep("Treatment",3)), levels = c("Control","Treatment"))
condition
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
dds
res = results(dds, contrast=c("condition", "Treatment", "Control"))
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="results.csv")
table(res$padj<0.01)
diff_gene_deseq2 <-subset(res, padj < 0.01 & abs(log2FoldChange) > 2)
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_results.csv")
up_gene_deseq2 <-subset(res, padj < 0.01 & log2FoldChange > 2)
dim(up_gene_deseq2)
head(up_gene_deseq2)
write.csv(up_gene_deseq2, file= "up_results.csv")
down_gene_deseq2 <-subset(res, padj < 0.01 & log2FoldChange < -2)
dim(down_gene_deseq2)
head(down_gene_deseq2)
write.csv(down_gene_deseq2, file= "down_results.csv")
