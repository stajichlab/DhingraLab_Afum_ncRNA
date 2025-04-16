#!/usr/bin/env Rscript

library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(fdrtool)
library(geneplotter)
library(EDASeq)
library(tidyverse)

my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))

#countdata <- read.table("results/read_count_A1163", header=TRUE, row.names=1)
countdata_all <- read_tsv("results/read_count_A1163",skip=1)  %>% 
  rename_all(~(stringr::str_replace_all(.,"\\.Aligned\\.out\\.bam", "") %>%
               stringr::str_replace_all(.,"results/STAR_A1163/","")))

#countdata <- countdata_all %>% dplyr::select(-c(Chr,Start,End,Strand,Length))
# keep WT and remove 30C
WTcount <- countdata_all %>% column_to_rownames(var="Geneid") %>% dplyr::select(contains('WT_')) %>% dplyr::select(!contains('_30_'))
DELcount <- countdata_all %>% column_to_rownames(var="Geneid") %>% dplyr::select(contains('Del_')) %>% dplyr::select(!contains('_30_'))

samples <- read_csv("metadata.csv")

# WT first
WTsamples <- samples %>% dplyr::filter(Genotype=="WT" & Condition != 30)
WTexprnames <- WTsamples %>% dplyr::select(c(SAMPLE))
# remove 30C

# check that experimental columns match in order
# all(exprnames %in% colnames(countdata))
all(WTexprnames == colnames(WTcount))

#countdata <- countdata[,exprnames]
#all(exprnames == colnames(countdata))

# DEseq2 analyses
# ignoring geno for these pairwise
#geno = factor( samples$Genotype)
rep = factor( WTsamples$Replicate)
treatment = factor (WTsamples$Condition)
#genotype = geno,
sampleTable <- data.frame(replicate = rep,
                          condition = treatment)
# check this is right order!
rownames(sampleTable) = colnames(WTcount)

dds <- DESeqDataSetFromMatrix(countData = WTcount,
                              colData   = sampleTable, 
                              design    = ~ condition )

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RNASeq_WT_Temp.pdf")

plotDispEsts(dds)

multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))

MA.idx = t(combn(1:4, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(dds, normalized = T), 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
}

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]


#datCollapsed <- collapseReplicates(dds, groupby=dds$genotype,run=dds$replicate,renameCols=TRUE)

df2 <- as.data.frame(colData(dds)[,c("condition")])
rownames(df2) = colnames(WTcount)
colnames(df2) = c("Temperature")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Top Expression")

topVar <- head(order(rowVars(assay(vsd)),
                     decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]

pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD WT Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD WT Top Expression")

pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD WT Top Expression")

topVar <- head(order(rowVars(assay(rld)),
                     decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD WT Most different")


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pdf("plots/WT_PCA_expresion.pdf")
pcaData <- plotPCA(vsd, intgroup=c("condition","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=treatment,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
mat<-log.norm.counts[topVarGenes,]
mat<-mat -rowMeans(mat)

pdf("plots/WT_RNASeq_heatmap_allTemp.pdf")
pheatmap(mat,method="complete",main = "TopVar normalized", show_rownames = F,
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )

topVarGenes1 <- order(-rowVars(assay(rld)))[0:2000]
mat1 <- assay(rld)[ topVarGenes1, ]
mat1<- mat1 - rowMeans(mat1)

pheatmap(mat1, method="complete",
         main = "Unsupervised WT 2000 genes ",
         show_rownames = F,annotation_legend = FALSE, 
         legend=T, cluster_cols=TRUE)

dds <- DESeq(dds)
res <- results(dds)
#res    
#res <- results(dds, contrast=c("condition","25","37"))
#res

resultsNames(dds)
#resLFC <- lfcShrink(dds, coef="condition_37_vs_25", type="apeglm")
resLFC <- lfcShrink(dds, contrast=c("condition","25","37"), type="ashr")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
plotMA(res05, ylim=c(-2,2),main="25 vs 37")
plotMA(resLFC, ylim=c(-2,2),main="25 vs 37")

resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
summary(resSig)
write.csv(resSig,"reports/WT_25_vs_37.csv")
write.csv(fpm(dds),"reports/WT_FPM.csv")

# DRAW 25 VS 37 HEATMAP

df2 <- as.data.frame(colData(dds)[,c("condition")])
                     
rownames(df2) = colnames(WTcount)
colnames(df2) = c("Temperature")
df2$Temperature = paste0(df2$Temperature, "C")

vsdvals_25_37 <- as_tibble(assay(vsd)) %>% dplyr::select(contains("_25_") | contains("_37_")) %>% 
  add_column(GeneID = rownames(assay(vsd)))

vsdvals_25_37f <- vsdvals_25_37 %>% filter(GeneID %in% rownames(resSig)) %>% column_to_rownames(var="GeneID")
pheatmap(vsdvals_25_37f, cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/WT_RNAseq_25_37_heatmap.pdf",
         main="WT 25 vs 37 Expression")

vsdvals_25_37f <- vsdvals_25_37f %>% arrange(WT_25_1)
pheatmap(vsdvals_25_37f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/WT_RNAseq_25_37_heatmap_25ordered.pdf",
         main="WT 25 vs 37 Expression ordered by 25C")


res <- results(dds, contrast=c("condition","25","42"))
res

resultsNames(dds)
#resLFC <- lfcShrink(dds, coef="condition_25_vs_42", type="apeglm")
resLFC <- lfcShrink(dds, contrast=c("condition","25","42"), type="ashr")

resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01, na.rm=TRUE)

plotMA(res, ylim=c(-2,2),main="25 vs 42")
plotMA(resLFC, ylim=c(-2,2),main="25 vs 42")

resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
summary(resSig)
write.csv(resSig,"reports/WT_25_vs_42.csv")

# DRAW 25 VS 42 HEATMAP

vsdvals_25_42 <- as_tibble(assay(vsd)) %>% dplyr::select(contains("_25_") | contains("_42_")) %>% 
  add_column(GeneID = rownames(assay(vsd)))

vsdvals_25_42f <- vsdvals_25_42 %>% filter(GeneID %in% rownames(resSig)) %>% column_to_rownames(var="GeneID")

df2 <- as.data.frame(colData(dds)[,c("condition")])
rownames(df2) = colnames(WTcount)
colnames(df2) = c("Temperature")
pheatmap(vsdvals_25_42f, cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,filename="plots/WT_RNAseq_25_42_heatmap.pdf",
         cluster_cols=FALSE, annotation_col=df2,main="WT 25 vs 42")

vsdvals_25_42f <- vsdvals_25_42f %>% arrange(WT_25_1)
pheatmap(vsdvals_25_42f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/WT_RNAseq_25_42_heatmap_ordered.pdf",
         main="WT 25 vs 42 Expression ordered by 25C")


#res <- results(dds, contrast=c("condition","37","42"))
#resultsNames(dds)

resLFC <- lfcShrink(dds, contrast=c("condition","37","42"), type="ashr")

resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01, na.rm=TRUE)
plotMA(res, ylim=c(-2,2),main="37 vs 42")
plotMA(resLFC, ylim=c(-2,2),main="37 vs 42")

# Get diff expressed
resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/WT_37_vs_42.csv")

df2 <- as.data.frame(colData(dds)[,c("condition")])

rownames(df2) = colnames(WTcount)
colnames(df2) = c("Temperature")
df2$Temperature = paste0(df2$Temperature, "C")

vsdvals_37_42 <- as_tibble(assay(vsd)) %>% dplyr::select(contains("_37_") | contains("_42_")) %>% 
  add_column(GeneID = rownames(assay(vsd)))


vsdvals_37_42f <- vsdvals_37_42 %>% filter(GeneID %in% rownames(resSig)) %>%column_to_rownames(var="GeneID")
pheatmap(vsdvals_37_42f, cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/WT_RNAseq_37_42_heatmap.pdf",
         main="WT 37 vs 42 Expression")

vsdvals_37_42f <- vsdvals_37_42f %>% arrange(WT_37_1)
pheatmap(vsdvals_37_42f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/WT_RNAseq_37_42_heatmap_ordered.pdf",
         main="WT 37 vs 42 Expression order by 37C")



# Now DEL genotype
DELsamples <- samples %>% dplyr::filter(Genotype=="Del" & Condition != 30)
DELexprnames <- DELsamples %>% dplyr::select(c(SAMPLE))
# remove 30C

# check that experimental columns match in order
# all(exprnames %in% colnames(countdata))
all(DELexprnames == colnames(DELcount))

#countdata <- countdata[,exprnames]
#all(exprnames == colnames(countdata))

# DEseq2 analyses
# ignoring geno for these pairwise
#geno = factor( samples$Genotype)
rep = factor( DELsamples$Replicate)
treatment = factor (DELsamples$Condition)
#genotype = geno,
sampleTable <- data.frame(replicate = rep,
                          condition = treatment)
# check this is right order!
rownames(sampleTable) = colnames(DELcount)

dds <- DESeqDataSetFromMatrix(countData = DELcount,
                              colData   = sampleTable, 
                              design    = ~ condition )

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RNASeq_DEL_Temp.pdf")

plotDispEsts(dds)

multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))

MA.idx = t(combn(1:4, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(dds, normalized = T), 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
}

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]


#datCollapsed <- collapseReplicates(dds, groupby=dds$genotype,run=dds$replicate,renameCols=TRUE)

df2 <- as.data.frame(colData(dds)[,c("condition")])
rownames(df2) = colnames(DELcount)
colnames(df2) = c("Temperature")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="DEL VSD Top Expression")

topVar <- head(order(rowVars(assay(vsd)),
                     decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]

pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="DEL VSD Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="DEL RLD Top Expression")

pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="DEL RLD Top Expression")

topVar <- head(order(rowVars(assay(rld)),
                     decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="DEL RLD Most different")


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pdf("plots/DEL_PCA_expresion.pdf")
pcaData <- plotPCA(vsd, intgroup=c("condition","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=treatment,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
mat<-log.norm.counts[topVarGenes,]
mat<-mat -rowMeans(mat)

pheatmap(mat,method="complete",main = "DEL TopVar normalized", show_rownames = F, filename="plots/DEL_RNASeq_heatmap_allTemp.pdf",
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )

topVarGenes1 <- order(-rowVars(assay(rld)))[0:2000]
mat1 <- assay(rld)[ topVarGenes1, ]
mat1<- mat1 - rowMeans(mat1)

pheatmap(mat1, method="complete",
         main = "Unsupervised 2000 genes ",
         show_rownames = F,annotation_legend = FALSE, filename="plots/DEL_RNASeq_heatmap_2000_allTemp.pdf",
         legend=T, cluster_cols=TRUE)

dds <- DESeq(dds)
res <- results(dds)
#res    
#res <- results(dds, contrast=c("condition","25","37"))
#res

resultsNames(dds)
#resLFC <- lfcShrink(dds, coef="condition_37_vs_25", type="apeglm")
resLFC <- lfcShrink(dds, contrast=c("condition","25","37"), type="ashr")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
plotMA(res05, ylim=c(-2,2),main="DEL 25 vs 37")
plotMA(resLFC, ylim=c(-2,2),main="DEL 25 vs 37")

resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
summary(resSig)
write.csv(resSig,"reports/DEL_25_vs_37.csv")
write.csv(fpm(dds),"reports/DEL_FPM.csv")

# DRAW 25 VS 37 HEATMAP

df2 <- as.data.frame(colData(dds)[,c("condition")])

rownames(df2) = colnames(DELcount)
colnames(df2) = c("Temperature")
df2$Temperature = paste0(df2$Temperature, "C")

vsdvals_25_37 <- as_tibble(assay(vsd)) %>% dplyr::select(contains("_25_") | contains("_37_")) %>% 
  add_column(GeneID = rownames(assay(vsd)))

vsdvals_25_37f <- vsdvals_25_37 %>% filter(GeneID %in% rownames(resSig)) %>% column_to_rownames(var="GeneID")
pheatmap(vsdvals_25_37f, cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2, filename="plots/DEL_RNAseq_25_37_heatmap.pdf",
         main="DEL 25 vs 37 Expression")

vsdvals_25_37f <- vsdvals_25_37f %>% arrange(Del_25_1)
pheatmap(vsdvals_25_37f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/DEL_RNAseq_25_37_heatmap_25Cordered.pdf",
         main="DEL 25 vs 37 Expression ordered by 25C")


#res <- results(dds, contrast=c("condition","25","42"))
#res

resultsNames(dds)
#resLFC <- lfcShrink(dds, coef="condition_42_vs_25", type="apeglm")
resLFC <- lfcShrink(dds, contrast=c("condition","25","42"), type="ashr")

resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01, na.rm=TRUE)
plotMA(res, ylim=c(-2,2),main="DEL 42 vs 25")
plotMA(resLFC, ylim=c(-2,2),main="DEL 42 vs 25")

# DRAW 25 VS 42 HEATMAP

vsdvals_25_42 <- as_tibble(assay(vsd)) %>% dplyr::select(contains("_25_") | contains("_42_")) %>% 
  add_column(GeneID = rownames(assay(vsd)))

vsdvals_25_42f <- vsdvals_25_42 %>% filter(GeneID %in% rownames(resSig)) %>% column_to_rownames(var="GeneID")

df2 <- as.data.frame(colData(dds)[,c("condition")])
rownames(df2) = colnames(DELcount)
colnames(df2) = c("Temperature")
pheatmap(vsdvals_25_42f, cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,filename="plots/DEL_RNAseq_25_42_heatmap.pdf",
         cluster_cols=FALSE, annotation_col=df2,main="DEL 25 vs 42")

vsdvals_25_42f <- vsdvals_25_42f %>% arrange(Del_25_1)
pheatmap(vsdvals_25_42f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/DEL_RNAseq_25_42_heatmap_25ordered.pdf",
         main="DEL 25 vs 42 Expression ordered by 25C")


# Get diff expressed
resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/DEL_25_vs_42.csv")

res <- results(dds, contrast=c("condition","37","42"))
resultsNames(dds)

resLFC <- lfcShrink(dds, contrast=c("condition","37","42"), type="ashr")

resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01, na.rm=TRUE)
plotMA(res, ylim=c(-2,2),main="DEL 37 vs 42")
plotMA(resLFC, ylim=c(-2,2),main="DEL 37 vs 42")

# Get diff expressed
resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/DEL_37_vs_42.csv")

df2 <- as.data.frame(colData(dds)[,c("condition")])

rownames(df2) = colnames(DELcount)
colnames(df2) = c("Temperature")
df2$Temperature = paste0(df2$Temperature, "C")

vsdvals_37_42 <- as_tibble(assay(vsd)) %>% dplyr::select(contains("_37_") | contains("_42_")) %>% 
  add_column(GeneID = rownames(assay(vsd)))


vsdvals_37_42f <- vsdvals_37_42 %>% filter(GeneID %in% rownames(resSig)) %>%column_to_rownames(var="GeneID")
pheatmap(vsdvals_37_42f, cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/DEL_RNAseq_37_42_heatmap.pdf",
         main="DEL 37 vs 42 Expression")

vsdvals_37_42f <- vsdvals_37_42f %>% arrange(Del_37_1)
pheatmap(vsdvals_37_42f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/DEL_RNAseq_37_42_heatmap_37ordered.pdf",
         main="DEL 37 vs 42 Expression order by 37C")



# now same temperature where genotype is the variable

Temp25count <- countdata_all %>% column_to_rownames(var="Geneid") %>% dplyr::select(contains('_25_'))

Temp25Samples <- samples %>% dplyr::filter(Condition == 25)
Temp25exprnames <- Temp25Samples %>% dplyr::select(c(SAMPLE))

# check that experimental columns match in order
# all(exprnames %in% colnames(countdata))
all(Temp25exprnames == colnames(Temp25count))

#countdata <- countdata[,exprnames]
#all(exprnames == colnames(countdata))

# DEseq2 analyses
# ignoring temp for these pairwise
geno = factor( Temp25Samples$Genotype)
rep = factor( Temp25Samples$Replicate)
treatment = factor (Temp25Samples$Condition)
sampleTable <- data.frame(replicate = rep,
                          genotype = geno)
# check this is right order!
rownames(sampleTable) = colnames(Temp25count)

dds <- DESeqDataSetFromMatrix(countData = Temp25count,
                              colData   = sampleTable, 
                              design    = ~ genotype )

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RNASeq_25Temp_GenoCompare.pdf")

plotDispEsts(dds)

multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))

MA.idx = t(combn(1:4, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(dds, normalized = T), 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
}

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]


#datCollapsed <- collapseReplicates(dds, groupby=dds$genotype,run=dds$replicate,renameCols=TRUE)

df2 <- as.data.frame(colData(dds)[,c("genotype")])
rownames(df2) = colnames(Temp25count)
colnames(df2) = c("Temperature")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 25 VSD Top Expression")

topVar <- head(order(rowVars(assay(vsd)),
                     decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]

pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 25 VSD Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 25 RLD Top Expression")

pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 25 RLD Top Expression")

topVar <- head(order(rowVars(assay(rld)),
                     decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 25 RLD Most different")


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pdf("plots/RNASeq_25Temp_PCA_expresion.pdf")
pcaData <- plotPCA(vsd, intgroup=c("genotype","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
mat<-log.norm.counts[topVarGenes,]
mat<-mat -rowMeans(mat)

pheatmap(mat,method="complete",main = "Temp 25 TopVar normalized", show_rownames = F, filename="plots/Temp25_RNASeq_heatmap_allTemp.pdf",
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )

topVarGenes1 <- order(-rowVars(assay(rld)))[0:2000]
mat1 <- assay(rld)[ topVarGenes1, ]
mat1<- mat1 - rowMeans(mat1)

pheatmap(mat1, method="complete",
         main = "Unsupervised 2000 genes ",
         show_rownames = F,annotation_legend = FALSE, filename="plots/Temp25_RNASeq_heatmap_2000_allTemp.pdf",
         legend=T, cluster_cols=TRUE)

dds <- DESeq(dds)
res <- results(dds)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="genotype_WT_vs_Del", type="apeglm")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
plotMA(res05, ylim=c(-2,2),main="Temp 25 WT vs Del")
plotMA(resLFC, ylim=c(-2,2),main="Temp 25 WT vs Del")

resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
summary(resSig)
write.csv(resSig,"reports/Temp25_WT_vs_Del.csv")
write.csv(fpm(dds),"reports/Temp25_FPM.csv")

# DRAW 25 WT vs Del HEATMAP

df2 <- as.data.frame(colData(dds)[,c("genotype")])

rownames(df2) = colnames(Temp25count)
colnames(df2) = c("Genotype")

vsdvals <- as_tibble(assay(vsd)) %>% add_column(GeneID = rownames(assay(vsd)))

vsdvals_f <- vsdvals %>% filter(GeneID %in% rownames(resSig)) %>% column_to_rownames(var="GeneID")
pheatmap(vsdvals_f, cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2, filename="plots/Temp25_RNAseq_WT_Del_heatmap.pdf",
         main="Temp 25 WT vs Del Expression")

vsdvals_f <- vsdvals_f %>% arrange(WT_25_1)
pheatmap(vsdvals_f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/Temp25_RNAseq_WT_Del_heatmap_WTordered.pdf",
         main="Temp 25 WT vs Del Expression Ordered by WT")


# now same temperature where genotype is the variable -- 37

Temp37count <- countdata_all %>% column_to_rownames(var="Geneid") %>% dplyr::select(contains('_37_'))

Temp37Samples <- samples %>% dplyr::filter(Condition == 37)
Temp37exprnames <- Temp25Samples %>% dplyr::select(c(SAMPLE))

# check that experimental columns match in order
# all(exprnames %in% colnames(countdata))
all(Temp37exprnames == colnames(Temp37count))

geno = factor( Temp37Samples$Genotype)
rep = factor( Temp37Samples$Replicate)
treatment = factor (Temp37Samples$Condition)
sampleTable <- data.frame(replicate = rep,
                          genotype = geno)
# check this is right order!
rownames(sampleTable) = colnames(Temp37count)

dds <- DESeqDataSetFromMatrix(countData = Temp37count,
                              colData   = sampleTable, 
                              design    = ~ genotype )

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RNASeq_37Temp_GenoCompare.pdf")

plotDispEsts(dds)

multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))

MA.idx = t(combn(1:4, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(dds, normalized = T), 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
}

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]

df2 <- as.data.frame(colData(dds)[,c("genotype")])
rownames(df2) = colnames(Temp37count)
colnames(df2) = c("Temperature")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 37 VSD Top Expression")

topVar <- head(order(rowVars(assay(vsd)),
                     decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]

pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 37 VSD Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 37 RLD Top Expression")

pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 37 RLD Top Expression")

topVar <- head(order(rowVars(assay(rld)),
                     decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 37 RLD Most different")


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pdf("plots/RNASeq_37Temp_PCA_expresion.pdf")
pcaData <- plotPCA(vsd, intgroup=c("genotype","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
mat<-log.norm.counts[topVarGenes,]
mat<-mat -rowMeans(mat)

pheatmap(mat,method="complete",main = "Temp 37 TopVar normalized", show_rownames = F, filename="plots/Temp37_RNASeq_heatmap_allTemp.pdf",
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )

topVarGenes1 <- order(-rowVars(assay(rld)))[0:2000]
mat1 <- assay(rld)[ topVarGenes1, ]
mat1<- mat1 - rowMeans(mat1)

pheatmap(mat1, method="complete",
         main = "Unsupervised 2000 genes ",
         show_rownames = F,annotation_legend = FALSE, filename="plots/Temp37_RNASeq_heatmap_2000_allTemp.pdf",
         legend=T, cluster_cols=TRUE)

dds <- DESeq(dds)
res <- results(dds)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="genotype_WT_vs_Del", type="apeglm")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
plotMA(res05, ylim=c(-2,2),main="Temp 37 WT vs Del")
plotMA(resLFC, ylim=c(-2,2),main="Temp 37 WT vs Del")

resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
summary(resSig)
write.csv(resSig,"reports/Temp37_WT_vs_Del.csv")
write.csv(fpm(dds),"reports/Temp37_FPM.csv")

# DRAW 37 WT vs Del HEATMAP

df2 <- as.data.frame(colData(dds)[,c("genotype")])

rownames(df2) = colnames(Temp37count)
colnames(df2) = c("Genotype")

vsdvals <- as_tibble(assay(vsd)) %>% add_column(GeneID = rownames(assay(vsd)))

vsdvals_f <- vsdvals %>% filter(GeneID %in% rownames(resSig)) %>% column_to_rownames(var="GeneID")
pheatmap(vsdvals_f, cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2, filename="plots/Temp37_RNAseq_WT_Del_heatmap.pdf",
         main="Temp 37 WT vs Del Expression")

vsdvals_f <- vsdvals_f %>% arrange(WT_37_1)
pheatmap(vsdvals_f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/Temp37_RNAseq_WT_Del_heatmap_WTordered.pdf",
         main="Temp 37 WT vs Del Expression Ordered by WT")


# now same temperature where genotype is the variable -- 42

Temp42count <- countdata_all %>% column_to_rownames(var="Geneid") %>% dplyr::select(contains('_42_'))

Temp42Samples <- samples %>% dplyr::filter(Condition == 42)
Temp42exprnames <- Temp42Samples %>% dplyr::select(c(SAMPLE))

# check that experimental columns match in order
# all(exprnames %in% colnames(countdata))
all(Temp42exprnames == colnames(Temp42count))

#countdata <- countdata[,exprnames]
#all(exprnames == colnames(countdata))

# DEseq2 analyses
# ignoring temp for these pairwise
geno = factor( Temp42Samples$Genotype)
rep = factor( Temp42Samples$Replicate)
treatment = factor (Temp42Samples$Condition)
sampleTable <- data.frame(replicate = rep,
                          genotype = geno)
# check this is right order!
rownames(sampleTable) = colnames(Temp42count)

dds <- DESeqDataSetFromMatrix(countData = Temp42count,
                              colData   = sampleTable, 
                              design    = ~ genotype )

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RNASeq_42Temp_GenoCompare.pdf")

plotDispEsts(dds)

multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))

MA.idx = t(combn(1:4, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(dds, normalized = T), 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
}

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]


df2 <- as.data.frame(colData(dds)[,c("genotype")])
rownames(df2) = colnames(Temp42count)
colnames(df2) = c("Temperature")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 42 VSD Top Expression")

topVar <- head(order(rowVars(assay(vsd)),
                     decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]

pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 42 VSD Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 42 RLD Top Expression")

pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 42 RLD Top Expression")

topVar <- head(order(rowVars(assay(rld)),
                     decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Temp 42 RLD Most different")


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pdf("plots/RNASeq_42Temp_PCA_expresion.pdf")
pcaData <- plotPCA(vsd, intgroup=c("genotype","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
mat<-log.norm.counts[topVarGenes,]
mat<-mat -rowMeans(mat)

pheatmap(mat,method="complete",main = "Temp 42 TopVar normalized", show_rownames = F, filename="plots/Temp42_RNASeq_heatmap_allTemp.pdf",
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )

topVarGenes1 <- order(-rowVars(assay(rld)))[0:2000]
mat1 <- assay(rld)[ topVarGenes1, ]
mat1<- mat1 - rowMeans(mat1)

pheatmap(mat1, method="complete",
         main = "Unsupervised 2000 genes ",
         show_rownames = F,annotation_legend = FALSE, filename="plots/Temp42_RNASeq_heatmap_2000_allTemp.pdf",
         legend=T, cluster_cols=TRUE)

dds <- DESeq(dds)
res <- results(dds)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="genotype_WT_vs_Del", type="apeglm")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
plotMA(res05, ylim=c(-2,2),main="Temp 42 WT vs Del")
plotMA(resLFC, ylim=c(-2,2),main="Temp 42 WT vs Del")

resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) >= 1 & baseMean >= 5)
resSig <- resSig[order(resSig$pvalue),]
summary(resSig)
write.csv(resSig,"reports/Temp42_WT_vs_Del.csv")
write.csv(fpm(dds),"reports/Temp42_FPM.csv")

# DRAW 42 WT vs Del HEATMAP

df2 <- as.data.frame(colData(dds)[,c("genotype")])

rownames(df2) = colnames(Temp42count)
colnames(df2) = c("Genotype")

vsdvals <- as_tibble(assay(vsd)) %>% add_column(GeneID = rownames(assay(vsd)))

vsdvals_f <- vsdvals %>% filter(GeneID %in% rownames(resSig)) %>% column_to_rownames(var="GeneID")
pheatmap(vsdvals_f, cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2, filename="plots/Temp42_RNAseq_WT_Del_heatmap.pdf",
         main="Temp 42 WT vs Del Expression")

vsdvals_f <- vsdvals_f %>% arrange(WT_42_1)
pheatmap(vsdvals_f, cluster_rows=FALSE, show_rownames=FALSE,
         fontsize_row = 4,fontsize_col = 7,
         cluster_cols=FALSE, 
         annotation_col=df2,filename="plots/Temp42_RNAseq_WT_Del_heatmap_WTordered.pdf",
         main="Temp 42 WT vs Del Expression Ordered by WT")




