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

my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))

countdata <- read.table("results/read_count_A1163", header=TRUE, row.names=1)
colnames(countdata) <- gsub("results\\.STAR_A1163\\.", "", colnames(countdata))
colnames(countdata) <- gsub("\\.Aligned\\.out\\.bam", "", colnames(countdata))

countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)
head(countdata)

samples <- read.csv("metadata.csv",header=TRUE)
exprnames <- samples["SAMPLE"]
#exprnames <- sub("","",exprnames,perl=FALSE)

# check that experimental columns match in order
all(exprnames %in% colnames(countdata))
all(exprnames == colnames(countdata))
# reorder the columns
#countdata <- countdata[,exprnames]
all(exprnames == colnames(countdata))

# DEseq2 analyses
geno = factor( samples$Genotype)
rep = factor( samples$Replicate)
treatment = factor (samples$Condition)

sampleTable <- data.frame(genotype = geno,
                          replicate = rep,
                          condition = treatment)
# check this is right order!
rownames(sampleTable) = colnames(countdata)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData   = sampleTable, 
                              design    = genotype ~ condition )

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
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RNASeq_Temp_Genotype.pdf")

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

df2 <- as.data.frame(colData(dds)[,c("condition","genotype")])
rownames(df2) = colnames(countdata)
colnames(df2) = c("Temperature","Genotype")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Top Expression")

topVar <- head(order(rowVars(assay(vsd)),
                     decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]

pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")

pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")

topVar <- head(order(rowVars(assay(rld)),
                     decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Most different")



sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pdf("plots/PCA_expresion.pdf")
pcaData <- plotPCA(vsd, intgroup=c("genotype","condition","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# FIX ME HERE
ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=treatment,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=replicate,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=replicate,label=treatment)) +
  geom_point(size=3) + geom_text(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()


norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
mat<-log.norm.counts[topVarGenes,]
mat<-mat -rowMeans(mat)

pheatmap(mat,method="complete",main = "TopVar normalized", show_rownames = F,
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )

topVarGenes1 <- order(-rowVars(assay(rld)))[0:1000]
mat1 <- assay(rld)[ topVarGenes1, ]
mat1<- mat1 - rowMeans(mat1)

pheatmap(mat1, method="complete",
         main = "Unsupervised 1000 genes ",
         show_rownames = F,annotation_legend = FALSE, 
         legend=T, cluster_cols=TRUE)

dds <- DESeq(dds)
res <- results(dds)
res    
res <- results(dds, contrast=c("condition","25","30"))
res

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_30_vs_25", type="apeglm")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01, na.rm=TRUE)
plotMA(res, ylim=c(-2,2),main="30 vs 25")
plotMA(resLFC, ylim=c(-2,2),main="30 vs 25")

resSig <- subset(resLFC, padj < 0.05)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/all_30_vs_25.csv")
write.csv(fpm(dds),"reports/all_FPM.csv")

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_37_vs_25", type="apeglm")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01, na.rm=TRUE)
plotMA(res, ylim=c(-2,2),main="37 vs 25")
plotMA(resLFC, ylim=c(-2,2),main="37 vs 25")

resSig <- subset(resLFC, padj < 0.05)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/all_37_vs_25.csv")

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_42_vs_25", type="apeglm")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01, na.rm=TRUE)
plotMA(res, ylim=c(-2,2),main="42 vs 25")
plotMA(resLFC, ylim=c(-2,2),main="42 vs 25")

# Get diff expressed
resSig <- subset(resLFC, padj < 0.05)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/all_42_vs_25.csv")

topresSig <- as.data.frame(head(resSig$log2FoldChange,50))
colnames(topresSig) = c("log2FoldChange")
pheatmap(topresSig, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,
         fontsize_row = 6,fontsize_col = 7,width=5,cellwidth=10,
         cluster_cols=FALSE, main="Log(2) fold change of 42 vs 25")


topChange <- head(resLFC[order(resLFC$log2FoldChange,decreasing=FALSE),],100)
topChangePlot =  data.frame(log2FoldChange = topChange$log2FoldChange)
rownames(topChangePlot) = rownames(topChange)
pheatmap(topChangePlot,method="complete",main = "Fold Change Plot Down", show_rownames = T,show_colnames=F,
         cluster_rows=FALSE,cluster_cols=FALSE,cexRow=0.3,legend=T, cexRow=0.3,
         fontsize_row = 6)
