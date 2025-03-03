#!/usr/bin/env Rscript

library(tidyverse)
library(fs)
library(cowplot)
library(DESeq2)
library(tximport)
library(dplyr)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(fdrtool)
library(geneplotter)
library(EDASeq)
library(dplyr)
unpaired_folder <- "results/ncRNA_bbmap_rpkm"
rpkmfiles <- fs::dir_ls(unpaired_folder,regexp = "\\.merge_FR\\.rpkm\\.txt$")

rpkmdata <- rpkmfiles %>% map_dfr(read_tsv,.id = "source", skip=4) %>% 
  mutate(Name = str_split_i(`#Name`,' ',1),
         sourcedata = str_split_i(str_split_i(source,'/',3),'\\.',1)) %>% 
  dplyr::select(c(sourcedata,Length,Bases,Coverage,Name,Reads,RPKM,Frags,FPKM))

ncRNA <- rpkmdata %>% filter(grepl('Afu-',Name) )

libcount <- rpkmdata %>% group_by(sourcedata) %>% summarize(sum(Reads)) %>% 
  mutate(CT=`sum(Reads)`) %>% dplyr::select(c(sourcedata,CT))

# fix at 
bygene <- ncRNA %>% dplyr::select(c(sourcedata,Length,Name,Reads)) %>% filter(grepl("^WT",sourcedata)) %>%
  pivot_wider(names_from = sourcedata, values_from=Reads)

rownames(bygene) <-bygene$Name


p0 <- ggplot(data=ncRNA %>% filter(Name != "Afu-254") %>% filter(grepl("^WT",sourcedata)), 
             aes(x=sourcedata, y=FPKM, fill=Name)) +
  geom_bar(stat="identity", width=0.75) + 
  theme_cowplot(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p0

ggsave("plots/ncRNA_WT_merge_FPKM.pdf",p0,width=8,height=8)


p0A <- ggplot(data=ncRNA %>% filter(Name == "Afu-182")  %>% filter(grepl("^WT",sourcedata)), 
             aes(x=sourcedata, y=FPKM, fill=Name)) +
  geom_bar(stat="identity",width=0.75) + 
  theme_cowplot(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

p0A

ggsave("plots/Afu-182_WT_merge_FPKM.pdf",p0A,width=8,height=8)

my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))


# fix at 
WT_bygene <- rpkmdata %>% dplyr::select(c(sourcedata,Name,Reads)) %>% filter(grepl("^WT",sourcedata)) %>%
  pivot_wider(names_from = sourcedata, values_from=Reads)

countdata <- as.matrix(WT_bygene %>% dplyr::select(-c(Name)))
rownames(countdata) <- WT_bygene$Name
head(countdata)

samples <- read.csv("metadata.csv",header=TRUE) %>% filter(grepl("^WT",SAMPLE))
exprnames <- samples["SAMPLE"]

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
                              design    = ~condition )

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
  as.data.frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/Af293_DESeq2.pdf")
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

p<- ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
p
ggsave("plots/hexbin.pdf",p)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]


#datCollapsed <- collapseReplicates(dds, groupby=dds$genotype,run=dds$replicate,renameCols=TRUE)

df2 <- as.data.frame(colData(dds)[,c("condition")])
rownames(df2) = colnames(countdata)
colnames(df2) = c("Temperature")
p<-pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Top Expression")

ggsave("plots/heatmap_VSD.pdf",p,width=8,height=8)
# just get the ncRNA
select2 = substr( rownames(counts(dds,normalized=TRUE)),1,4) == "Afu-"
ncRNA_only <- assay(vsd)[select2,]
select3 = rownames(ncRNA_only) != "Afu-254"
ncRNA_only2 <- ncRNA_only[select3,]
p<-pheatmap(ncRNA_only2, cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2, main="VSD ncRNA Expression")
ggsave("plots/heatmap_ncRNA_VSD.pdf",p,width=8,height=8)
topVar <- head(order(rowVars(assay(vsd)),
                     decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]

p<-pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Most different")
ggsave("plots/heatmap_VSD_mostdifferent.pdf",p,width=8,height=8)

p<-pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")
ggsave("plots/heatmap_RLD_topexpression.pdf",p,width=8,height=8)

p<-pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=TRUE, annotation_col=df2,main="RLD Top Expression")
ggsave("plots/heatmap_RLD_cluster_topexpression.pdf",p,width=8,height=8)
topVar <- head(order(rowVars(assay(rld)),
                     decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
p<-pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Most different")
ggsave("plots/heatmap_RLD_cluster_mostdiff.pdf",p,width=8,height=8)


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, sep="_")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

p<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
ggsave("plots/compare_similarity.pdf",p)
norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

topVarGenes <- order(-rowVars(log.norm.counts)[0:100])
mat<-log.norm.counts[topVarGenes,]
mat<-mat -rowMeans(mat)

p<-pheatmap(mat,method="complete",main = "TopVar normalized", show_rownames = F,
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )
ggsave("plots/heatmap_TopVariation_normalized.pdf",p,width=16,height=16)
topVarGenes1 <- order(-rowVars(assay(rld)))[0:1000]
mat1 <- assay(rld)[ topVarGenes1, ]
mat1<- mat1 - rowMeans(mat1)

p<-pheatmap(mat1, method="complete",
         main = "Unsupervised 1000 genes ",
         show_rownames = F,annotation_legend = FALSE, 
         legend=T, cluster_cols=TRUE)
ggsave("plots/heatmap_Unsupervised_top1k_normalized.pdf",p,width=16,height=16)
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

pdf("plots/MA_plots_af293_bbmap_combined.pdf")

plotMA(res, ylim=c(-2,2),main="30 vs 25")
plotMA(resLFC, ylim=c(-2,2),main="30 vs 25")

resultFPM <- as_tibble(resLFC)
resultFPM$Name = rownames(resLFC)
resSig = resultFPM %>% filter(baseMean > 10 & (padj < 0.01 & abs(log2FoldChange) > 2)) %>% 
  dplyr::select(c(Name,baseMean,log2FoldChange,lfcSE,pvalue,padj)) %>% arrange(log2FoldChange,padj)

write.csv(resSig,"reports/all_30_vs_25.csv")


resSig <- subset(resLFC, padj < 0.05)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/bbmap_all_30_vs_25.csv")
write.csv(fpm(dds),"reports/bbmap_all_FPM.csv")

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

resultFPM <- as_tibble(resLFC)
resultFPM$Name = rownames(resLFC)
resSig = resultFPM %>% filter(baseMean > 10 & (padj < 0.01 & abs(log2FoldChange) > 2)) %>% 
  dplyr::select(c(Name,baseMean,log2FoldChange,lfcSE,pvalue,padj)) %>% arrange(log2FoldChange,padj)

write.csv(resSig,"reports/all_37_vs_25.csv")

resSig <- subset(resLFC, padj < 0.05)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/all_LFC_37_vs_25.csv")

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

resultFPM <- as_tibble(resLFC)
resultFPM$Name = rownames(resLFC)
resSig = resultFPM %>% filter(baseMean > 10 & (padj < 0.01 | abs(log2FoldChange) > 2)) %>% 
  dplyr::select(c(Name,baseMean,log2FoldChange,lfcSE,pvalue,padj)) %>% arrange(log2FoldChange,padj)
write.csv(resSig,"reports/all_42_vs_25.csv")

# Get diff expressed
resSig <- subset(resLFC, padj < 0.05)
resSig <- resSig[order(resSig$pvalue),]
write.csv(resSig,"reports/all_LFC_42_vs_25.csv")

topresSig <- as.data.frame(head(resSig$log2FoldChange,50))
colnames(topresSig) = c("log2FoldChange")

p <-pheatmap(topresSig, show_rownames=TRUE,show_colnames=FALSE,cluster_rows=FALSE,
         fontsize_row = 6,fontsize_col = 7,width=5,cellwidth=10,
         cluster_cols=FALSE, main="Log(2) fold change of 42 vs 25")
ggsave("plots/42_vs_25_logfold_heatmap.pdf",p,width=16,height=16)
