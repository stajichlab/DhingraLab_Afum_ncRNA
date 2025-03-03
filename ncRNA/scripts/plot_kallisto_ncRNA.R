#!/usr/bin/env Rscript

library(tidyverse)
library(vroom)
library(fs)
library(cowplot)

folder <- "results/kallisto_alltranscripts"
experiments <- fs::dir_ls(folder,recurse=TRUE,glob="*.tsv")

countdata <- experiments %>% map_dfr(read_tsv,.id = "source") %>% 
  mutate(sourceup = str_split_i(source,'/',4),
         kmer = str_split_i(str_split_i(source,'/',3),'_',2)) %>% 
           select(-c(source))

ncRNA <- countdata %>% filter(grepl('Afu-',target_id) ) %>% filter(kmer == 19)

p <- ggplot(data=ncRNA, aes(x=sourceup, y=log(tpm)/log(10), fill=target_id)) +
  geom_bar(stat="identity", width=.75) + facet_wrap(~kmer) + 
  ggtitle("Af293 kallisto RNASeq expression") + ylab("log10(tpm)") + theme_cowplot(12) + theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) 

ggsave("ncRNA_kallisto.kmer_plot.pdf",p,width=16,height=12)

p2 <- ggplot(data=ncRNA %>% filter(target_id == 'Afu-182' & kmer == 19 ), aes(x=sourceup, y=tpm, fill=target_id)) +
  geom_bar(stat="identity", width=.75) + theme_cowplot(12) + theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) 


ggsave("ncRNA_kallisto.Afu-182.kmer_19.pdf",p2,width=16,height=12) + theme_cowplot(12)
