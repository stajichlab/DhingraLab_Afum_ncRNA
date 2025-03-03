#!/usr/bin/env Rscript

library(tidyverse)
library(fs)
library(cowplot)

unpaired_folder <- "results/ncRNA_bbmap_rpkm"
fwdfiles <- fs::dir_ls(unpaired_folder,regexp = "\\.single_1\\.rpkm\\.txt$")
revfiles <- fs::dir_ls(unpaired_folder,regexp = "\\.single_2\\.rpkm\\.txt$")

fwddata <- fwdfiles %>% map_dfr(read_tsv,.id = "source", skip=4) %>% 
  mutate(Name = str_split_i(`#Name`,' ',1),
         sourcedata = str_split_i(str_split_i(source,'/',3),'\\.',1)) %>% 
#         sourcedata = paste0(str_split_i(str_split_i(source,'/',3),'\\.',1),"_FWD")) %>% 
  select(c(sourcedata,Length,Bases,Coverage,Name,Reads,RPKM,Frags,FPKM))

ncRNA_fwd <- fwddata %>% filter(grepl('Afu-',Name) ) %>% mutate(READ="FWD")

revdata <- revfiles %>% map_dfr(read_tsv,.id = "source", skip=4) %>% 
  mutate(Name = str_split_i(`#Name`,' ',1),
#         sourcedata = paste0(str_split_i(str_split_i(source,'/',3),'\\.',1),"_REV")) %>% 
sourcedata = str_split_i(str_split_i(source,'/',3),'\\.',1)) %>% 
  select(c(sourcedata,Length,Bases,Coverage,Name,Reads,RPKM,Frags,FPKM))

bylib_len_fwd <- fwddata %>% group_by(sourcedata) %>% summarize(sum(Reads)) %>% 
  mutate(FWDCT=`sum(Reads)`) %>% select(c(sourcedata,FWDCT))
bylib_len_rev <- revdata %>% group_by(sourcedata) %>% summarize(sum(Reads)) %>%
  mutate(REVCT=`sum(Reads)`) %>% select(c(sourcedata,REVCT))

libcount <- bylib_len_fwd %>% left_join(bylib_len_rev) %>% mutate(TOTAL=FWDCT + REVCT)

# fix at 
bygene_fwd <- ncRNA_fwd %>% select(c(sourcedata,Length,Name,Reads)) %>% 
  pivot_wider(names_from = sourcedata, values_from=Reads)
bygene_rev <- ncRNA_rev %>% select(c(sourcedata,Length,Name,Reads)) %>% 
  pivot_wider(names_from = sourcedata, values_from=Reads)

rownames(bygene_fwd) <-bygene_fwd$Name
rownames(bygene_rev) <-bygene_rev$Name

bind_rows(bygene_fwd, bygene_rev ) %>% select( -c(Length)) %>%
  group_by(Name) %>% 
  summarise_all(sum) %>% 

1000 * (3411/1488) / 18.5
ncRNA_rev <- revdata %>% filter(grepl('Afu-',Name) ) %>% mutate(READ="REV") %>% 
  filter(grepl('WT_',sourcedata))

ncRNA <- bind_rows(ncRNA_fwd,ncRNA_rev) %>% 
  filter(grepl('WT_',sourcedata))

#ncRNA %>% 
#  # evaluate following calls for each value in the rowname column
#  group_by(sourcedata,Name) %>% 
#  # add all non-grouping variables
#  summarise_all(sum)

p0 <- ggplot(data=ncRNA %>% filter(Name != "Afu-254"), 
             aes(x=sourcedata, y=FPKM, fill=Name)) +
  geom_bar(stat="identity", width=0.75) + 
  theme_cowplot(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(rows=vars(READ))

p0

p0A <- ggplot(data=ncRNA %>% filter(Name == "Afu-182"), 
             aes(x=sourcedata, y=FPKM, fill=Name)) +
  geom_bar(stat="identity",width=0.75) + 
  theme_cowplot(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(rows=vars(READ))

p0A

ggsave("plots/Afu-182_fwd_rev.pdf",p0A,width=8,height=8)
#p1 <- ggplot(data=ncRNA_fwd %>% filter(Name == "Afu-182"), aes(x=sourcedata, y=FPKM, fill=Name)) +
p1 <- ggplot(data=ncRNA_fwd %>% filter(Name != "Afu-254"), aes(x=sourcedata, y=FPKM, fill=Name)) +
  geom_bar(stat="identity", width=1) + 
  theme_cowplot(12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

p2 <- ggplot(data=ncRNA_rev %>% filter(Name != "Afu-254"), aes(x=sourcedata, y=FPKM, fill=Name)) +
  geom_bar(stat="identity", width=.75) + theme_cowplot(12) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  

plot_grid(p1, p2, labels = c('FWD', 'REV'), label_size = 12)

