#!/usr/bin/env Rscript

library(dplyr)
library(tidyverse)
library(ggplot2)
library(UpSetR)

# afu-182 deletion comparison
folder="reports"
Temp42_WT_Del <- read_csv(file.path(folder,"Temp42_WT_vs_Del.csv"), col_names=TRUE)
Temp42names = Temp42_WT_Del %>% rename(GeneID="...1") %>% dplyr::select(GeneID) %>% add_column(Temperature=42)

Temp37_WT_Del <- read_csv(file.path(folder,"Temp37_WT_vs_Del.csv"), col_names=TRUE)
Temp37names = Temp37_WT_Del %>% rename(GeneID="...1")  %>% dplyr::select(GeneID) %>% add_column(Temperature=37)

Temp25_WT_Del <- read_csv(file.path(folder,"Temp25_WT_vs_Del.csv"), col_names=TRUE)
Temp25names = Temp25_WT_Del %>% rename(GeneID="...1") %>% dplyr::select(GeneID) %>% add_column(Temperature=25)

#WT_af182Del = bind_rows(Temp42names,Temp37names,Temp25names)

# without descriptor
Temp42_WT_Del <- read_csv(file.path(folder,"Temp42_WT_vs_Del.csv"), col_names=TRUE)
Temp42names = Temp42_WT_Del %>% rename(GeneID="...1") %>% dplyr::select(GeneID) 

Temp37_WT_Del <- read_csv(file.path(folder,"Temp37_WT_vs_Del.csv"), col_names=TRUE)
Temp37names = Temp37_WT_Del %>% rename(GeneID="...1")  %>% dplyr::select(GeneID) 

Temp25_WT_Del <- read_csv(file.path(folder,"Temp25_WT_vs_Del.csv"), col_names=TRUE)
Temp25names = Temp25_WT_Del %>% rename(GeneID="...1") %>% dplyr::select(GeneID) 



WT_af182Del = list(Temp42 = Temp42names$GeneID,
                   Temp37 = Temp37names$GeneID,
                   Temp25 = Temp25names$GeneID
                   )

upset(fromList(WT_af182Del), order.by = "freq",  
      mainbar.y.label = "Genes",
      sets.x.label    = "Significant Genes Per Comparison" )

pdf("plots/UpSet_WT_vs_afu182delta.pdf")
upset(fromList(WT_af182Del), order.by = "freq",  
      mainbar.y.label = "Genes",
      sets.x.label    = "Significant Genes Per Comparison")

dev.off()


# WT temperature comparison
WT_25_v_42 <- read_csv(file.path(folder,"WT_25_vs_42.csv"), col_names=TRUE)
WT25_42names = WT_25_v_42 %>% rename(GeneID="...1") %>% dplyr::select(GeneID)

WT_25_v_37 <- read_csv(file.path(folder,"WT_25_vs_37.csv"), col_names=TRUE)
WT25_37names = WT_25_v_37 %>% rename(GeneID="...1") %>% dplyr::select(GeneID)

WT_37_v_42 <- read_csv(file.path(folder,"WT_37_vs_42.csv"), col_names=TRUE)
WT37_42names = WT_37_v_42 %>% rename(GeneID="...1") %>% dplyr::select(GeneID)


WT_TempCompare = list(Temp25_42 = WT25_42names$GeneID,
                   Temp25_37 = WT25_37names$GeneID,
                   Temp37_42 = WT37_42names$GeneID)

upset(fromList(WT_TempCompare), order.by = "freq",  
      mainbar.y.label = "Genes",
      sets.x.label    = "Significant Genes Per Comparison" )

pdf("plots/UpSet_WT_Temp_compare.pdf")
upset(fromList(WT_TempCompare), order.by = "freq",  
      mainbar.y.label = "Genes",
      sets.x.label    = "Significant Genes Per Comparison")
dev.off()

