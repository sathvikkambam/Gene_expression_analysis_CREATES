library(tidyverse)
library(ggplot2)
library(DESeq2)
library("Biobase")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("matrixStats")
library("lattice")
library("gtools")
library("dplyr")


files <- list.files("data/Lung_Cancer_Data_Test", pattern="*.tsv", full.names=TRUE)
col_name = c("gene_id","gene_name","gene_type","unstranded","stranded_first","stranded_second","tpm_unstranded","fpkm_unstranded","fpkm_uq_unstranded")
lung_cancer <- read_tsv(files, comment="#", col_names = col_name, skip=6)

cd151_exp <- lung_cancer |> filter(gene_name=="CD151")
cd151_exp |> ggplot(aes(x=tpm_unstranded)) + geom_density(color="red")

lung_normal <- read_tsv("data/Normal_Lung_Test/gene_tpm_2022-06-06_v10_lung.gct", skip=2)
cd151_norm <- lung_normal |> filter(Description == "CD151")

filesN <- list.files("data/Norm_LC_2", pattern="*.tsv", full.names=TRUE)
lung_norm_2 <- read_tsv(filesN, comment="#", col_names = col_name, skip=6)
cd151_norm_2 <- lung_norm_2 |> filter(gene_name=="CD151")

