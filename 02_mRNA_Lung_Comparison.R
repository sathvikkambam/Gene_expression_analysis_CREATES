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

#Batch effect needs to be removed
#DESeq2??

#Loads tsv files of normal lung cells and appends vertically

files1 <- list.files("data/Normal_Lung", pattern="*.tsv", full.names=TRUE)
col_name = c("gene_id","gene_name","gene_type","unstranded","stranded_first","stranded_second","tpm_unstranded","fpkm_unstranded","fpkm_uq_unstranded")
lung_normal <- read_tsv(files1, comment="#", col_names = col_name, skip=6)

#Isolates columns of interest and labels each gene by the sample it was taken from
unstranded_norm <- lung_normal |> select("gene_id", "gene_name", "unstranded") |> mutate(
  sample = ((row_number()-1) %/% 60660) + 1
)

#Turns it into a table of unstranded RNA counts
raw_norm <- unstranded_norm |> select(!gene_name) |> pivot_wider(names_from=gene_id, values_from = unstranded)
raw_norm <- raw_norm |> select(!sample)
norm_matrix_raw <- as.matrix(t(raw_norm))
colnames(norm_matrix_raw) <- paste("N", 1:311 |> lapply(toString), sep = "")


#Isolates only the CD151 data from normal lung cells
cd151_normal <- lung_normal |> filter(gene_name == "CD151")
cd151_normal <- cd151_normal |> select(unstranded, tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded)
cd151_normal <- cd151_normal |> mutate(
  status = "normal"
)

#Extracts data from tsv files of lung cancer cells
files2 <- list.files("data/Cancer_Lung", pattern="*.tsv", full.names=TRUE)
col_name = c("gene_id","gene_name","gene_type","unstranded","stranded_first","stranded_second","tpm_unstranded","fpkm_unstranded","fpkm_uq_unstranded")
lung_cancer <- read_tsv(files2, comment="#", col_names = col_name, skip=6)

#Labels all genes with sample number
unstranded_cancer <- lung_cancer |> select("gene_id", "gene_name", "unstranded") |> mutate(
  sample = ((row_number()-1) %/% 60660) + 1
)

#Gets raw count matrix
raw_cancer <- unstranded_cancer |> select(!gene_name) |> pivot_wider(names_from=gene_id, values_from = unstranded)
raw_cancer <- raw_cancer |> select(!sample)
norm_matrix_cancer <- as.matrix(t(raw_cancer))
colnames(norm_matrix_cancer) <- paste("C", 1:338 |> lapply(toString), sep = "")

#Combines cancer and normal data
norm_matrix <- cbind(norm_matrix_raw, norm_matrix_cancer)

#Changes column data
coldata <- data.frame(1:649)
rownames(coldata) <- colnames(norm_matrix)
colnames(coldata) <- "condition"
coldata <- coldata |> mutate(
  condition = if_else(row_number() <= 311, "normal", "cancer")
)

coldata$condition <- factor(coldata$condition)

#Parallel Cores for faster running
#Need to set parallel = TRUE
library("BiocParallel")
register(SnowParam(4))

#Uses DESeq2 method to analyze the matrix data
dds <- DESeqDataSetFromMatrix(countData = norm_matrix,
                              colData = coldata,
                              design = ~ condition)
#Adding gene name data
featureData <- data.frame(gene_name=lung_normal[1:60660, "gene_name"])
mcols(dds) <- DataFrame(mcols(dds), featureData)


#Prefilters to speed up method
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#Sets normal patients as the control
dds$condition <- relevel(dds$condition, ref = "normal")

#Runs all the analysis
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds)

resLFC <- lfcShrink(dds, coef="condition_cancer_vs_normal", type="apeglm", parallel = TRUE)
plotMA(resLFC, ylim=c(-10,10))

#To get gene name, run:
#(unstranded_norm |> filter(gene_name=="CD109"))[1, ]$gene_id
cPlot <- plotCounts(dds, 
                    gene="ENSG00000177697.19", 
                    intgroup="condition", 
                    returnData = TRUE
                    )
cPlot |> ggplot(aes(x=condition, y=count, color=condition)) + geom_boxplot() + 
  geom_point(position=position_jitter(w=0.1,h=0), alpha = 0.2)

#Gets all CD151 data from cancer cells
cd151_cancer <- lung_cancer |> filter(gene_name == "CD151")
cd151_cancer <- cd151_cancer |> select(unstranded, tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded)
cd151_cancer <- cd151_cancer |> mutate(
  status = "cancer"
)

#Joins cd151 data toether for plotting
cd151 <- bind_rows(cd151_normal, cd151_cancer)
cd151 |> ggplot(aes(y=unstranded, x=status, color=status)) + geom_boxplot()


cd151$status <- factor(cd151$status, levels=c("normal", "cancer"))
cd151 |> ggplot(aes(y=tpm_unstranded, x=status, color=status)) + geom_boxplot()


#Need to remove batch effect
#Calculates the probability of statistical difference between gene expression
t.test(cd151_cancer$tpm_unstranded, cd151_normal$tpm_unstranded)
