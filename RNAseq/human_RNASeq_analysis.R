# This file contains the analysis for human RNAseq data

# setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(PCAtools)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
})

dir.create("outputs", showWarnings = F, recursive = T)

# loading and processing data ------------------------------------------------------------

hallmarks_pathways <- msigdbr(species = "Homo sapiens", category = 'H', subcategory = NULL)
hallmarks_pathways <- split(x = hallmarks_pathways$ensembl_gene, f = hallmarks_pathways$gs_name)

ct_mtx <- read.table("data/human_samples_gene_count.txt", quote = "", header = T, fill = T, sep = "\t")
gene_symbols <- data.frame(gene_symbol = ct_mtx$gene_name) %>% `rownames<-`(ct_mtx$gene_id)

ct_mtx <- ct_mtx %>%
  select(-c("gene_name", "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family")) %>%
  column_to_rownames('gene_id')

metadata <- read.table("data/human_rnaseq_samples_info.txt", sep = '\t', header = T) %>%
  column_to_rownames('sample_id')
metadata <- metadata[colnames(ct_mtx),]
metadata$condition <- gsub("NF", "Normal_fat", metadata$condition)
metadata$condition <- gsub("WD", "WDLPS", metadata$condition)
metadata$condition <- gsub("DD", "DDLPS", metadata$condition)
metadata$recurrence <- gsub("NR", "No_Rec", metadata$recurrence)
metadata$recurrence <- gsub("R", "Rec", metadata$recurrence)

# DESeq object all samples ------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = ct_mtx,
                              colData = metadata,
                              rowData = gene_symbols,
                              design = ~ 1)
dds$condition <- factor(dds$condition, levels = c('Normal_fat', 'WDLPS', 'DDLPS'))

saveRDS(dds, "outputs/all_human_samples_dds.rds")

# WDLPS vs fat ------------------------------------------------------------

wdlps_idx <- dds$condition %in% c('Normal_fat', 'WDLPS')
dds_wdlps <- dds[,wdlps_idx]
dds_wdlps$condition <- factor(dds_wdlps$condition, levels = c('Normal_fat', 'WDLPS'))
design(dds_wdlps) <- model.matrix(~ condition, data = colData(dds_wdlps))
dds_wdlps <- DESeq(dds_wdlps)

# > DGE analysis ----------------------------------------------------------

wdlps_v_fat <- lfcShrink(dds_wdlps, coef = 'conditionWDLPS', type = 'ashr') %>%
  data.frame() %>%
  merge(gene_symbols, by = 0) %>%
  filter(!is.na(padj)) %>%
  rename(gene_id = Row.names) %>%
  mutate(status = ifelse(log2FoldChange > 0.5 & padj < 0.05, 1,
                         ifelse(log2FoldChange < -0.5 & padj < 0.05, -1, 0)))

write.csv(wdlps_v_fat, "outputs/human_wdlps_v_fat_dge.csv", row.names = F, quote = F)

# > GSEA ------------------------------------------------------------------

wdlps_v_fat_rank <- wdlps_v_fat %>%
  pull(log2FoldChange) %>%
  `names<-`(wdlps_v_fat$gene_id) %>%
  sort()

wdlps_v_fat_gsea <- fgsea(pathways = hallmarks_pathways,
                          stats = wdlps_v_fat_rank)

saveRDS(wdlps_v_fat_gsea, file = "outputs/human_wdlps_v_fat_gsea.rds")
write.csv(wdlps_v_fat_gsea[,-8], file = "outputs/human_wdlps_v_fat_gsea.csv", quote = FALSE, row.names = F)

# DDLPS vs fat ------------------------------------------------------------

ddlps_idx <- dds$condition %in% c("Normal_fat","DDLPS")
dds_ddlps <- dds[,ddlps_idx]
dds_ddlps$condition <- factor(dds_ddlps$condition, levels = c("Normal_fat","DDLPS"))
design(dds_ddlps) <- model.matrix(~ condition, data = colData(dds_ddlps))
dds_ddlps <- DESeq(dds_ddlps)

# > DGE analysis ----------------------------------------------------------

ddlps_v_fat <- lfcShrink(dds_ddlps, coef = 'conditionDDLPS', type = 'ashr') %>%
  data.frame() %>%
  merge(gene_symbols, by = 0) %>%
  filter(!is.na(padj)) %>%
  rename(gene_id = Row.names) %>%
  mutate(status = ifelse(log2FoldChange > 0.5 & padj < 0.05, 1,
                         ifelse(log2FoldChange < -0.5 & padj < 0.05, -1, 0)))

write.csv(ddlps_v_fat, "outputs/human_ddlps_v_fat_dge.csv", row.names = F, quote = F)

# > GSEA ------------------------------------------------------------------

ddlps_v_fat_rank <- ddlps_v_fat %>%
  pull(log2FoldChange) %>%
  `names<-`(ddlps_v_fat$gene_id) %>%
  sort()

ddlps_v_fat_gsea <- fgsea(pathways = hallmarks_pathways,
                          stats = ddlps_v_fat_rank)

saveRDS(ddlps_v_fat_gsea, file = "outputs/human_ddlps_v_fat_gsea.rds")
write.csv(ddlps_v_fat_gsea[,-8], file = "outputs/human_ddlps_v_fat_gsea.csv", quote = FALSE, row.names = F)


