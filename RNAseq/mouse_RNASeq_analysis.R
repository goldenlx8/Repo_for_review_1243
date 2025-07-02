# This file contains the mouse RNAseq analysis

# setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(PCAtools)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(EnhancedVolcano)
})

dir.create("outputs", showWarnings = F, recursive = T)

# loading data ------------------------------------------------------------

hallmarks_pathways <- msigdbr(species = "Mus musculus", category = 'H', subcategory = NULL)
hallmarks_pathways <- split(x = hallmarks_pathways$ensembl_gene, f = hallmarks_pathways$gs_name)

um_core_samples <- read.table("data/UM_core_samples.txt", header = T, quote = "", sep = "\t")
novogene_samples <- read.table("data/Novogene_samples_UM_core_Processing.txt", header = T, quote = "", sep = "\t")

ct_mtx <- merge(um_core_samples, novogene_samples, by = c("gene_id", "entrezgene_id", "external_gene_name", "description"))
gene_symbols <- data.frame(gene_symbol = ct_mtx$external_gene_name) %>% `rownames<-`(ct_mtx$gene_id)

ct_mtx <- ct_mtx %>%
  select(-c("entrezgene_id", "external_gene_name", "description")) %>%
  column_to_rownames('gene_id') %>%
  `colnames<-`(gsub(pattern = "^X", replacement = "", x = colnames(.))) %>%
  `colnames<-`(gsub(pattern = "\\.", replacement = "\\-", x = colnames(.)))

metadata <- read.table("data/mouse_rnaseq_samples_info.txt", sep = '\t', header = T) %>%
  select(sample_id, sex, mouse_id, condition, tissue, TILs, source) %>%
  column_to_rownames('sample_id')

# DESeq object all samples ------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = ct_mtx,
                              colData = metadata,
                              rowData = gene_symbols,
                              design = ~ 1)
saveRDS(dds, "outputs/all_mouse_samples_dds.rds")
write.csv(merge(assay(vst(dds)),gene_symbols, by = 0),
          "outputs/mouse_normalized_count_mtx.csv",
          row.names = F)

## Restrict the analysis to normal fat, WDLPS, DDLPS, spotaneous tumors
## Remove N1059 samples as ...
idx <- dds$condition %in% c("Normal_fat","WDLPS","DDLPS") &
  dds$tissue %in% c("Spont Tumor p53PTEN", "Normal_fat") &
  colnames(dds) != 'N1059'
dds <- dds[,idx]
dds$condition <- factor(dds$condition, levels = c('Normal_fat', 'WDLPS', 'DDLPS'))
saveRDS(dds, "outputs/subset_mouse_dds.rds")

# WDLPS vs fat ------------------------------------------------------------

wdlps_idx <- dds$condition %in% c("Normal_fat","WDLPS") &
  dds$tissue %in% c("Spont Tumor p53PTEN", "Normal_fat")
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

write.csv(wdlps_v_fat, "outputs/mouse_wdlps_v_fat_dge.csv", row.names = F, quote = F)

mouse_ann <- data.frame(gene = wdlps_v_fat$gene_id,
                        gene_symbol = wdlps_v_fat$gene_symbol,
                        pval_adj = wdlps_v_fat$padj)

EnhancedVolcano(toptable = wdlps_v_fat,
                lab = wdlps_v_fat$gene_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1,
                pCutoff = 0.05,
                drawConnectors = T, labSize = 3,
                title = 'WDLPS vs. Normal Fat', subtitle = NULL)

### REMOVE LABELS FROM VOLCANO PLOT & MOVE LEGEND
EnhancedVolcano(toptable = wdlps_v_fat,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1,
                pCutoff = 0.05,
                drawConnectors = T,
                legendPosition = "right",
                title = 'WDLPS vs. Normal Fat', subtitle = NULL) + 
  theme(plot.title = element_text(hjust = 0.5))

# > GSEA ------------------------------------------------------------------

wdlps_v_fat_rank <- wdlps_v_fat %>%
  pull(log2FoldChange) %>%
  `names<-`(wdlps_v_fat$gene_id) %>%
  sort()

wdlps_v_fat_gsea <- fgsea(pathways = hallmarks_pathways,
                          stats = wdlps_v_fat_rank)

saveRDS(wdlps_v_fat_gsea, file = "outputs/mouse_wdlps_v_fat_gsea.rds")
write.csv(wdlps_v_fat_gsea[,-8],
          file = "outputs/mouse_wdlps_v_fat_gsea.csv", quote = FALSE, row.names = F)


# DDLPS vs fat ------------------------------------------------------------

ddlps_idx <- dds$condition %in% c("Normal_fat","DDLPS") &
  dds$tissue %in% c("Spont Tumor p53PTEN", "Normal_fat")
dds_ddlps <- dds[,ddlps_idx]
dds_ddlps$condition <- factor(dds_ddlps$condition, levels = c('Normal_fat', 'DDLPS'))
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

write.csv(ddlps_v_fat, "outputs/mouse_ddlps_v_fat_dge.csv", row.names = F, quote = F)

EnhancedVolcano(toptable = ddlps_v_fat,
                lab = ddlps_v_fat$gene_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1,
                pCutoff = 0.05,
                drawConnectors = T, labSize = 3,
                title = 'DDLPS vs. Normal Fat', subtitle = NULL)


### REMOVE LABELS FROM VOLCANO PLOT & MOVE LEGEND
EnhancedVolcano(toptable = ddlps_v_fat,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1,
                pCutoff = 0.05,
                drawConnectors = T,
                legendPosition = "right",
                title = 'DDLPS vs. Normal Fat', subtitle = NULL) +
  theme(plot.title = element_text(hjust = 0.5))


### CHANGE VP TO ONLY HIGHLIGHTED GENES 
goi = c('Cdk4', 'Cdk6', 'Cebpa', 'Mdm2', 'Hmga2')
gl = ifelse(ddlps_v_fat$gene_symbol %in% goi, 'cyan', 'black')

GOI.VP = EnhancedVolcano(toptable = ddlps_v_fat,
               lab = ddlps_v_fat$gene_symbol,
               selectLab = c('Cdk4', 'Cdk6', 'Cebpa', 'Mdm2', 'Hmga2'),
               x = 'log2FoldChange',
               y = 'padj',
               FCcutoff = 1,
               pCutoff = 0.05,
               drawConnectors = T, labSize = 4,
               boxedLabels = T, legendPosition = "right",
               labFace = "bold",
               arrowheads = FALSE, widthConnectors = 1.25,
               directionConnectors = "both", colConnectors = "black",
               title = 'DDLPS vs. Normal Fat', subtitle = NULL) + 
  theme(plot.title = element_text(hjust = 0.5))

GOI.VP + geom_point(data = ddlps_v_fat[ddlps_v_fat$gene_symbol %in% goi, ], 
                    aes(x = log2FoldChange, y = -log10(padj)), 
                    color = "cyan", size = 2.5, shape = 1, stroke = 2)

with(ddlps_v_fat[ddlps_v_fat$gene_symbol %in% goi, ],
                 points(log2FoldChange, -log10(padj), cex = 2, lwd = 2, col = 'black', pch = 1))

# > GSEA ------------------------------------------------------------------

ddlps_v_fat_rank <- ddlps_v_fat %>%
  pull(log2FoldChange) %>%
  `names<-`(ddlps_v_fat$gene_id) %>%
  sort()

ddlps_v_fat_gsea <- fgsea(pathways = hallmarks_pathways,
                          stats = ddlps_v_fat_rank)

saveRDS(ddlps_v_fat_gsea, file = "outputs/mouse_ddlps_v_fat_gsea.rds")
write.csv(ddlps_v_fat_gsea[,-8],
          file = "outputs/mouse_ddlps_v_fat_gsea.csv", quote = FALSE, row.names = F)

