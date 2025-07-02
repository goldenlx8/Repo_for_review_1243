# This analysis contains the reanalysis of the DDLPS v FAT MS GENE SIGNATURE on HUMAN DDLPS & NF

# SET DIRECTORY & LOAD PACKAGES -------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(DESeq2)
library(tidyverse)
library(Glimma)
library(PCAtools)
library(edgeR)
library(fgsea)
library(msigdbr)
library(EnhancedVolcano)

# LOADING DATA ------------------------------------------------------------

hallmarks_pathways = msigdbr(species = "Homo sapiens", category = 'H', subcategory = NULL)
hallmarks_pathways = split(x = hallmarks_pathways$ensembl_gene, f = hallmarks_pathways$gs_name)

ct_mtx = read.table("data/human_samples_gene_count.txt", quote = "", header = T, fill = T, sep = "\t")
gene_symbols = data.frame(gene_symbol = ct_mtx$gene_name) %>% `rownames<-`(ct_mtx$gene_id)

ct_mtx = ct_mtx %>%
  dplyr::select(-c("gene_name", "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family")) %>%
  column_to_rownames('gene_id')

metadata = read.table("data/human_rnaseq_samples_info.txt", sep = '\t', header = T) %>%
  column_to_rownames('sample_id')
metadata = metadata[colnames(ct_mtx),]
metadata$condition = gsub("NF", "Normal_fat", metadata$condition)
metadata$condition = gsub("WD", "WDLPS", metadata$condition)
metadata$condition = gsub("DD", "DDLPS", metadata$condition)
metadata$recurrence = gsub("NR", "No_Rec", metadata$recurrence)
metadata$recurrence = gsub("R", "Rec", metadata$recurrence)

# DESeq OBJECT ALL SAMPLES ------------------------------------------------
dds = DESeqDataSetFromMatrix(countData = ct_mtx,
                              colData = metadata,
                              rowData = gene_symbols,
                              design = ~ 1)
dds$condition = factor(dds$condition, levels = c('Normal_fat', 'WDLPS', 'DDLPS'))

saveRDS(dds, "outputs/all_human_samples_dds.rds")
write.csv(merge(assay(vst(dds)),gene_symbols, by = 0),
          "outputs/human_normalized_count_mtx.csv",
          row.names = F)

# DDLPS vs FAT ------------------------------------------------------------

ddlps_idx = dds$condition %in% c("Normal_fat","DDLPS")
dds_ddlps = dds[,ddlps_idx]
dds_ddlps$condition = factor(dds_ddlps$condition, levels = c("Normal_fat","DDLPS"))
design(dds_ddlps) = model.matrix(~ condition, data = colData(dds_ddlps))
dds_ddlps = DESeq(dds_ddlps)

# LOADING MORE PACKAGES ---------------------------------------------------

library(BiocManager)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathview)
library(gage)
library(gageData)
### REQUIRED PACKAGES
library(matrixStats)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(GSVA)
library(DESeq2)
library(tidyverse)
library(PCAtools)
library(edgeR)
library(fgsea)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(nichenetr)

# MOUSE DDLPS v FAT GENE SIGNATURE --> HUMAN DDLPS & NF -------------------

### LOADING DATA 
#mouse_wdlps_v_fat_dge = read.csv("outputs/mouse_wdlps_v_fat_dge.csv")
mouse_ddlps_v_fat_dge = read.csv("outputs/mouse_ddlps_v_fat_dge.csv")
#human_wdlps_v_fat_dge = read.csv("outputs/human_wdlps_v_fat_dge.csv")
#human_ddlps_v_fat_dge = read.csv("outputs/human_ddlps_v_fat_dge.csv")


### MS DDLPS GENE SIGNATURE - CONVERT MS GENE SYMBOLS TO HUMAN TO USE ON HUMAN SAMPLES

MS.DDLPS.V.FAT.UP = mouse_ddlps_v_fat_dge %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>%
  pull(gene_symbol) %>% convert_mouse_to_human_symbols()

MS.DDLPS.V.FAT.UP = MS.DDLPS.V.FAT.UP[!is.na(MS.DDLPS.V.FAT.UP)]

# CONVERT THE HUMANIZED MOUSE UP GENE LIST TO HUMAN ENSEMBL ID'S 
MS.DDLPS.V.FAT.UP = mapIds(org.Hs.eg.db, keys = MS.DDLPS.V.FAT.UP, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
MS.DDLPS.V.FAT.UP = MS.DDLPS.V.FAT.UP[!is.na(MS.DDLPS.V.FAT.UP)]

MS.DDLPS.V.FAT.DN = mouse_ddlps_v_fat_dge %>%
  filter(padj < 0.05 & log2FoldChange < -2) %>%
  pull(gene_symbol) %>% convert_mouse_to_human_symbols()

MS.DDLPS.V.FAT.DN = MS.DDLPS.V.FAT.DN[!is.na(MS.DDLPS.V.FAT.DN)]

# CONVERT THE HUMANIZED MOUSE DOWN GENE LIST TO HUMAN ENSEMBL ID'S 
MS.DDLPS.V.FAT.DN = mapIds(org.Hs.eg.db, keys = MS.DDLPS.V.FAT.DN, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
MS.DDLPS.V.FAT.DN = MS.DDLPS.V.FAT.DN[!is.na(MS.DDLPS.V.FAT.DN)]

GS.DDLPS.MS = list(MS.DDLPS.V.FAT.UP = MS.DDLPS.V.FAT.UP,
                   MS.DDLPS.V.FAT.DN = MS.DDLPS.V.FAT.DN) 

### REMOVING THE WDLPS SAMPLES FROM THE METADATA AFTER DDS_DDLPS IS CREATED - NEED FOR LATER VISUALIZATION
metadata = metadata[metadata$condition != "WDLPS", ]

### GSVA ENRICHMENT 

metadata.n = metadata %>% subset(select = -c(recurrence, TILS))
metadata.n$condition = gsub("Normal_fat", "Normal Fat", metadata.n$condition)
metadata.n = metadata.n %>% rename(Condition = condition)

GAUSSIAN.GSVA.ENRICH = gsva(gsvaParam(log(cpm(counts(dds_ddlps)) + 1), GS.DDLPS.MS, kcdf = 'Gaussian'))
COLOR.PALETTE = colorRampPalette(c('blue','white','deeppink'))

# USING SUBSETTED METADATA - CLEANER PLOT
ANN.COLOR = list(Condition = c(DDLPS = "purple2", `Normal Fat` = "seagreen1"))

pheatmap(GAUSSIAN.GSVA.ENRICH[,order(metadata.n$Condition)], annotation_col = metadata.n, 
         cluster_rows = F, cluster_cols = F, main = "MS DDLPS Gene Signature GSVA Enrichment in Human DDLPS & Normal Fat", 
         gaps_col = 8, color = COLOR.PALETTE(10), annotation_colors = ANN.COLOR)
