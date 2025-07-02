title: "3D PCA HUMAN & MOUSE LPS v NORMAL FAT"
author: "EMMA KENNA (ekenna@umich.edu)"

# LOADING PACKAGES --------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(DESeq2)
library(tidyverse)
library(Glimma)
library(msigdbr)
library(fgsea)
library(PCAtools)
library(DT)
library(edgeR)
library(plotly)
library(stats)
library(biomaRt)
library(babelgene)
library(ashr)
library(ggplot2)

# LOADING DATA ------------------------------------------------------------
### HUMAN DATA 

hallmarks_pathways = msigdbr(species = "Homo sapiens", category = 'H', subcategory = NULL)
hallmarks_pathways = split(x = hallmarks_pathways$ensembl_gene, f = hallmarks_pathways$gs_name)

ct_mtx = read.table("data/human_samples_gene_count.txt", quote = "", header = T, fill = T, sep = "\t")
gene_symbols = data.frame(gene_symbol = ct_mtx$gene_name) %>% `rownames<-`(ct_mtx$gene_id)

ct_mtx = ct_mtx %>%
  dplyr::select(-c("gene_name", "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family")) %>%
  column_to_rownames('gene_id')

metadata = read.table("data/human_rnaseq_samples_info_species.txt", sep = '\t', header = T) %>%
  column_to_rownames('sample_id')
metadata = metadata %>% dplyr::select(-c("recurrence", "TILS"))
metadata = metadata[colnames(ct_mtx),]

### MOUSE DATA

m_hallmarks_pathways = msigdbr(species = "Mus musculus", category = 'H', subcategory = NULL)
m_hallmarks_pathways = split(x = m_hallmarks_pathways$ensembl_gene, f = m_hallmarks_pathways$gs_name)

um_core_samples = read.table("data/UM_core_samples.txt", header = T, quote = "", sep = "\t")
novogene_samples = read.table("data/Novogene_samples_UM_core_Processing.txt", header = T, quote = "", sep = "\t")

m_ct_mtx = merge(um_core_samples, novogene_samples, by = c("gene_id", "entrezgene_id", "external_gene_name", "description"))
m_gene_symbols = data.frame(gene_symbol = m_ct_mtx$external_gene_name) %>% `rownames<-`(m_ct_mtx$gene_id)

m_ct_mtx = m_ct_mtx %>%
  dplyr::select(-c("entrezgene_id", "external_gene_name", "description")) %>%
  column_to_rownames('gene_id') %>%
  `colnames<-`(gsub(pattern = "^X", replacement = "", x = colnames(.))) %>%
  `colnames<-`(gsub(pattern = "\\.", replacement = "\\-", x = colnames(.)))

m_metadata = read.table("data/mouse_rnaseq_samples_info_species.txt", sep = '\t', header = T) %>%
  column_to_rownames('sample_id')
m_metadata = m_metadata %>% dplyr::select(condition, species)
m_metadata = m_metadata[colnames(m_ct_mtx),]

# DESeq OBJECT ALL SAMPLES ------------------------------------------------

### HUMAN 

H.dds = DESeqDataSetFromMatrix(countData = ct_mtx,
                               colData = metadata,
                               rowData = gene_symbols,
                               design = ~ 1)
H.dds$condition = factor(H.dds$condition, levels = c('Normal_Fat', 'WDLPS', 'DDLPS'))
saveRDS(H.dds, "outputs/all_human_samples_dds_species.rds")

### MOUSE 

M.dds = DESeqDataSetFromMatrix(countData = m_ct_mtx,
                               colData = m_metadata,
                               rowData = m_gene_symbols,
                               design = ~ 1)
idx = M.dds$condition %in% c("Normal_Fat","WDLPS","DDLPS") &
  colnames(M.dds) != 'N1059'& colnames(M.dds) != 'N1011' & colnames(M.dds) != 'N1018' &
  colnames(M.dds) != 'N1076' & colnames(M.dds) != 'N1128'
M.dds = M.dds[,idx]
M.dds$condition = factor(M.dds$condition, levels = c('Normal_Fat', 'WDLPS', 'DDLPS'))
saveRDS(M.dds, "outputs/subset_mouse_dds_species.rds")

# MERGING HUMAN & MOUSE RDS FILES -----------------------------------------

### LOAD IN THE TWO RDS FILES IF NOT ALREADY LOADED IN

H.DF = readRDS('outputs/all_human_samples_dds_species.rds')
M.DF = readRDS('outputs/subset_mouse_dds_species.rds')

meta1 = as.data.frame(colData(H.DF)) %>% rownames_to_column(., var = "sample_id")
meta2 = as.data.frame(colData(M.DF)) %>% rownames_to_column(., var = "sample_id")

## ADDING SPECIES TO CONDITION COLUMN - NEED FOR 4 COLOR VISUALIZATION LATER
meta1$condition = as.character(meta1$condition)
meta1$condition[meta1$condition == "DDLPS"] = "Human DDLPS"
meta1$condition[meta1$condition == "Normal_Fat"] = "Human Normal Fat"

meta2$condition = as.character(meta2$condition)
meta2$condition[meta2$condition == "DDLPS"] = "Mouse DDLPS"
meta2$condition[meta2$condition == "Normal_Fat"] = "Mouse Normal Fat"

merged.meta = full_join(meta1, meta2, by = c("condition", "species", "sample_id")) 

### PULL COUNT DATA

m.counts = counts(M.DF)
h.counts = counts(H.DF)

m.counts.sub = m.counts[-c(1,2), ]
mouse.genes = rownames(m.counts.sub)
orth = babelgene::orthologs(genes = mouse.genes, species = "mouse", human = FALSE)

common.genes = intersect(orth$human_ensembl, rownames(h.counts))

c.indx = orth$human_ensembl %in% common.genes
orth.common = orth[c.indx, ]
orth.common.sub = orth.common %>% dplyr::select(ensembl, human_ensembl)

indx = rownames(h.counts) %in% common.genes
common.h.counts = h.counts[indx, ]

m.counts.es = as.data.frame(m.counts.sub) %>% rownames_to_column(., var = "ensembl")
orth.m.counts = m.counts.es %>% full_join(orth.common.sub, by = join_by(ensembl))
common.m.counts = orth.m.counts %>% semi_join(orth.common.sub, by = join_by(ensembl)) %>% 
  dplyr::select(-ensembl)
common.m.counts = column_to_rownames(common.m.counts, var = "human_ensembl")

### REMOVING DUPLICATE ENSEMBL ROWS
DR.common.m.counts = common.m.counts %>% rowwise() %>% 
  mutate(avg_expr = mean(c_across(-human_ensembl), na.rm = TRUE)) %>% 
  ungroup() %>% group_by(human_ensembl) %>% 
  slice_max(avg_expr, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% dplyr::select(-avg_expr)

filt.common.m.counts = column_to_rownames(DR.common.m.counts, var = "human_ensembl")
m.common.m.counts = as.matrix(filt.common.m.counts)

combined.counts = cbind(m.common.m.counts, common.h.counts)

### SAVE THE LARGE MATRIX FOR LATER USE - IF NEEDED
saveRDS(combined.counts, "outputs/combined_counts.rds")

# HUMAN & MOUSE PCA PLOTTING ----------------------------------------------

### NORMALIZE & LOG TRANSFORM COUNT MATRIX - DIFFERENT STRATEGY USED
### ADJUSTING META DATA FILE FOR PCA ANALYSIS


merged.meta = merged.meta %>% column_to_rownames(., var = "sample_id")
merged.meta = merged.meta[colnames(combined.counts),]

### USE DESeq2 TO TRANSFORM COUNTS FOR PCA

HM.DDS = DESeqDataSetFromMatrix(countData = combined.counts, 
                                colData = merged.meta, 
                                design = ~ species + condition)

### NORMALIZE COUNT MATRIX USING VST

HM.VST = assay(vst(HM.DDS, blind = TRUE))

### PERFORMING PCA ON THE COMBINED COUNTS

HM.RV = rowVars(HM.VST)
HM.SELECT = order(HM.RV, decreasing = TRUE)[seq_len(min(500, length(HM.DDS)))]
HM.PCA = pca(HM.VST[HM.SELECT, ], metadata = colData(HM.DDS))

## EXTRACTING DATA FOR PLOTTING

HM.PCA.DF = as.data.frame(HM.PCA$rotated)
HM.PCA.DF$sample = rownames(HM.PCA.DF)
merged.meta$sample = rownames(merged.meta)
ANN.HM.PCA = merge(HM.PCA.DF, merged.meta, by = "sample")
NW.ANN.PCA = ANN.HM.PCA[ANN.HM.PCA$condition != "WDLPS",]

### 3D PCA - CONDITION 4 COLORS 

four.color = plot_ly(data = NW.ANN.PCA,
                     x = ~PC2, y = ~PC3, z = ~PC4,
                     color = ~condition,
                     colors = c("deeppink", "cornflowerblue", "purple", "black"),
                     text = ~sample,
                     type = "scatter3d",
                     mode = "markers",
                     marker = list(size = 5)) %>%
  layout(title = "3D PCA Plot: Histological Subtypes across Species")

four.color = four.color %>% layout(legend = list(font = list(size = 16), x = 0.8, y = 0.5,
                                                 xanchor = "left", yanchor = "middle"))
four.color
