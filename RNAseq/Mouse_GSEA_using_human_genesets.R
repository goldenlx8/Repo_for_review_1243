# This file contains the analysis for the mouse GSEA plots from human gene sets

# setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressPackageStartupMessages({
library(fgsea)
library(tidyverse)
library(nichenetr)
library(cowplot)
})

# Using Human DE signature ------------------------------------------------

## Loading data
mouse_wdlps_v_fat_dge <- read.csv("outputs/mouse_wdlps_v_fat_dge.csv")
mouse_ddlps_v_fat_dge <- read.csv("outputs/mouse_ddlps_v_fat_dge.csv")
human_wdlps_v_fat_dge <- read.csv("outputs/human_wdlps_v_fat_dge.csv")
human_ddlps_v_fat_dge <- read.csv("outputs/human_ddlps_v_fat_dge.csv")

## Human WDLPS GSEA on mouse data
human_wdlps_v_fat_up <- human_wdlps_v_fat_dge %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(gene_symbol) %>% convert_human_to_mouse_symbols()
human_wdlps_v_fat_up <- human_wdlps_v_fat_up[!is.na(human_wdlps_v_fat_up)]
human_wdlps_v_fat_dn <- human_wdlps_v_fat_dge %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  pull(gene_symbol) %>% convert_human_to_mouse_symbols()
human_wdlps_v_fat_dn <- human_wdlps_v_fat_dn[!is.na(human_wdlps_v_fat_dn)]

human_wdlps_pathways <- list(human_wdlps_v_fat_up = human_wdlps_v_fat_up,
                             human_wdlps_v_fat_dn = human_wdlps_v_fat_dn)
wdlps_v_fat_rank <- mouse_wdlps_v_fat_dge %>%
  pull(log2FoldChange) %>%
  `names<-`(mouse_wdlps_v_fat_dge$gene_symbol) %>%
  sort()

wdlps_v_fat_gsea <- fgsea(pathways = human_wdlps_pathways,
                          stats = wdlps_v_fat_rank)

wdlps_up <- plotEnrichment(human_wdlps_pathways[["human_wdlps_v_fat_up"]],
                           wdlps_v_fat_rank) + ggtitle('Human WDLPS vs Fat Upreg Geneset (NES = 1.8, padj = 0.01)')
wdlps_dn <- plotEnrichment(human_wdlps_pathways[["human_wdlps_v_fat_dn"]],
                           wdlps_v_fat_rank) + ggtitle('Human WDLPS vs Fat Downreg Geneset (NES = 1.7, padj = 0.01)')
plot_grid(wdlps_up, wdlps_dn, nrow = 2)

## Human DDLPS GSEA on mouse data
human_ddlps_v_fat_up <- human_ddlps_v_fat_dge %>%
  filter(padj < 0.05 & log2FoldChange > 2) %>%
  pull(gene_symbol) %>% convert_human_to_mouse_symbols()
human_ddlps_v_fat_up <- human_ddlps_v_fat_up[!is.na(human_ddlps_v_fat_up)]
human_ddlps_v_fat_dn <- human_ddlps_v_fat_dge %>%
  filter(padj < 0.05 & log2FoldChange < -2) %>%
  pull(gene_symbol) %>% convert_human_to_mouse_symbols()
human_ddlps_v_fat_dn <- human_ddlps_v_fat_dn[!is.na(human_ddlps_v_fat_dn)]

human_ddlps_pathways <- list(human_ddlps_v_fat_up = human_ddlps_v_fat_up,
                             human_ddlps_v_fat_dn = human_ddlps_v_fat_dn)
ddlps_v_fat_rank <- mouse_ddlps_v_fat_dge %>%
  pull(log2FoldChange) %>%
  `names<-`(mouse_ddlps_v_fat_dge$gene_symbol) %>%
  sort()

ddlps_v_fat_gsea <- fgsea(pathways = human_ddlps_pathways,
                          stats = ddlps_v_fat_rank)

ddlps_up <- plotEnrichment(human_ddlps_pathways[["human_ddlps_v_fat_up"]],
                           ddlps_v_fat_rank) + ggtitle('Human DDLPS vs Fat Upreg Geneset (NES = 2.5, padj = 2e-50)')
ddlps_dn <- plotEnrichment(human_ddlps_pathways[["human_ddlps_v_fat_dn"]],
                           ddlps_v_fat_rank) + ggtitle('Human DDLPS vs Fat Downreg Geneset (NES = -2.5, padj = 9.5e-35)')
plot_grid(ddlps_up, ddlps_dn, nrow = 2)

# Comparing GSEA scores ---------------------------------------------------

mouse_wdlps_v_fat_gsea <- read.csv("outputs/mouse_wdlps_v_fat_gsea.csv")
human_wdlps_v_fat_gsea <- read.csv('outputs/human_wdlps_v_fat_gsea.csv')

mouse_NES <- mouse_wdlps_v_fat_gsea %>% dplyr::select(pathway, NES) %>% rename(Mouse_NES = NES)
human_NES <- human_wdlps_v_fat_gsea %>% dplyr::select(pathway, NES) %>% rename(Human_NES = NES)
merged_NES <- merge(mouse_NES, human_NES, by = "pathway")

ggplot(merged_NES, aes(x = Mouse_NES, y = Human_NES, text = pathway)) +
geom_point() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_label_repel(aes(label = pathway, size = 2), nudge_x = 0.1, nudge_y = 0.1, size = 2)

