# setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(PCAtools)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(mMCPcounter)
  library(pheatmap)
  library(escape)
  library(nichenetr)
  library(GSVA)
  library(cowplot)
  library("AnnotationDbi")
  library("org.Mm.eg.db")
})

dir.create("outputs", showWarnings = F, recursive = T)

# loading data ------------------------------------------------------------

subset_mouse_dds <- readRDS("outputs/subset_mouse_dds.rds")
subset_mouse_dds <- estimateSizeFactors(subset_mouse_dds)
ct_mtx <- counts(subset_mouse_dds)
normalized_ct_mtx <- log(counts(subset_mouse_dds, normalized = TRUE) + 1)
sample_names <- read.table("data/mouse_rnaseq_samples_info.txt", sep = "\t", header = T, row.names = 'sample_id') %>%
  dplyr::select("sample_id.old.") %>%
  `colnames<-`("old_sample_id")
metadata <- as.data.frame(colData(subset_mouse_dds)) %>%
  merge(sample_names, by = 0) %>%
  column_to_rownames("Row.names")

# GSVA ------------------------------------------------------------

data("mMCPcounter_signatures_GCRm38",
     envir = sys.frame(sys.nframe()),
     package = "mMCPcounter")
write.csv(mMCPcounter_signatures_GCRm38, "outputs/mMCPcounter_signatures_GCRm38.csv")

mMCPcounter_signatures_GCRm38 <- mMCPcounter_signatures_GCRm38 %>%
  dplyr::select(Denomination, ENSEMBL.ID)

GeneSets <- split(mMCPcounter_signatures_GCRm38,
                  mMCPcounter_signatures_GCRm38$Denomination)
GeneSets <- lapply(GeneSets, function(x){pull(x,'ENSEMBL.ID')})

gaussian_gsva_enrich <- gsva(gsvaParam(log(cpm(counts(subset_mouse_dds)) + 1), GeneSets, kcdf = 'Gaussian'))
color_palette <- colorRampPalette(c('blue','white','red'))
pheatmap(gaussian_gsva_enrich[unique(mMCPcounter_signatures_GCRm38$Denomination),order(metadata$TILs)],
         color = color_palette(100),
         annotation_col = metadata[,c('TILs', 'condition')],
         cluster_rows = T, cluster_cols = F,
         main = 'GSVA Enrichemt (using mMCP Signatures)', gaps_col = c(3, 13),
         labels_col = metadata$old_sample_id[order(metadata$TILs)])
