setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressPackageStartupMessages({
  library(pheatmap)
  library(janitor)
  library(tidyverse)
  library(readxl)
  library(GenomicRanges)
})

samples_info <- read_xlsx("data/samples_info.xlsx") %>%
  as.data.frame() %>%
  clean_names() %>%
  filter(wes_rna_or_both != 'RNA') %>%
  mutate(spon_vs_cell_line_vs_cldt = ifelse(tumor_vs_control %in% c("Normal liver", "Normal fat"), tumor_vs_control, spon_vs_cell_line_vs_cldt)) %>%
  dplyr::select(c('sample_id', "mouse_id", "sex", "tumor_vs_control", "spon_vs_cell_line_vs_cldt"))

tumor_samples_names <- samples_info %>% filter(!tumor_vs_control %in% c("Normal_fat", "Normal_liver", "AAV8_fat")) %>% pull(sample_id)
spont_tumor <- samples_info %>% filter(!tumor_vs_control %in% c("Normal_fat", "Normal_liver", "AAV8_fat") & spon_vs_cell_line_vs_cldt == 'Spont Tumor p53PTEN') %>% arrange(desc(tumor_vs_control)) %>% pull(sample_id)
cellLine <- samples_info %>% filter(!tumor_vs_control %in% c("Normal_fat", "Normal_liver", "AAV8_fat") & spon_vs_cell_line_vs_cldt != 'Spont Tumor p53PTEN') %>% arrange(desc(tumor_vs_control)) %>% pull(sample_id)

tumor_cnr <- lapply(list.files("outputs/CNVKit/", pattern = '[0-9].cnr$', full.names = T), read.delim, header = T, sep = "\t") %>%
  `names<-`(gsub(".cnr","", list.files("outputs/CNVKit/", pattern = '[0-9].cnr$')))

pten_exon_5_logR <- lapply(tumor_samples_names, function(x){
  cnr <- GRanges(tumor_cnr[[x]])
  pten_gr <- GRanges("chr19:32777261-32777499")
  pten_overlap <- findOverlaps(cnr, pten_gr)
  pten_cnr <- cnr[queryHits(pten_overlap),]
  log2 <- weighted.mean(pten_cnr$log2, pten_cnr$weight)
  event <- ifelse(log2 < -0.25, "loss", ifelse(log2 > +0.2, "gain", "neutral"))
  return(data.frame(log2, event))
}) %>% `names<-`(tumor_samples_names) %>%
  bind_rows(.id = 'sample_id') %>%
  mutate(gene = 'Pten_Exon5') %>%
  dplyr::select(sample_id, gene, log2, event)

# pten_logR <- lapply(tumor_samples_names, function(x){
#   cnr <- GRanges(tumor_cnr[[x]])
#   pten_gr <- GRanges("chr19:32734897-32803560")
#   pten_overlap <- findOverlaps(cnr, pten_gr)
#   pten_cnr <- cnr[queryHits(pten_overlap),]
#   log2 <- weighted.mean(pten_cnr$log2, pten_cnr$weight)
#   event <- ifelse(log2 < -0.25, "loss", ifelse(log2 > +0.2, "gain", "neutral"))
#   return(data.frame(log2, event))
# }) %>% `names<-`(tumor_samples_names) %>%
#   bind_rows(.id = 'sample_id') %>%
#   mutate(gene = 'Pten') %>%
#   dplyr::select(sample_id, gene, log2, event)

trp53_logR <- lapply(tumor_samples_names, function(x){
  cnr <- GRanges(tumor_cnr[[x]])
  trp53_gr <- GRanges("chr11:69471185-69482699")
  trp53_overlap <- findOverlaps(cnr, trp53_gr)
  trp53_cnr <- cnr[queryHits(trp53_overlap),]
  log2 <- weighted.mean(trp53_cnr$log2, trp53_cnr$weight)
  event <- ifelse(log2 < -0.25, "loss", ifelse(log2 > +0.2, "gain", "neutral"))
  return(data.frame(log2, event))
}) %>% `names<-`(tumor_samples_names) %>%
  bind_rows(.id = 'sample_id') %>%
  mutate(gene = 'Trp53') %>%
  dplyr::select(sample_id, gene, log2, event)

plot_df <- rbind(pten_exon_5_logR, trp53_logR)

plot_df_spont <- plot_df %>% filter(sample_id %in% spont_tumor) %>% mutate(sample_id = factor(sample_id, levels = spont_tumor))
plot_df_spont$sample_id <- factor(plot_df_spont$sample_id, levels = c("N1307", "N1637", "N1690", "N1255", "N1369","N1314", "N1345", "N1010", "N1019", "N1062", "N1108", "N1285", "N1319", "N1692"))
log2_spont <- ggplot(plot_df_spont, aes(sample_id, gene)) +
  geom_tile(aes(fill = log2), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = 'white', high  = "red", midpoint = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), axis.title.y = element_blank(), axis.title.x = element_blank())

plot_df_celline <- plot_df %>% filter(sample_id %in% cellLine) %>% mutate(sample_id = factor(sample_id, levels = cellLine))
log2_cellline <- ggplot(plot_df_celline, aes(sample_id, gene)) +
  geom_tile(aes(fill = log2), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = 'white', high  = "red", midpoint = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), axis.title.y = element_blank(), axis.title.x = element_blank())

mdm2_logR <- lapply(tumor_samples_names, function(x){
  cnr <- GRanges(tumor_cnr[[x]])
  mdm2_gr <- GRanges("chr10:117524780-117546663")
  mdm2_overlap <- findOverlaps(cnr, mdm2_gr)
  mdm2_cnr <- cnr[queryHits(mdm2_overlap),]
  log2 <- weighted.mean(mdm2_cnr$log2, mdm2_cnr$weight)
  event <- ifelse(log2 < -0.25, "loss", ifelse(log2 > +0.2, "gain", "neutral"))
  return(data.frame(log2, event))
}) %>% `names<-`(tumor_samples_names) %>%
  bind_rows(.id = 'sample_id') %>%
  mutate(gene = 'Mdm2') %>%
  dplyr::select(sample_id, gene, log2, event)
frs2_logR <- lapply(tumor_samples_names, function(x){
  cnr <- GRanges(tumor_cnr[[x]])
  frs2_gr <- GRanges("chr10:116905332-116984415")
  frs2_overlap <- findOverlaps(cnr, frs2_gr)
  frs2_cnr <- cnr[queryHits(frs2_overlap),]
  log2 <- weighted.mean(frs2_cnr$log2, frs2_cnr$weight)
  event <- ifelse(log2 < -0.25, "loss", ifelse(log2 > +0.2, "gain", "neutral"))
  return(data.frame(log2, event))
}) %>% `names<-`(tumor_samples_names) %>%
  bind_rows(.id = 'sample_id') %>%
  mutate(gene = 'Frs2') %>%
  dplyr::select(sample_id, gene, log2, event)

hmga2_logR <- lapply(tumor_samples_names, function(x){
  cnr <- GRanges(tumor_cnr[[x]])
  hmga2_gr <- GRanges("chr10:120197180-120312374")
  hmga2_overlap <- findOverlaps(cnr, hmga2_gr)
  hmga2_cnr <- cnr[queryHits(hmga2_overlap),]
  log2 <- weighted.mean(hmga2_cnr$log2, hmga2_cnr$weight)
  event <- ifelse(log2 < -0.25, "loss", ifelse(log2 > +0.2, "gain", "neutral"))
  return(data.frame(log2, event))
}) %>% `names<-`(tumor_samples_names) %>%
  bind_rows(.id = 'sample_id') %>%
  mutate(gene = 'Hmga2') %>%
  dplyr::select(sample_id, gene, log2, event)

cdk4_logR <- lapply(tumor_samples_names, function(x){
  cnr <- GRanges(tumor_cnr[[x]])
  cdk4_gr <- GRanges("chr10:126899403-126903789")
  cdk4_overlap <- findOverlaps(cnr, cdk4_gr)
  cdk4_cnr <- cnr[queryHits(cdk4_overlap),]
  log2 <- weighted.mean(cdk4_cnr$log2, cdk4_cnr$weight)
  event <- ifelse(log2 < -0.25, "loss", ifelse(log2 > +0.2, "gain", "neutral"))
  return(data.frame(log2, event))
}) %>% `names<-`(tumor_samples_names) %>%
  bind_rows(.id = 'sample_id') %>%
  mutate(gene = 'Cdk4') %>%
  dplyr::select(sample_id, gene, log2, event)

plot_df <- rbind(mdm2_logR, frs2_logR, hmga2_logR, cdk4_logR)

plot_df_spont <- plot_df %>% filter(sample_id %in% spont_tumor) %>% mutate(sample_id = factor(sample_id, levels = spont_tumor))
log2_spont <- ggplot(plot_df_spont, aes(sample_id, gene)) +
  geom_tile(aes(fill = log2), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = 'white', high  = "red", midpoint = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), axis.title.y = element_blank(), axis.title.x = element_blank())

plot_df_celline <- plot_df %>% filter(sample_id %in% cellLine) %>% mutate(sample_id = factor(sample_id, levels = cellLine))
log2_cellline <- ggplot(plot_df_celline, aes(sample_id, gene)) +
  geom_tile(aes(fill = log2), colour = "white") +
  scale_fill_gradient2(low = "blue", mid = 'white', high  = "red", midpoint = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), axis.title.y = element_blank(), axis.title.x = element_blank())

# event <- ggplot(plot_df, aes(sample_id, gene)) +
#   geom_tile(aes(fill = factor(event)), colour = "white") +
#   scale_fill_manual(values = c('neutral' = 'white', 'gain' = 'red', 'loss' = 'darkblue')) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

