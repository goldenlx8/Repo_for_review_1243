setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Importing libraries
suppressPackageStartupMessages({
  library(karyoploteR)
  library(regioneR)
  library(zoo)
  library(janitor)
  library(readxl)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(tidyverse)
  library(fields)
})

# Setting up output dir
outputDir <- "outputs/Pten_Trp53_plots"
dir.create(outputDir, showWarnings = F, recursive = T)

# Importing samples info
samples_info <- read_xlsx("data/samples_info.xlsx") %>%
  as.data.frame() %>%
  clean_names() %>%
  filter(wes_rna_or_both != 'RNA') %>%
  mutate(spon_vs_cell_line_vs_cldt = ifelse(tumor_vs_control %in% c("Normal_fat", "Normal_liver", "AAV8_fat"), tumor_vs_control, spon_vs_cell_line_vs_cldt)) %>%
  dplyr::select(c('sample_id', "mouse_id", "sex", "tumor_vs_control", "spon_vs_cell_line_vs_cldt"))

spont_tumor <- samples_info %>% filter(spon_vs_cell_line_vs_cldt == "Spont Tumor p53PTEN") %>% arrange(desc(tumor_vs_control)) %>% pull(sample_id)
cellLine <- samples_info %>% filter(spon_vs_cell_line_vs_cldt %in% c("Cell-line 1011", "Cell-line 1018", "CLDT 1011", "CLDT 1018")) %>% pull(sample_id)
tumor_samples_names <- samples_info %>% filter(!tumor_vs_control %in% c("Normal_fat", "Normal_liver", "AAV8_fat")) %>% pull(sample_id)
color_palette <- colorRampPalette(c("blue", "red"))(10)

# Importing coverage bed files
coverage <- lapply(samples_info$sample_id, function(x){
  message(x)
  gr <- GRanges(read.delim(paste0("outputs/bed_files/", x,".bed"), sep = "\t", header = F) %>%
                  dplyr::select(V1, V2, V3) %>%
                  `colnames<-`(c("chromosome", "start", "end")))
}) %>% `names<-`(samples_info$sample_id)

# Importing CNVKit cnr files
cnr <- lapply(samples_info %>% filter(!tumor_vs_control %in% c("Normal_fat", "Normal_liver", "AAV8_fat")) %>% pull(sample_id), function(x){
  gr <- GRanges(read.delim(paste0("outputs/CNVKit_revised_run/", x, ".cnr"), sep = "\t", header = T))
}) %>% `names<-`(samples_info %>% filter(!tumor_vs_control %in% c("Normal_fat", "Normal_liver", "AAV8_fat")) %>% pull(sample_id))

target_bed <- GRanges(read.table("outputs/CNVKit_revised_run/S0276129_Regions_mm39.target.bed") %>%
                        `colnames<-`(c("chromosome", "start", "end", "gene")))

cnr <- lapply(cnr, function(x){
  GenomicRanges::merge(x, target_bed)
})
all_cnr <- bind_rows(lapply(cnr, data.frame), .id = 'sample_id')

# Trp53 -------------------------------------------------------------------

# > Spontanous Tumor ------------------------------------------------------

# > Setup -----------------------------------------------------------------

spont_tumor_cnr <- all_cnr %>% filter(sample_id %in% spont_tumor) %>% GRanges()
spont_tumor_trp53_overlaps <- findOverlaps(spont_tumor_cnr, GRanges("chr11:69468307-69485577"))
spont_tumor_trp53_cnr <- spont_tumor_cnr[queryHits(spont_tumor_trp53_overlaps),]
spont_tumor_ymin = min(spont_tumor_trp53_cnr$log2)
spont_tumor_ymax = max(spont_tumor_trp53_cnr$log2)

# > Plot Trp53 Track ------------------------------------------------------

Trp53 <- plotKaryotype(genome = 'mm39',
                       plot.type = 2,
                       chromosomes = c("chr11"),
                       zoom = GRanges("chr11:69468307-69485577"))
kpAddBaseNumbers(Trp53, tick.dist = 2000,  tick.len = 10, add.units = T, cex = 0.5,
                 minor.ticks = T, minor.tick.len = 5, minor.tick.dist = 100)
Trp53_genes.data <- makeGenesDataFromTxDb(txdb=TxDb.Mmusculus.UCSC.mm39.knownGene, karyoplot=Trp53)
Trp53_genes.data <- addGeneNames(Trp53_genes.data)
Trp53_genes.data <- mergeTranscripts(Trp53_genes.data)

total.tracks <- 42
kpPlotGenes(Trp53, Trp53_genes.data, gene.name.position = "left", data.panel = 2, r0 = 0.3, r1=0.7)

# > N1307 -------------------------------------------------------------------

at <- autotrack(current.track = 1, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 2, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1307']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1307"), cex = 0.5)
at <- autotrack(current.track = 3, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1307']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1637 -------------------------------------------------------------------

at <- autotrack(current.track = 4, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 5, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1637']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1637"), cex = 0.5)
at <- autotrack(current.track = 6, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1637']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1690 -------------------------------------------------------------------

at <- autotrack(current.track = 7, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 8, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1690']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1690"), cex = 0.5)
at <- autotrack(current.track = 9, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1690']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1010 -------------------------------------------------------------------

at <- autotrack(current.track = 10, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 11, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1010']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1010"), cex = 0.5)
at <- autotrack(current.track = 12, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1010']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1019 -------------------------------------------------------------------

at <- autotrack(current.track = 13, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 14, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1019']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1019"), cex = 0.5)
at <- autotrack(current.track = 15, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1019']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1062 -------------------------------------------------------------------

at <- autotrack(current.track = 16, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 17, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1062']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1062"), cex = 0.5)
at <- autotrack(current.track = 18, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1062']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1108 -------------------------------------------------------------------

at <- autotrack(current.track = 19, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 20, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1108']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1108"), cex = 0.5)
at <- autotrack(current.track = 21, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1108']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1255 -------------------------------------------------------------------

at <- autotrack(current.track = 22, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 23, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1255']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1255"), cex = 0.5)
at <- autotrack(current.track = 24, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1255']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1285 -------------------------------------------------------------------

at <- autotrack(current.track = 25, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 26, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1285']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1285"), cex = 0.5)
at <- autotrack(current.track = 27, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1285']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1314 -------------------------------------------------------------------

at <- autotrack(current.track = 28, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 29, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1314']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1314"), cex = 0.5)
at <- autotrack(current.track = 30, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1314']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)


# > N1319 -------------------------------------------------------------------

at <- autotrack(current.track = 31, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 32, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1319']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1319"), cex = 0.5)
at <- autotrack(current.track = 33, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1319']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)


# > N1345 -------------------------------------------------------------------

at <- autotrack(current.track = 34, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 35, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1345']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1345"), cex = 0.5)
at <- autotrack(current.track = 36, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1345']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1369 -------------------------------------------------------------------

at <- autotrack(current.track = 37, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 38, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1369']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1369"), cex = 0.5)
at <- autotrack(current.track = 39, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1369']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1692 -------------------------------------------------------------------

at <- autotrack(current.track = 40, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 41, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1692']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1692"), cex = 0.5)
at <- autotrack(current.track = 42, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_trp53_cnr[spont_tumor_trp53_cnr$sample_id == 'N1692']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > Legend ------------------------------------------------------------------

image.plot(legend.only = TRUE,
           zlim = c(spont_tumor_ymin, spont_tumor_ymax),
           col = color_palette,
           legend.args = list(text = "log2 ratio", side = 4, line = 2.5, cex = 1),
           smallplot = c(0.89, 0.9, 0.9, 0.98))


# CellLine ----------------------------------------------------------------

# > Setup -----------------------------------------------------------------

cellLine_cnr <- all_cnr %>% filter(sample_id %in% cellLine) %>% GRanges()
cellLine_trp53_overlaps <- findOverlaps(cellLine_cnr, GRanges("chr11:69468307-69485577"))
cellLine_trp53_cnr <- cellLine_cnr[queryHits(cellLine_trp53_overlaps),]
cellLine_ymin = min(cellLine_trp53_cnr$log2)
cellLine_ymax = max(cellLine_trp53_cnr$log2)

# > Plot Trp53 Track ------------------------------------------------------

Trp53 <- plotKaryotype(genome = 'mm39',
                       plot.type = 2,
                       chromosomes = c("chr11"),
                       zoom = GRanges("chr11:69468307-69485577"))
kpAddBaseNumbers(Trp53, tick.dist = 2000,  tick.len = 10, add.units = T, cex = 0.5,
                 minor.ticks = T, minor.tick.len = 5, minor.tick.dist = 100)
Trp53_genes.data <- makeGenesDataFromTxDb(txdb=TxDb.Mmusculus.UCSC.mm39.knownGene, karyoplot=Trp53)
Trp53_genes.data <- addGeneNames(Trp53_genes.data)
Trp53_genes.data <- mergeTranscripts(Trp53_genes.data)

total.tracks <- 12
kpPlotGenes(Trp53, Trp53_genes.data, gene.name.position = "left", data.panel = 2, r0 = 0.3, r1=0.7)

# > N1011 -------------------------------------------------------------------

at <- autotrack(current.track = 1, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 2, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1011']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1011"), cex = 0.5)
at <- autotrack(current.track = 3, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- cellLine_trp53_cnr[cellLine_trp53_cnr$sample_id == 'N1011']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = cellLine_ymin, ymax = cellLine_ymax)

# > N1018 -------------------------------------------------------------------

at <- autotrack(current.track = 4, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 5, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1018']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1018"), cex = 0.5)
at <- autotrack(current.track = 6, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- cellLine_trp53_cnr[cellLine_trp53_cnr$sample_id == 'N1018']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = cellLine_ymin, ymax = cellLine_ymax)

# > N1076 -------------------------------------------------------------------

at <- autotrack(current.track = 7, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 8, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1076']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1076"), cex = 0.5)
at <- autotrack(current.track = 9, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- cellLine_trp53_cnr[cellLine_trp53_cnr$sample_id == 'N1076']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = cellLine_ymin, ymax = cellLine_ymax)

# > N1128 -------------------------------------------------------------------

at <- autotrack(current.track = 10, total.tracks = total.tracks, margin = 0.2)
kpRect(Trp53, chr = 'chr11',x0 = 69468307, x1 = 69485577,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 11, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(Trp53, data= coverage[['N1128']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(Trp53, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1128"), cex = 0.5)
at <- autotrack(current.track = 12, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- cellLine_trp53_cnr[cellLine_trp53_cnr$sample_id == 'N1128']
kpHeatmap(Trp53,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = cellLine_ymin, ymax = cellLine_ymax)

# > Legend ------------------------------------------------------------------

image.plot(legend.only = TRUE,
           zlim = c(cellLine_ymin, cellLine_ymax),
           col = color_palette,
           legend.args = list(text = "log2 ratio", side = 4, line = 2.5, cex = 1),
           smallplot = c(0.89, 0.9, 0.9, 0.98))

# Pten --------------------------------------------------------------------
# > Setup -----------------------------------------------------------------

spont_tumor_cnr <- all_cnr %>% filter(sample_id %in% spont_tumor) %>% GRanges()
spont_tumor_pten_overlaps <- findOverlaps(spont_tumor_cnr, GRanges("chr19:32717731-32820726"))
spont_tumor_pten_cnr <- spont_tumor_cnr[queryHits(spont_tumor_pten_overlaps),]
pten_cnr <- spont_tumor_cnr[queryHits(spont_tumor_pten_overlaps),]
spont_tumor_ymin = min(spont_tumor_pten_cnr$log2)
spont_tumor_ymax = max(spont_tumor_pten_cnr$log2)

# > Spontanous Tumor ------------------------------------------------------
# > Plot Pten Track ------------------------------------------------------

pten <- plotKaryotype(genome = 'mm39',
                      plot.type = 2,
                      chromosomes = c("chr19"),
                      zoom = GRanges("chr19:32717731-32820726"))
kpAddBaseNumbers(pten, tick.dist = 10000,  tick.len = 10, add.units = T, cex = 0.5,
                 minor.ticks = T, minor.tick.len = 5, minor.tick.dist = 1000)
pten_genes.data <- makeGenesDataFromTxDb(txdb=TxDb.Mmusculus.UCSC.mm39.knownGene, karyoplot=pten)
pten_genes.data <- addGeneNames(pten_genes.data)
pten_genes.data <- mergeTranscripts(pten_genes.data)

total.tracks <- 42
kpPlotGenes(pten, pten_genes.data, gene.name.position = "left", data.panel = 2, r0 = 0.3, r1=0.7)

# > N1307 -------------------------------------------------------------------

at <- autotrack(current.track = 1, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 2, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1307']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1307"), cex = 0.5)
at <- autotrack(current.track = 3, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1307']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1637 -------------------------------------------------------------------

at <- autotrack(current.track = 4, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 5, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1637']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1637"), cex = 0.5)
at <- autotrack(current.track = 6, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1637']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1690 -------------------------------------------------------------------

at <- autotrack(current.track = 7, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 8, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1690']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1690"), cex = 0.5)
at <- autotrack(current.track = 9, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1690']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1010 -------------------------------------------------------------------

at <- autotrack(current.track = 10, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 11, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1010']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1010"), cex = 0.5)
at <- autotrack(current.track = 12, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1010']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1019 -------------------------------------------------------------------

at <- autotrack(current.track = 13, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 14, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1019']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1019"), cex = 0.5)
at <- autotrack(current.track = 15, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1019']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1062 -------------------------------------------------------------------

at <- autotrack(current.track = 16, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 17, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1062']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1062"), cex = 0.5)
at <- autotrack(current.track = 18, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1062']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1108 -------------------------------------------------------------------

at <- autotrack(current.track = 19, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 20, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1108']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1108"), cex = 0.5)
at <- autotrack(current.track = 21, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1108']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1255 -------------------------------------------------------------------

at <- autotrack(current.track = 22, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 23, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1255']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1255"), cex = 0.5)
at <- autotrack(current.track = 24, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1255']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1285 -------------------------------------------------------------------

at <- autotrack(current.track = 25, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 26, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1285']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1285"), cex = 0.5)
at <- autotrack(current.track = 27, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1285']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1314 -------------------------------------------------------------------

at <- autotrack(current.track = 28, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 29, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1314']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1314"), cex = 0.5)
at <- autotrack(current.track = 30, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1314']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)


# > N1319 -------------------------------------------------------------------

at <- autotrack(current.track = 31, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 32, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1319']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1319"), cex = 0.5)
at <- autotrack(current.track = 33, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1319']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)


# > N1345 -------------------------------------------------------------------

at <- autotrack(current.track = 34, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 35, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1345']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1345"), cex = 0.5)
at <- autotrack(current.track = 36, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1345']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1369 -------------------------------------------------------------------

at <- autotrack(current.track = 37, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 38, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1369']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1369"), cex = 0.5)
at <- autotrack(current.track = 39, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1369']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > N1692 -------------------------------------------------------------------

at <- autotrack(current.track = 40, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 41, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1692']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1692"), cex = 0.5)
at <- autotrack(current.track = 42, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- spont_tumor_pten_cnr[spont_tumor_pten_cnr$sample_id == 'N1692']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = spont_tumor_ymin, ymax = spont_tumor_ymax)

# > Legend ------------------------------------------------------------------

image.plot(legend.only = TRUE,
           zlim = c(spont_tumor_ymin, spont_tumor_ymax),
           col = color_palette,
           legend.args = list(text = "log2 ratio", side = 4, line = 2.5, cex = 1),
           smallplot = c(0.89, 0.9, 0.9, 0.98))


# CellLine ----------------------------------------------------------------

# > Setup -----------------------------------------------------------------

cellLine_cnr <- all_cnr %>% filter(sample_id %in% cellLine) %>% GRanges()
cellLine_pten_overlaps <- findOverlaps(cellLine_cnr, GRanges("chr19:32717731-32820726"))
cellLine_pten_cnr <- cellLine_cnr[queryHits(cellLine_pten_overlaps),]
cellLine_ymin = min(cellLine_pten_cnr$log2)
cellLine_ymax = max(cellLine_pten_cnr$log2)

# > Plot pten Track ------------------------------------------------------

pten <- plotKaryotype(genome = 'mm39',
                      plot.type = 2,
                      chromosomes = c("chr11"),
                      zoom = GRanges("chr19:32717731-32820726"))
kpAddBaseNumbers(pten, tick.dist = 10000,  tick.len = 10, add.units = T, cex = 0.5,
                 minor.ticks = T, minor.tick.len = 5, minor.tick.dist = 1000)
pten_genes.data <- makeGenesDataFromTxDb(txdb=TxDb.Mmusculus.UCSC.mm39.knownGene, karyoplot=pten)
pten_genes.data <- addGeneNames(pten_genes.data)
pten_genes.data <- mergeTranscripts(pten_genes.data)

total.tracks <- 12
kpPlotGenes(pten, pten_genes.data, gene.name.position = "left", data.panel = 2, r0 = 0.3, r1=0.7)

# > N1011 -------------------------------------------------------------------

at <- autotrack(current.track = 1, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 2, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1011']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1011"), cex = 0.5)
at <- autotrack(current.track = 3, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- cellLine_pten_cnr[cellLine_pten_cnr$sample_id == 'N1011']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = cellLine_ymin, ymax = cellLine_ymax)

# > N1018 -------------------------------------------------------------------

at <- autotrack(current.track = 4, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 5, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1018']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1018"), cex = 0.5)
at <- autotrack(current.track = 6, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- cellLine_pten_cnr[cellLine_pten_cnr$sample_id == 'N1018']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = cellLine_ymin, ymax = cellLine_ymax)

# > N1076 -------------------------------------------------------------------

at <- autotrack(current.track = 7, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 8, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1076']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1076"), cex = 0.5)
at <- autotrack(current.track = 9, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- cellLine_pten_cnr[cellLine_pten_cnr$sample_id == 'N1076']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = cellLine_ymin, ymax = cellLine_ymax)

# > N1128 -------------------------------------------------------------------

at <- autotrack(current.track = 10, total.tracks = total.tracks, margin = 0.2)
kpRect(pten, chr = 'chr19',x0 = 32717731, x1 = 32820726,
       y0 = 1, y1 = 1.1, r0 = at$r0, r1=at$r1,
       data.panel = 1, col = 'darkgrey', border = 'darkgrey')
at <- autotrack(current.track = 11, total.tracks = total.tracks, margin = 0.2)
kpPlotCoverage(pten, data= coverage[['N1128']], r0 = at$r0, r1=at$r1, data.panel = 1)
kpAddLabels(pten, data.panel = 1, r0 = at$r0, r1=at$r1, labels = c("N1128"), cex = 0.5)
at <- autotrack(current.track = 12, total.tracks = total.tracks, margin = 0.2)
sample_cnr_data <- cellLine_pten_cnr[cellLine_pten_cnr$sample_id == 'N1128']
kpHeatmap(pten,
          data = sample_cnr_data,
          r0 = at$r0, r1=at$r1,
          y = sample_cnr_data$log2, data.panel = 1,
          colors = color_palette,
          ymin = cellLine_ymin, ymax = cellLine_ymax)

# > Legend ------------------------------------------------------------------

image.plot(legend.only = TRUE,
           zlim = c(cellLine_ymin, cellLine_ymax),
           col = color_palette,
           legend.args = list(text = "log2 ratio", side = 4, line = 2.5, cex = 1),
           smallplot = c(0.89, 0.9, 0.9, 0.98))




