setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(tidyverse)
library(janitor)
library(GenomicRanges)


samples_info <- readxl::read_xlsx("data/samples_info.xlsx") %>%
  clean_names()

cns_files <- list.files("outputs/CNVKit/", pattern = 'call.cns$', full.names = T)
cns <- lapply(cns_files, read.table, header = T) %>%
  `names<-`(gsub(".call.cns", "", list.files("outputs/CNVKit/", pattern = 'call.cns$')))

samples_info <- samples_info %>% filter(sample_id %in% names(cns))
samples_info <- samples_info %>% filter(tumor_vs_control != "Normal fat")

spont_tumor <- samples_info %>%
  filter(spon_vs_cell_line_vs_cldt == 'Spont Tumor p53PTEN') %>%
  pull(sample_id)

cell_lines <- samples_info %>%
  filter(spon_vs_cell_line_vs_cldt != 'Spont Tumor p53PTEN') %>%
  pull(sample_id)

cns_spont_tumor <- cns[spont_tumor]

cns_spont_tumor_mutational_burden <- lapply(cns_spont_tumor, function(x){
  freq <- data.frame(table(x$cn)) %>%
    filter(Var1 != 2) %>%
    `colnames<-`(c("CN", "Freq")) %>%
    mutate(CN = as.numeric(as.character(CN))) %>%
    mutate(Event = ifelse(CN < 2, "Loss", "Gain"))
}) %>% bind_rows(.id = "sample")

ggplot(cns_spont_tumor_mutational_burden, aes(x = sample, y = Freq, fill = Event)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("No. of Unbalanced Segments (Spontaneous Tumor)")

cns_spont_cellLine <- cns[cell_lines]

cns_cellLine_mutational_burden <- lapply(cns_spont_cellLine, function(x){
  freq <- data.frame(table(x$cn)) %>%
    filter(Var1 != 2) %>%
    `colnames<-`(c("CN", "Freq")) %>%
    mutate(CN = as.numeric(as.character(CN))) %>%
    mutate(Event = ifelse(CN < 2, "Loss", "Gain"))
}) %>% bind_rows(.id = "sample")

ggplot(cns_cellLine_mutational_burden, aes(x = sample, y = Freq, fill = Event)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("No. of Unbalanced Segments (Cell Line)")


## Calculating percent altered
cns_fil <- lapply(cns, function(x){x %>% filter(cn != 2)})
cns_gr <- lapply(cns_fil, GRanges)
cns_gr <- lapply(cns_gr, function(x){x$fraction  <- (width(x) / 2.7e9) ; return(x)})
percent_altered <- lapply(cns_gr, function(x){sum(x$fraction)}) %>% unlist() %>% data.frame() %>% rownames_to_column() %>% `colnames<-`(c("sample_id", "percent_altered"))
metadata <- samples_info %>%
  filter(spon_vs_cell_line_vs_cldt == 'Spont Tumor p53PTEN') %>%
  select(sample_id, tumor_vs_control, ti_ls)
percent_altered <- merge(percent_altered, metadata, by = 'sample_id')
percent_altered$condition <- paste0(percent_altered$tumor_vs_control, " (", percent_altered$ti_ls, ")")


ggplot(percent_altered, aes(x = condition, y = percent_altered)) +
  geom_boxplot() +
  geom_point() +
  xlab("") + ylab("% Genome Altered") +
  theme_bw()

percent_altered <- percent_altered[order(percent_altered$percent_altered),]
percent_altered$sample_id <- factor(percent_altered$sample_id, levels = percent_altered$sample_id)
ggplot(percent_altered, aes(x = sample_id, y = percent_altered, fill = condition)) +
  geom_bar(stat = 'identity') +
  xlab("") + ylab("% Genome Altered") +
  labs(fill = "Grade") +
  coord_flip()

percent_altered$condition <- factor(percent_altered$condition, levels = c("WDLPS (Low TILs)", "DDLPS (High TILs)", "DDLPS (Low TILs)"))
percent_altered <- percent_altered[order(percent_altered$condition),]
percent_altered$sample_id <- factor(percent_altered$sample_id, levels = percent_altered$sample_id)
ggplot(percent_altered, aes(x = sample_id, y = percent_altered, fill = condition)) +
  geom_bar(stat = 'identity') +
  xlab("") + ylab("% Genome Altered") +
  labs(fill = "Grade") +
  coord_flip()

### CREATE COMPARISONS FOR STAT_COMPARE_MEANS()
TEST.CNV = list(c("DDLPS (Low TILs)", "WDLPS (Low TILs)"), c("DDLPS (High TILs)", "DDLPS (Low TILs)"), c("DDLPS (High TILs)", "WDLPS (Low TILs)"))

ggplot(percent_altered, aes(x = condition, y = percent_altered)) +
  geom_boxplot() +
  geom_point() +
  xlab("") + ylab("% Genome Altered") +
  theme_bw()

### STATISTICS BOXPLOT FOR % GENOME ALTERED
ggplot(percent_altered, aes(x = condition, y = percent_altered)) + geom_boxplot(alpha = 0.8, fill = c("red", "dodgerblue", "blue3")) +
  stat_summary(fun = mean) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) +
  geom_point() + stat_compare_means(comparisons = TEST.CNV, label =  "p.signif") + theme_minimal() +
  labs(x = "", y = "% Genome Altered") + theme(axis.text.x = element_text(face="bold", size = 11.5), axis.title.y = element_text(face="bold", size = 12))

### CHANGING FIGURE AESTHETICS
ggplot(percent_altered, aes(x = sample_id, y = percent_altered, fill = condition)) +
  geom_bar(stat = 'identity') +
  xlab("") + ylab("% Genome Altered") +
  labs(fill = "Grade") +
  coord_flip() + theme_minimal() + scale_fill_manual(values = c("red", "dodgerblue", "blue3")) + 
  theme(axis.title.x = element_text(face = "bold", size = 11.5, vjust = -1), 
        axis.text.y = element_text(face = "bold", size = 11.5, color = "black"), 
        legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold", size = 9.25))
