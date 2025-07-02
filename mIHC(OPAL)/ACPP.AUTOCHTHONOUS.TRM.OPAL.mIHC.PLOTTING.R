# title: "MS OPAL SPONTANEOUS TRM CELL % ANALYSIS"
# author: "EMMA KENNA (ekenna@umich.edu)"

### LIBRARIES
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(huxtable)
library(viridis)

### LOADING DATA
TRM.SPONT.LPS.CP = read.csv(file = '2024.01.23.SUM.CELL.PERCENTS.TRM.csv', sep = ",", header = TRUE)
TEST = list(c("DDLPSlow", "WDLPS"), c("DDLPShigh", "DDLPSlow"), c("DDLPShigh", "WDLPS"))
TRM.SPONT.LPS.CP$HISTOLOGY = factor(TRM.SPONT.LPS.CP$HISTOLOGY, levels = c("WDLPS", "DDLPSlow", "DDLPShigh"))

### CD4+ [CD3+/CD8a-/PD1-] PLOTTING & STATISTICS 
TRM.SPONT.LPS.CP %>% select(CD3..CD8..PD1., HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = CD3..CD8..PD1.)) + stat_summary(fun = mean, geom = "bar", alpha = 0.7, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, label = 'p.signif') + labs(x = "", y = "% CD3+CD8a-PD1- / Total Cells", face = "bold", title = "CD4+ T Cells") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), axis.title.y = element_text(face="bold", size = 15)) 

### CD8+ [CD3+/CD8a+/PD1-] PLOTTING & STATISTICS 
TRM.SPONT.LPS.CP %>% select(CD3..CD8..PD1..2, HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = CD3..CD8..PD1..2)) + stat_summary(fun = mean, geom = "bar", alpha = 0.6, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, label = 'p.signif') + labs(x = "", y = "% CD3+CD8a+PD1- / Total Cells", title = "CD8+ T Cells") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), axis.title.y = element_text(face="bold", size = 15)) 

### PLOTTING WHERE PD1+ OUT OF ALL CD4+ OR CD8+
### LOADING DATA
TRM.SPONT.LPS.CP.RC = read.csv(file = '24.01.23.SUM.CELL.PERCENTS.TRM.RECALC.csv', sep = ",", header = TRUE)
TEST = list(c("DDLPSlow", "WDLPS"), c("DDLPShigh", "DDLPSlow"), c("DDLPShigh", "WDLPS"))
TRM.SPONT.LPS.CP.RC$HISTOLOGY = factor(TRM.SPONT.LPS.CP.RC$HISTOLOGY, levels = c("WDLPS", "DDLPSlow", "DDLPShigh"))

### CD4+PD1+ [CD3+/CD8a-/PD1+] OUT OF CD4+ & STATISTICS 
TRM.SPONT.LPS.CP.RC %>% select(CD3..CD8..PD1., HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = CD3..CD8..PD1.)) + stat_summary(fun = mean, geom = "bar", alpha = 0.8, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, label = "p.signif", method.args = list(exact = FALSE)) + labs(x = "", y = "% CD3+CD8a-PD1+ / CD3+CD8a- Cells", title = "CD4+ PD1+ T Cells") + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), axis.title.y = element_text(face="bold", size = 15)) 

### CD8+PD1+ [CD3+/CD8a+/PD1+] OUT OF CD8+ & STATISTICS 
TRM.SPONT.LPS.CP.RC %>% select(CD3..CD8..PD1..1, HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = CD3..CD8..PD1..1)) + stat_summary(fun = mean, geom = "bar", alpha = 0.8, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, method.args = list(exact = FALSE)) + labs(x = "", y = "% CD3+CD8a+PD1+ / CD3+CD8a+ Cells", title = "CD8+ PD1+ T Cells") + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), axis.title.y = element_text(face="bold", size = 15)) 

# TO SEE THE NUMERICAL P-VALUES REMOVE THE label = "p.signif" argument from the stat_compare_means_function 

### CD8+ PHENOTYPE STACKED BAR PLOT
### CD8+ ONLY DATA SET
CD8.TRM.SLPS.CP = TRM.SPONT.LPS.CP %>% dplyr::select(-CD3., -CD3..CD8..PD1., -CD3..CD8..PD1..1, -CD3..CD8..CD103..CD69..PD1..2, -CD3..CD8..CD103..CD69..PD1..3)

### CD8+ DATA COLUMN RENAMING 
CD8.TRM.SLPS.CP = CD8.TRM.SLPS.CP %>% dplyr::rename(., `CD3+CD8+` = CD3..CD8..PD1..2, `CD3+CD8+PD1+` = CD3..CD8..PD1..3, `CD3+CD8+CD103+` = CD3..CD8..CD103..CD69..PD1., 
                                                    `CD3+CD8+CD103+PD1+` = CD3..CD8..CD103..CD69..PD1..1, `CD3+CD8+CD69+` = CD3..CD8..CD103..CD69..PD1..4, `CD3+CD8+CD69+PD1+` = CD3..CD8..CD103..CD69..PD1..5) 

### TAKING OUT JUST CD3+CD8+ 
CD8.TRM.SLPS.CP = CD8.TRM.SLPS.CP %>% dplyr::select(-`CD3+CD8+`)

### PIVOT CD8+ DATA LONGER
CD8.CP.PHENO = CD8.TRM.SLPS.CP %>% pivot_longer(., !HISTOLOGY & !SAMPLE_ID, names_to = "SUBTYPE", values_to = "PERCENT")

### CD8+ STACKED BAR PLOT
CD8.CP.PHENO %>% ggplot(., aes(x = HISTOLOGY, y = PERCENT, fill = SUBTYPE)) + geom_bar(position = "fill", stat = "identity") + 
  labs(x = "", y = "% of CD8+ Cells", title = "CD8+ Subtype Distribution Across Histologies", fill = "CD8+ Subtype") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), axis.title.y = element_text(face="bold", size = 14)) + scale_fill_viridis(option = "G", discrete = T)

### CD8+ PHENOTYPE DISTRIBUTION STATISTICS
library(BayesFactor)
CD8.PHENO.NS = CD8.CP.PHENO %>% select(-SAMPLE_ID)
CD8.PHENO.MATRIX = as.matrix(CD8.PHENO.NS)
chisq.test(table(CD8.PHENO.MATRIX))