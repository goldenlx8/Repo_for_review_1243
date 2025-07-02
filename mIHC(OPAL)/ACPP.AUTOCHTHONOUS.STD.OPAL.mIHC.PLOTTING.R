# title: "MS LPS OPAL SPONTANEOUS CELL % ANALYSIS"
# author: "EMMA KENNA (ekenna@umich.edu)"

### LIBRARIES
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(huxtable)

### LOADING DATA
SPONT_LPS_CP = read.csv(file = '2024.05.13.SUM.CELL.PERCENTS.csv', sep = ",", header = TRUE)
TEST = list(c("DDLPSlow", "WDLPS"), c("DDLPShigh", "DDLPSlow"), c("DDLPShigh", "WDLPS"))
SPONT_LPS_CP$HISTOLOGY = factor(SPONT_LPS_CP$HISTOLOGY, levels = c("WDLPS", "DDLPSlow", "DDLPShigh"))

### CD3+ PLOTTING & STATISTICS
SPONT_LPS_CP %>% select(CD3., HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = CD3.)) + stat_summary(fun = mean, geom = "bar", alpha = 0.8, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, label = "p.signif") + 
  labs(x = "", y = "% CD3+ / Total Cells", face = "bold", title = "CD3+ T Cells") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), axis.title.y = element_text(face="bold", size = 15)) 

### CD3+/CD8a+ PLOTTING & STATISTICS
SPONT_LPS_CP %>% select(CD3..CD8a., HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = CD3..CD8a.)) + stat_summary(fun = mean, geom = "bar", alpha = 0.8, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, label = "p.signif") + labs(x = "", y = "% CD3+CD8a+ / Total Cells", title = "CD8a+ T Cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), 
        axis.title.y = element_text(face="bold", size = 15)) 

### CD3+/FOXP3+ PLOTTING & STATISTICS
SPONT_LPS_CP %>% select(CD3..FOXP3., HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = CD3..FOXP3.)) + stat_summary(fun = mean, geom = "bar", alpha = 0.8, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, label = "p.signif") + 
  labs(x = "", y = "% CD3+FOXP3+ / Total Cells", title = "Regulatory T Cells") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), 
        axis.title.y = element_text(face="bold", size = 15)) 

### B CELL [CD19+] PLOTTING & STATISTICS 
SPONT_LPS_CP %>% select(CD19., HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = CD19.)) + stat_summary(fun = mean, geom = "bar", alpha = 0.8, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, label = "p.signif") + 
  labs(x = "", y = "% CD19+ / Total Cells", title = "B Cells") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), 
        axis.title.y = element_text(face="bold", size = 15)) 

### MACROPHAGE [F480+] PLOTTING & STATISTICS
SPONT_LPS_CP %>% select(F480., HISTOLOGY) %>%
  ggplot(., aes(x = HISTOLOGY, y = F480.)) + stat_summary(fun = mean, geom = "bar", alpha = 0.8, fill = c("red", "dodgerblue", "blue3")) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.3) + 
  geom_jitter() + stat_compare_means(comparisons = TEST, label = "p.signif") + 
  labs(x = "", y = "% F480+ / Total Cells", title = "Macrophages") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), axis.text.x = element_text(face="bold", size = 15), 
        axis.title.y = element_text(face="bold", size = 15)) 

# TO SEE THE NUMERICAL P-VALUES REMOVE THE label = "p.signif" argument from the stat_compare_means_function 