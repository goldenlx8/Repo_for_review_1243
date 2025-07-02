# This script was used to analyze and generate the figures for the spectral flow cytometry - lymphocyte and myeloid populations


# Data Processing and setting up my environment
#Load libraries:
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)

# Load data 
# Lymphocytes 
data_lymph <- read.csv(file = "LPSmousemodel_Flow_Transposed_Lymphocyte_pops.csv - Sheet1.csv", header = T)
# Myeloid
data_myeloid <- read.csv(file = "LPSmousemodel_Flow_Transposed_Myeloid_pops.csv - Sheet1.csv", header = T)

# Clean lymph data (throw out samples with < 10% live cells)
data_lymph$percent_live <- (data_lymph$Live / data_lymph$Lymphocytes) * 100
filtered_data_lymph <- data_lymph  %>% filter(percent_live >= 10)
bad_miceIDs_lymph <- anti_join(data_lymph, filtered_data_lymph)
# MouseIDs lost after filtering for lymphocyte data = 1138, 1179, 1184, 1185, 1186

# Clean myeloid data (throw out samples with < 10% live cells)
data_myeloid$percent_live <- (data_myeloid$Live / data_myeloid$AllCells) * 100
filtered_data_myeloid <- data_myeloid  %>% filter(percent_live >= 10)
bad_miceIDs_myeloid <- anti_join(data_myeloid, filtered_data_myeloid)
# MouseIDs lost after filtering =  1179, 1184, 1186, 1188

# Calculate percentages for populations of interest
filtered_data_lymph$CD3_percent_of_CD45 <- (filtered_data_lymph$CD3 / filtered_data_lymph$CD45) * 100
filtered_data_lymph$CD19_percent_of_CD45 <- (filtered_data_lymph$CD19 / filtered_data_lymph$CD45) * 100
filtered_data_lymph$CD8_percent_of_CD3 <- (filtered_data_lymph$CD8 / filtered_data_lymph$CD3) * 100
filtered_data_lymph$CD4_percent_of_CD3 <- (filtered_data_lymph$CD4 / filtered_data_lymph$CD3) * 100
filtered_data_lymph$PD1_percent_of_CD8 <- (filtered_data_lymph$PD1 / filtered_data_lymph$CD8) * 100
filtered_data_lymph$LAG3_percent_of_CD8 <- (filtered_data_lymph$LAG3 / filtered_data_lymph$CD8) * 100
filtered_data_lymph$CD69_percent_of_CD8 <- (filtered_data_lymph$CD69 / filtered_data_lymph$CD8) * 100
filtered_data_lymph$CD103_percent_of_CD8 <- (filtered_data_lymph$CD103 / filtered_data_lymph$CD8) * 100
filtered_data_lymph$CD69CD103_percent_of_CD8 <- (filtered_data_lymph$CD69CD103 / filtered_data_lymph$CD8) * 100

filtered_data_myeloid$CD11b_percent_of_CD45 <- (filtered_data_myeloid$CD11b / filtered_data_myeloid$CD45) * 100
filtered_data_myeloid$CD11bF480_percent_of_CD45 <- (filtered_data_myeloid$CD11bF480 / filtered_data_myeloid$CD45) * 100
filtered_data_myeloid$F480high_percent_of_Live <- (filtered_data_myeloid$F480high / filtered_data_myeloid$Live) * 100

# Set up statistics
TEST_controls <- list(c("N1011", "N1018"), c("N1011", "N1343"), c("N1018", "N1343"))

# Set up custom theme for figures
custom_theme <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 15), axis.title.y = element_text(face = 'bold', size = 10), axis.text.x = element_text(face = "bold", size = 15))

### Lymphocytes ###
# CD8+ T cells % of CD3
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD8_percent_of_CD3)) + geom_point() + labs(x = '', y = '% of CD3+ T cells', face = 'bold', title = "CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red")) 
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD8_percent_of_CD3)) + geom_point() + labs(x = '', y = '% of CD3+ T cells', face = 'bold', title = "CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red")) 

# CD4+ T cells % of CD3
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD4_percent_of_CD3)) + geom_point() + labs(x = '', y = '% of CD3+ T cells', face = 'bold', title = "CD4 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD4_percent_of_CD3)) + geom_point() + labs(x = '', y = '% of CD3+ T cells', face = 'bold', title = "CD4 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))

# CD3+ T cells % of CD45
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD3_percent_of_CD45)) + geom_point() + labs(x = '', y = '% of CD45+ cells', face = 'bold', title = "CD3 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD3_percent_of_CD45)) + geom_point() + labs(x = '', y = '% of CD45+ cells', face = 'bold', title = "CD3 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))

# CD19+ B cells % of CD45
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD19_percent_of_CD45)) + geom_point() + labs(x = '', y = '% of CD45+ cells', face = 'bold', title = "CD19 B cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD19_percent_of_CD45)) + geom_point() + labs(x = '', y = '% of CD45+ cells', face = 'bold', title = "CD19 B cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))

# CD69 % of CD8
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD69_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "CD69+ CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD69_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "CD69+ CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))

# CD103 % of CD8
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD103_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "CD103+ CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD103_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "CD103+ CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))

# CD69/CD103 % of CD8
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD69CD103_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "CD69+/CD103+ CD8 TRMs", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = CD69CD103_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "CD69+/CD103+ CD8 TRMs", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))

# PD1 % of CD8
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = PD1_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "PD1+ CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = PD1_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "PD1+ CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))

# LAG3 % of CD8
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = LAG3_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "LAG3+ CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))
ggplot(data = filtered_data_lymph, aes(x = CellLine, y = LAG3_percent_of_CD8)) + geom_point() + labs(x = '', y = '% of CD8+ T cells', face = 'bold', title = "LAG3+ CD8 T cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red"))

### Myeloid ###
# CD11b % of CD45
ggplot(data = filtered_data_myeloid, aes(x = CellLine, y = CD11b_percent_of_CD45)) + geom_point() + labs(x = '', y = '% of CD45+ cells', face = 'bold', title = "CD11b+ Myeloid Cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red")) 
ggplot(data = filtered_data_myeloid, aes(x = CellLine, y = CD11b_percent_of_CD45)) + geom_point() + labs(x = '', y = '% of CD45+ cells', face = 'bold', title = "CD11b+ Myeloid Cells", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red")) 

# CD11b+/F480+ out of CD45
ggplot(data = filtered_data_myeloid, aes(x = CellLine, y = CD11bF480_percent_of_CD45)) + geom_point() + labs(x = '', y = '% of CD45+ cells', face = 'bold', title = "CD11b+/F480+ Macrophages", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red")) 
ggplot(data = filtered_data_myeloid, aes(x = CellLine, y = CD11bF480_percent_of_CD45)) + geom_point() + labs(x = '', y = '% of CD45+ cells', face = 'bold', title = "CD11b+/F480+ Macrophages", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red")) 

#F480high out of live
ggplot(data = filtered_data_myeloid, aes(x = CellLine, y = F480high_percent_of_Live)) + geom_point() + labs(x = '', y = '% of live cells', face = 'bold', title = "F480high TAMs", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls) + custom_theme + scale_fill_manual(values = c("blue", "black", "red")) 
ggplot(data = filtered_data_myeloid, aes(x = CellLine, y = F480high_percent_of_Live)) + geom_point() + labs(x = '', y = '% of live cells', face = 'bold', title = "F480high TAMs", size = 15)  + geom_boxplot(aes(fill = CellLine), alpha = 0.5) + stat_compare_means(comparisons = TEST_controls, label = 'p.signif') + custom_theme + scale_fill_manual(values = c("blue", "black", "red")) 

