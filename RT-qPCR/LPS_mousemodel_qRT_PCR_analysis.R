#This script was used to analyze and generate the figures for qRT-PCR in the paper

# Load necessary libraries
library(dplyr)
library(ggbreak)
library(ggplot2)
library(RColorBrewer)

# Read the data
qrt_data_raw <- read.csv("LPSmousemodelpaper_qRTPCR_data_raw - Sheet1.csv", header = TRUE)

# Replace "Undetermined" with 40
qrt_data <- qrt_data_raw %>%
  mutate(Cq = if_else(Cq == "Undetermined", 40, as.numeric(Cq)))

# Filter for GAPDH and calculate the average Ct per sample
gapdh_data <- qrt_data %>%
  filter(Target == "GAPDH") %>%
  group_by(Sample) %>%
  summarize(GAPDH_Ct = mean(Cq, na.rm = TRUE))

# Merge the GAPDH Ct values with the main data
qrt_data <- qrt_data %>%
  left_join(gapdh_data, by = "Sample") %>%
  mutate(Delta_Ct = Cq - GAPDH_Ct)

# Trp53/Pten analysis:
# Filter for Trp53/Pten data 
trp53_pten_data <- qrt_data %>%
  filter(Target == "Trp53" | Target == "PTEN")

# Set reference for Trp53 and PTEN analysis (3T3 D9)
ref_dCt_Trp53Pten <- trp53_pten_data %>%
  filter(Sample == "3T3 D9") %>%
  group_by(Target) %>%
  summarize(ref_dCt_Trp53Pten = mean(Delta_Ct, na.rm = TRUE))

# Calculate ΔΔCt
trp53_pten_data <- trp53_pten_data %>%
  left_join(ref_dCt_Trp53Pten, by = "Target") %>%
  mutate(Delta_Delta_Ct = Delta_Ct - ref_dCt_Trp53Pten)

# Calculate relative quantification (RQ)
trp53_pten_data <- trp53_pten_data %>%
  mutate(RQ = 2^(-Delta_Delta_Ct))

# Adipocyte marker analysis:
# Split adipocyte markers data from adipocyte marker data (AdipoQ, Cebpa, FabP4, Pparg)
adipomarkers_data <- qrt_data %>%
  filter(Target == "AdipoQ" | Target == "Cebpa" | Target == "Pparg" | Target == "FabP4")

#Set reference for adipocyte analysis (3T3 D0)
ref_dCt_adipocyte <- adipomarkers_data %>%
  filter(Sample == "3T3 D0") %>%
  group_by(Target) %>%
  summarize(ref_dCt_adipocyte = mean(Delta_Ct, na.rm = TRUE))

# Calculate ΔΔCt
adipomarkers_data <- adipomarkers_data %>%
  left_join(ref_dCt_adipocyte, by = "Target") %>%
  mutate(Delta_Delta_Ct = Delta_Ct - ref_dCt_adipocyte)

# Calculate relative quantification (RQ)
adipomarkers_data <- adipomarkers_data %>%
  mutate(RQ = 2^(-Delta_Delta_Ct))

# Specify the desired order of samples
sample_order_adipo <- c("3T3 D0", "3T3 D2", "3T3 D9", "N1011", "N1018", "N1343") 
sample_order_trp53pten <- c("3T3 D9", "N1011", "N1018", "N1343") 
adipomarkers_data$Sample <- factor(adipomarkers_data$Sample, levels = sample_order_adipo)
trp53_pten_data$Sample <- factor(trp53_pten_data$Sample, levels = sample_order_trp53pten)

#Drop NAs
trp53_pten_data <- trp53_pten_data %>%
  filter(!is.na(Sample))

# Define the colors for each sample
sample_colors <- c(
  "3T3 D0" = "grey88",
  "3T3 D2" = "grey45",
  "3T3 D9" = "grey20",
  "N1011" = "blue",
  "N1018" = "grey",
  "N1343" = "red"
)

#Build summary stats
summary_stats_trp53pten <- trp53_pten_data %>%
  group_by(Sample, Target) %>%
  reframe(
    mean_RQ = mean(RQ, na.rm = TRUE),
    sd_RQ = sd(RQ, na.rm = TRUE),
    sem_RQ = sd_RQ / sqrt(n()),  # Standard Error of Mean (SEM)
    RQ = RQ
  )

summary_stats_adipo <- adipomarkers_data %>%
  group_by(Sample, Target) %>%
  reframe(
    mean_RQ = mean(RQ, na.rm = TRUE),
    sd_RQ = sd(RQ, na.rm = TRUE),
    sem_RQ = sd_RQ / sqrt(n()),  # Standard Error of Mean (SEM)
    RQ = RQ
  )


# Build figure function
plot_qpcr_bar_trp53pten <- function(gene_name, data, y_limits, y_breaks, y_cut) {
  ggplot(data %>% filter(Target == gene_name), aes(x = Sample, y = RQ, fill = Sample)) +
    geom_bar(stat = "summary", fun = "mean", position = position_dodge(), alpha = 0.9) +  # Mean bar
    geom_jitter(color = "black", width = 0.3, size = 0.1, alpha = 0.8) +  # Individual data points
    geom_errorbar(data = summary_stats_trp53pten %>% filter(Target == gene_name), 
                  aes(x = Sample, ymin = mean_RQ - sem_RQ, ymax = mean_RQ + sem_RQ), 
                  width = 0.4, color = "black") +  # Error bars (SEM)
    theme_minimal() +
    labs(title = NULL, 
         x = NULL, 
         y = "Fold Change (RQ)") +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    #scale_y_cut(breaks = y_cut) +
    theme(
      axis.text.x = element_text(family = "Arial", size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(family = "Arial", size = 6),
      axis.title.y = element_text(family = "Arial", size = 8),
      legend.position = "none"
    ) +
    scale_fill_manual(values = sample_colors) 
}

plot_qpcr_bar_adipo <- function(gene_name, data, y_limits, y_breaks, y_cut) {
  ggplot(data %>% filter(Target == gene_name), aes(x = Sample, y = RQ, fill = Sample)) +
    geom_bar(stat = "summary", fun = "mean", position = position_dodge(), alpha = 0.9) +  # Mean bar
    geom_jitter(color = "black", width = 0.3, size = 0.1, alpha = 0.8) +  # Individual data points
    geom_errorbar(data = summary_stats_adipo %>% filter(Target == gene_name), 
                  aes(x = Sample, ymin = mean_RQ - sem_RQ, ymax = mean_RQ + sem_RQ), 
                  width = 0.4, color = "black") +  # Error bars (SEM)
    theme_minimal() +
    labs(title = NULL, 
         x = NULL, 
         y = "Fold Change (RQ)") +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    #scale_y_cut(breaks = y_cut) +
    theme(
      axis.text.x = element_text(family = "Arial", size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(family = "Arial", size = 6),
      axis.title.y = element_text(family = "Arial", size = 8),
      legend.position = "none"
    ) +
    scale_fill_manual(values = sample_colors) 
}
#Generate plots
Pparg_plot <- plot_qpcr_bar_adipo("Pparg", adipomarkers_data, y_limits = c(0, 25.5), y_breaks = c(seq(0, 1.5, by = 0.5), seq(1.5, 25.5, by = 4)), y_cut = 1.5)
Cebpa_plot <- plot_qpcr_bar_adipo("Cebpa", adipomarkers_data, y_limits = c(0, 35), y_breaks = seq(0, 35,by = 5))
AdipoQ_plot <- plot_qpcr_bar_adipo("AdipoQ", adipomarkers_data, y_limits = c(0, 11000), y_breaks = c(seq(0, 10, by = 1), seq(10, 200, by = 50), seq(200, 11000, by = 2000)), y_cut = c(10,200))
FabP4_plot <- plot_qpcr_bar_adipo("FabP4", adipomarkers_data, y_limits = c(0, 200), y_breaks = c(seq(0, 2, by = 0.2), seq(10,200,by = 50)), y_cut = 2)
P53_plot <- plot_qpcr_bar_trp53pten("Trp53", trp53_pten_data, y_limits = c(0, 1.25), y_breaks = seq(0, 1.25, by = 0.2))
PTEN_plot <- plot_qpcr_bar_trp53pten("PTEN", trp53_pten_data, y_limits = c(0, 1.25), y_breaks = seq(0, 1.25, by = 0.2))

#View figures
Pparg_plot
Cebpa_plot
AdipoQ_plot
FabP4_plot
P53_plot
PTEN_plot

#Save figures
#ggsave("Pparg_qPCR_barplot.png", plot = Pparg_plot, width = 2, height = 2.75, dpi = 600)
#ggsave("Cebpa_qPCR_barplot.png", plot = Cebpa_plot, width = 2, height = 2.75, dpi = 600)
#ggsave("AdipoQ_qPCR_barplot.png", plot = AdipoQ_plot, width = 2, height = 2.75, dpi = 600)
#ggsave("FabP4_qPCR_barplot.png", plot = FabP4_plot, width = 2, height = 2.75, dpi = 600)
#ggsave("P53_qPCR_barplot.png", plot = P53_plot, width = 2.5, height = 4, dpi = 600)
#ggsave("PTEN_qPCR_barplot.png", plot = PTEN_plot, width = 2.5, height = 4, dpi = 600)








