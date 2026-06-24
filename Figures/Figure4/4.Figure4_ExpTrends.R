#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Mixed-effects models for in-vivo behavioural outcomes
#              (Rotarod RPM/latency, MNSS) across timepoints and treatment conditions.
#              Produces Figure 4 line plots.
# Usage: Rscript 4.Figure4_ExpTrends.R <input_csv> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 4.Figure4_ExpTrends.R <input_csv> [output_dir]")
}

input_csv  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "."

suppressPackageStartupMessages({
    library(tidyverse)
    library(lme4)
    library(lmerTest)
    library(ggplot2)
    library(emmeans)
    library(ggpubr)
})

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(input_csv)
# Coerce numeric columns
num_cols <- c(
    "OFVelocity_Day3","OFVelocity_Day7",
    "OFDistance_Day3","OFDistance_Day7",
    "OFInZoneOZ_Day3","OFInZoneOZ_Day7",
    "OFInZoneCorners_Day3","OFInZoneCorners_Day7",
    "OFInZoneBoarders_Day3","OFInZoneBoarders_Day7",
    "RotarodLatencyMeans_Day3","RotarodLatencyMeans_Day7",
    "RotarodLatencyMedians_Day3","RotarodLatencyMedians_Day7",
    "RotarodRPMMeans_Day3","RotarodRPMMeans_Day7",
    "RotarodRPMMedians_Day3","RotarodRPMMedians_Day7"
)
for (col in num_cols) if (col %in% colnames(data)) data[[col]] <- as.numeric(data[[col]])
data[, grepl("^X(\\.\\d+)?$", colnames(data))] <- NULL

condition_colors <- c("Control" = "black", "Everolimus" = "#CC79A7", "RA" = "#E69F00")
base_theme <- theme_bw() + theme(
    text         = element_text(family = "Helvetica", size = 6),
    axis.title   = element_text(size = 6),
    axis.text    = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 6),
    strip.text   = element_text(size = 6)
)

# ===================== RotarodRPMMedians =====================
df <- data %>%
    filter(Condition %in% c('Control','RA','Everolimus')) %>%
    select(Animal.Code, Condition, Gender, Dose,
           RotarodRPMMedians_BSL, RotarodRPMMedians_Day3, RotarodRPMMedians_Day7) %>%
    pivot_longer(cols = starts_with("RotarodRPMMedians"),
                 names_to = "Timepoint", values_to = "Score") %>%
    mutate(
        Timepoint = recode(Timepoint,
            RotarodRPMMedians_BSL  = "BSL",
            RotarodRPMMedians_Day3 = "Day3",
            RotarodRPMMedians_Day7 = "Day7"),
        Timepoint = factor(Timepoint, levels = c("BSL","Day3","Day7"), ordered = TRUE),
        Timepoint_Numeric = as.numeric(Timepoint),
        Gender    = as.factor(Gender),
        Condition = as.factor(Condition),
        Dose_scaled = scale(Dose)
    )

model <- glmer(Score ~ Timepoint_Numeric * Condition + Dose_scaled + (1 | Animal.Code),
               data = df %>% filter(Condition %in% c('Control','RA','Everolimus')),
               family = Gamma(link = 'log'))
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
print(summary(emtrends_model))

p_rpm <- ggplot(df, aes(x = Timepoint, y = Score, color = Condition, group = Condition)) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    labs(title = "Rotarod RPM Medians", x = "Timepoint", y = "RPM") +
    scale_color_manual(values = condition_colors) +
    base_theme

ggsave(file.path(output_dir, "RotarodRPMMedians_All_Lines.pdf"),
       plot = p_rpm, width = 6, height = 4, device = cairo_pdf)

# ===================== RotarodLatencyMedians =====================
df2 <- data %>%
    filter(Condition %in% c('Control','RA','Everolimus')) %>%
    select(Animal.Code, Condition, Gender, Dose, contains('RotarodLatencyMedians')) %>%
    pivot_longer(cols = starts_with("RotarodLatencyMedians"),
                 names_to = "Timepoint", values_to = "Score") %>%
    mutate(
        Timepoint = recode(Timepoint,
            RotarodLatencyMedians_BSL  = "BSL",
            RotarodLatencyMedians_Day3 = "Day3",
            RotarodLatencyMedians_Day7 = "Day7"),
        Timepoint = factor(Timepoint, levels = c("BSL","Day3","Day7"), ordered = TRUE),
        Timepoint_Numeric = as.numeric(Timepoint),
        Gender    = as.factor(Gender),
        Condition = as.factor(Condition),
        Dose_scaled = scale(Dose)
    )

model2 <- lmer(Score ~ Timepoint_Numeric * Condition + Dose_scaled + (1 | Animal.Code),
               data = df2 %>% filter(Condition %in% c('Control','RA','Everolimus')))
emtrends_model2 <- emtrends(model2, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
print(summary(emtrends_model2))

p_lat <- ggplot(df2, aes(x = Timepoint, y = Score, color = Condition, group = Condition)) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    labs(title = "Rotarod Latency Medians", x = "Timepoint", y = "Latency (s)") +
    scale_color_manual(values = condition_colors) +
    base_theme

ggsave(file.path(output_dir, "RotarodLatencyMedians_All_Lines.pdf"),
       plot = p_lat, width = 6, height = 4, device = cairo_pdf)

# ===================== MNSS =====================
df3 <- data %>%
    filter(Condition %in% c('Control','RA','Everolimus')) %>%
    select(Animal.Code, Condition, Gender, Dose, contains('MNSS')) %>%
    pivot_longer(cols = starts_with("MNSS"),
                 names_to = "Timepoint", values_to = "Score") %>%
    mutate(
        Timepoint = recode(Timepoint,
            MNSS_BSL  = "BSL", MNSS_2h = "2h",
            MNSS_Day1 = "Day1", MNSS_Day3 = "Day3", MNSS_Day7 = "Day7"),
        Timepoint = factor(Timepoint, levels = c("BSL","2h","Day1","Day3","Day7"), ordered = TRUE),
        Timepoint_Numeric = as.numeric(Timepoint),
        Gender    = as.factor(Gender),
        Condition = as.factor(Condition),
        Dose_scaled = scale(Dose)
    )

model3 <- lmer(Score ~ Timepoint_Numeric * Condition + Dose_scaled + (1 | Animal.Code),
               data = df3 %>% filter(Condition %in% c('Control','RA','Everolimus')))
emtrends_model3 <- emtrends(model3, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
print(summary(emtrends_model3))

p_mnss <- ggplot(df3, aes(x = Timepoint, y = Score, color = Condition, group = Condition)) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    labs(title = "MNSS Score", x = "Timepoint", y = "MNSS") +
    scale_color_manual(values = condition_colors) +
    base_theme

ggsave(file.path(output_dir, "MNSS_Score_All_Lines.pdf"),
       plot = p_mnss, width = 6, height = 4, device = cairo_pdf)

message("Saved outputs to: ", output_dir)
