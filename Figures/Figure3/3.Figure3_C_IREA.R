#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Cytokine activity analysis and differential enrichment score
#              comparison between MTC_123 and MTC_456 clusters (IREA output).
#              Produces Figure 3C (cytokine dot plot).
# Usage: Rscript 3.Figure3_C_IREA.R <irea_123_csv> <irea_456_csv> <output_pdf>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript 3.Figure3_C_IREA.R <irea_123_csv> <irea_456_csv> [output_pdf]")
}

irea_123_csv <- args[1]
irea_456_csv <- args[2]
output_pdf   <- if (length(args) >= 3) args[3] else "Figure3_C.pdf"

suppressPackageStartupMessages({
    library(ggplot2)
    library(cowplot)
    library(dplyr)
})

# Cytokine group color palette
TICK_Colors <- list(
  'IFNa1'='#20aa8b','IFNb'='#20aa8b','IFNe'='#20aa8b','IFNk'='#20aa8b',
  'IFNg'='#20aa8b','IFNl'='#20aa8b','IFNl1'='#20aa8b','IFNl2'='#20aa8b',
  'IL1a'='#922a6e','IL1b'='#922a6e','IL1ra'='#922a6e','IL18'='#922a6e',
  'IL33'='#922a6e','IL36a'='#922a6e','IL36RA'='#922a6e',
  'IL2'='#F6C76D','IL4'='#F6C76D','IL13'='#F6C76D','IL15'='#F6C76D',
  'IL7'='#F6C76D','TSLP'='#F6C76D','IL9'='#F6C76D','IL21'='#F6C76D',
  'IL3'='#8668c0','IL5'='#8668c0','GM-CSF'='#8668c0',
  'IL6'='#198591','IL11'='#198591','IL27'='#198591','IL30'='#198591',
  'IL31'='#198591','LIF'='#198591','OSM'='#198591','Cardiotrophin-1'='#198591',
  'Neuropoietin'='#198591','IL12'='#198591','IL23'='#198591','IL-Y'='#198591',
  'IL10'='#df70b0','IL19'='#df70b0','IL20'='#df70b0','IL22'='#df70b0','IL24'='#df70b0',
  'IL17A'='#096954','IL17B'='#096954','IL17C'='#096954','IL17D'='#096954',
  'IL17E'='#096954','IL17F'='#096954',
  'Flt3l'='#472d8c','IL34'='#472d8c','M-CSF'='#472d8c','G-CSF'='#472d8c',
  'SCF'='#472d8c','EGF'='#472d8c','VEGF'='#472d8c','FGF-basic'='#472d8c',
  'HGF'='#472d8c','IGF-I'='#472d8c',
  'LTA1-B2'='#53b8c4','LTA2-B1'='#53b8c4','TNFa'='#53b8c4','OX40L'='#53b8c4',
  'CD40L'='#53b8c4','FasL'='#53b8c4','CD27L'='#53b8c4','CD30L'='#53b8c4',
  '41BBL'='#53b8c4','TRAIL'='#53b8c4','RANKL'='#53b8c4','TWEAK'='#53b8c4',
  'APRIL'='#53b8c4','BAFF'='#53b8c4','LIGHT'='#53b8c4','TL1A'='#53b8c4','GITRL'='#53b8c4',
  'C3a'='#b7572e','C5a'='#b7572e',
  'Prolactin'='#91a160','Leptin'='#91a160','Adiponectin'='#91a160',
  'Resistin'='#91a160','TGF-beta-1'='#91a160','GDNF'='#91a160',
  'Persephin'='#91a160','Noggin'='#91a160','Decorin'='#91a160','TPO'='#91a160'
)

Control <- read.csv(irea_123_csv)
CASE    <- read.csv(irea_456_csv)
Control$Cond <- "MTC_123"
CASE$Cond    <- "MTC_456"
df <- rbind(Control, CASE)

df$Log10Qval <- -log10(df$padj + 1)
df$directed  <- df$Log10Qval * ifelse(df$Cond == 'MTC_123', -1, 1)
df$Significant <- ifelse(df$padj < 0.05, '+', "")

my_dataframe <- data.frame(
    Cytokine = names(TICK_Colors),
    Groups   = as.vector(unlist(TICK_Colors))
)
my_dataframe_sub <- my_dataframe %>% filter(Cytokine %in% df$Cytokine)
rownames(my_dataframe_sub) <- my_dataframe_sub$Cytokine

order_ticks  <- unique(my_dataframe_sub$Cytokine)
df$Cytokine  <- factor(df$Cytokine, levels = order_ticks)
df_ord       <- df[order(df$Cytokine), ]
df_ord$Groups <- my_dataframe_sub[as.character(df_ord$Cytokine), "Groups"]
df_ord$ES    <- df_ord$Enrichment.Score

Fig3C <- ggplot(df_ord, aes(x = Cytokine, y = Cond,
                             color = ES, fill = ES,
                             size = abs(ES), label = Significant)) +
    geom_point() +
    geom_text(alpha = 1, size = 4, colour = "white") +
    RotatedAxis() +
    theme_cowplot() +
    theme(
        text         = element_text(size = 12),
        axis.text.x  = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y  = element_text(size = 12),
        plot.title   = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text.x.bottom = element_text(colour = my_dataframe_sub$Groups)
    ) +
    scale_color_gradient2(low = 'navy', mid = 'white', high = 'firebrick3') +
    xlab("") + ylab("")

pdf(output_pdf, width = 16, height = 3.5)
print(Fig3C)
dev.off()

message("Saved: ", output_pdf)
