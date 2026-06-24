#!/usr/bin/env Rscript
# Author: Dimitrios Kyriakis
# Description: Supplementary Figure 1 — donor-level QC panels and cohort
#              demographic pie charts with anonymised donor IDs.
# Usage: Rscript SFigure1.R <seurat_rds> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript SFigure1.R <seurat_rds> [output_dir]")
}

input_rds  <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "."

suppressPackageStartupMessages({
    library(ggplot2)
    library(Seurat)
    library(cowplot)
    library(patchwork)
    library(dplyr)
})

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

mydata <- readRDS(input_rds)

# Anonymise donor IDs
Real_donor_ids <- sort(unique(mydata@meta.data$donor))
Crypto_ids <- sprintf("Donor_%02d", seq_along(Real_donor_ids))
names(Crypto_ids) <- Real_donor_ids
df <- mydata@meta.data
df$donor <- unlist(lapply(df$donor, function(x) Crypto_ids[x]))
mydata$donor <- df$donor

# ===================== Panel A — QC violin + metadata tiles =====================
g1 <- VlnPlot(mydata, features = c("nGenes","nCounts"), stack = TRUE, group.by = "donor") & theme_bw()
g2 <- VlnPlot(mydata, features = c("log.genes","log.umis"), stack = TRUE, group.by = "donor") & theme_bw()
g3 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
    geom_tile(aes(fill = Kit), colour = "black") + theme_bw() + coord_flip() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylab("Kit") + scale_fill_grey(start = 0.8, end = 0.2, na.value = NA) +
    theme(axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = unit(c(0,0,0,0), "cm"))
g4 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
    geom_tile(aes(fill = sex), colour = "black") + theme_bw() + coord_flip() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("sex") + scale_fill_manual(values = c("#E69F00","#56B4E9"), na.value = NA) +
    theme(axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = unit(c(0,0,0,0), "cm"))
g5 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
    geom_tile(aes(fill = race), colour = "black") + theme_bw() + coord_flip() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("race") + scale_fill_manual(values = c("#FFFE99","#F5C085","#BDADD4","#7FC97F"), na.value = NA) +
    theme(axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = unit(c(0,0,0,0), "cm"))
g6 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
    geom_tile(aes(fill = as.factor(mRS)), colour = "black") + theme_bw() + coord_flip() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("mRS") +
    scale_fill_manual(values = c("#4D4D4D","#999999","#E0E0E0","#F9DAC7","#EF8962","#B3242B"), na.value = NA) +
    theme(axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = unit(c(0,0,0,0), "cm"))
g7 <- ggplot(mydata@meta.data, aes(x = donor)) +
    geom_bar() + theme_bw() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    coord_flip()

sfig1A <- g1 + g3 + g4 + g5 + g6 + g7 +
    plot_layout(guides = "collect", widths = c(2, 0.2, 0.2, 0.2, 0.2, 1.5)) &
    theme(legend.position = 'right')

# ===================== Panel B — demographic pie charts =====================
Categories <- c('sex','race','BL','CL','dx','mRS')
df_unique <- df[, c('orig.ident','donor', Categories)] %>% distinct(donor, .keep_all = TRUE)
df_unique_ord <- df_unique[order(df_unique$donor), ]
rownames(df_unique_ord) <- df_unique_ord$donor

create_pie_chart <- function(category) {
    category_data <- df_unique_ord[[category]]
    summary_data  <- as.data.frame(table(category_data))
    colnames(summary_data) <- c(category, "Count")
    summary_data$Percentage <- round(100 * summary_data$Count / sum(summary_data$Count))
    p <- ggplot(summary_data, aes(x = "", y = Count, fill = .data[[category]])) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar(theta = "y") +
        geom_text(aes(label = paste0(Percentage, "%")),
                  position = position_stack(vjust = 0.5), size = 3) +
        labs(title = paste("Pie Chart:", category), fill = category) +
        theme_void() +
        theme(legend.position = "right", legend.text = element_text(size = 6))
    if (category == 'sex')  p <- p + scale_fill_manual(values = c("#E69F00","#56B4E9"), na.value = NA)
    if (category == 'race') p <- p + scale_fill_manual(values = c("#FFFE99","#F5C085","#BDADD4","#7FC97F"), na.value = NA)
    if (category == 'dx')   p <- p + scale_fill_manual(values = c("#E78AC3","#8DA1CB","#EF8C61"), na.value = NA)
    if (category == 'mRS')  p <- p + scale_fill_manual(values = c("#4D4D4D","#999999","#E0E0E0","#F9DAC7","#EF8962","#B3242B"), na.value = NA)
    p
}

plots  <- lapply(Categories, create_pie_chart)
sfig1B <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) / (plots[[5]] + plots[[6]])

pdf(file.path(output_dir, "SupFigure1_A.pdf"), width = 8, height = 6)
print(sfig1A + theme(legend.position = "right", legend.text = element_text(size = 5)))
dev.off()

pdf(file.path(output_dir, "SupFigure1_B.pdf"), width = 8, height = 6)
print(sfig1B)
dev.off()

pdf(file.path(output_dir, "SupFigure1.pdf"), width = 12, height = 6)
print(sfig1A + sfig1B + plot_layout(widths = c(2, 0.2, 0.2, 0.2, 0.2, 1.5, 5)))
dev.off()

message("Saved outputs to: ", output_dir)
