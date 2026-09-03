# ============================================================================
# Figure3_A_irea_cytokine_plot.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 3, Panel A — IREA cytokine-response enrichment for MTC_456 vs
# MTC_123. Enrichment score (ES) for each of the 86 Immune Dictionary
# cytokine responses; dot size/color = ES, white "+" = significant, family
# strip on top.
#
# Input:  TABLE_DIR/Figure3_IREA_Results.csv  (manually downloaded from IREA
#         after running Figure3_A_irea_export.R — see that script)
# Output: FIG_DIR/Figure3_A_irea_cytokine_plot.png
# ==============================================================================

source("../../setup/00_setup.R")
source("../../setup/00_theme_colors.R")

CASE <- read.csv(file.path(TABLE_DIR, "Figure3_IREA_Results.csv"))
CASE$Cytokine <- enc2utf8(as.character(CASE$Cytokine))

CASE$Cond <- "MTC456_vs_MTC123"
CASE$ES <- CASE$Enrichment.Score
CASE$Significant <- ifelse(CASE$padj < 0.05, "+", "")

my_dataframe <- data.frame(
  Cytokine = names(Cytokine_Family),
  Family   = as.vector(Cytokine_Family),
  stringsAsFactors = FALSE
)
my_dataframe$Groups <- Family_Colors[my_dataframe$Family]

my_dataframe_sub <- my_dataframe %>% filter(Cytokine %in% CASE$Cytokine)
rownames(my_dataframe_sub) <- my_dataframe_sub$Cytokine
order_ticks <- unique(my_dataframe_sub$Cytokine)

n_missing <- sum(!(CASE$Cytokine %in% names(Cytokine_Family)))
if (n_missing > 0) {
  warning(paste(n_missing, "cytokines in CASE do not match Cytokine_Family and will be dropped:",
                paste(setdiff(CASE$Cytokine, names(Cytokine_Family)), collapse = ", ")))
}

CASE$Cytokine <- factor(CASE$Cytokine, levels = order_ticks)
df_ord <- CASE[order(CASE$Cytokine), ]
df_ord$Groups <- my_dataframe_sub[as.character(df_ord$Cytokine), "Groups"]
df_ord$Family <- my_dataframe_sub[as.character(df_ord$Cytokine), "Family"]
df_ord$Family <- factor(df_ord$Family, levels = names(Family_Colors))
df_ord$ES <- df_ord$Enrichment.Score

strip <- ggplot(df_ord, aes(x = Cytokine, y = 1, fill = Family)) +
  geom_tile() +
  scale_fill_manual(values = Family_Colors, name = "Cytokine Families", drop = FALSE) +
  theme_void()
strip_notitle <- strip + theme(legend.position = "none")

p_A_irea <- ggplot(df_ord, aes(x = Cytokine, y = Cond, color = ES, size = abs(ES), label = Significant)) +
  geom_point(aes(fill = Family), shape = 22, size = 0, alpha = 0, color = NA) +
  geom_point() +
  geom_text(alpha = 1, size = 4, colour = "white") +
  scale_fill_manual(values = Family_Colors, name = "Cytokine Families", drop = FALSE) +
  RotatedAxis() + theme_cowplot() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1,
                                    colour = my_dataframe_sub[order_ticks, "Groups"]),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.x = unit(0.7, "cm"),
        legend.spacing.y = unit(0.2, "cm")) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(alpha = 1, size = 5, color = NA), order = 1),
    color = guide_colorbar(order = 2),
    size = guide_legend(order = 3)
  ) +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3") +
  xlab("") + ylab("")

p_A_final <- plot_grid(strip_notitle, p_A_irea, ncol = 1, align = "v", axis = "lr",
                        rel_heights = c(0.06, 1))

ggsave(file.path(FIG_DIR, "Figure3_A_irea_cytokine_plot.png"),
       p_A_final, width = 16, height = 5, dpi = 300, type = "cairo")
ggsave(file.path(FIG_DIR, "Figure3_A_irea_cytokine_plot.pdf"),
       p_A_final, width = 16, height = 5, dpi = 300, units = "in",device =grDevices::cairo_pdf)
cat("Saved: Figure3_A_irea_cytokine_plot.png\n")
