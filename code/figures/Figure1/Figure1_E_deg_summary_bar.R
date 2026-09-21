# ============================================================================
# Figure1_E_deg_summary_bar.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Figure 1, Panel E — significant DEG counts per cell type (favorable vs
# unfavorable outcome direction)
#
# Input:  TABLE_DIR/Figure1_master_DEGs_all_assays.csv  (from Figure1_master_DEGs.R)
# Output: TABLE_DIR/Figure1_E_deg_summary_counts.csv
#         FIG_DIR/Figure1_E_deg_summary_bar.(pdf|png)
# ==============================================================================

source("../../00_setup.R")
source("../../00_helpers_stats.R")

DEGs_file <- read.csv(file.path(TABLE_DIR, "Figure1_master_DEGs_all_assays.csv"))

deg_summary <- DEGs_file %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(direction = case_when(logFC > 0 ~ "Unfavorable", logFC < 0 ~ "Favorable")) %>%
  count(assay, direction) %>%
  complete(assay, direction = c("Favorable", "Unfavorable"), fill = list(n = 0)) %>%
  group_by(assay) %>%
  mutate(total_DEGs = sum(n),
         n_plot = ifelse(direction == "Favorable", -n, n),
         label = abs(n_plot)) %>%
  ungroup()

assay_order <- deg_summary %>%
  distinct(assay, total_DEGs) %>%
  arrange(total_DEGs) %>%
  pull(assay)

deg_summary <- deg_summary %>% mutate(assay = factor(assay, levels = assay_order))

p_E_bar <- ggplot(deg_summary, aes(x = assay, y = n_plot, fill = direction)) +
  geom_col(width = 0.65, color = NA) +
  geom_text(
    data = deg_summary %>% filter(label >= 2),
    aes(label = label, y = n_plot),
    hjust = ifelse(deg_summary$n_plot > 0, -0.25, 1.25),
    size = 3.5
  ) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  coord_flip() +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.15, 0.15))) +
  scale_fill_manual(values = c("Favorable" = "#4C78A8", "Unfavorable" = "#D95F5F")) +
  labs(x = NULL, y = "Significant DEGs (FDR < 5%)", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        axis.text = element_text(color = "black"),
        plot.margin = margin(10, 25, 10, 10)) +
  theme_natmed() +
  theme(legend.position = c(1.23, 0.25),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key = element_rect(fill = alpha("white", 0)))

save_table_safe(deg_summary, "Figure1_E_deg_summary_counts.csv")
save_plot_safe(p_E_bar, "Figure1_E_deg_summary_bar", width = 4, height = 4)
