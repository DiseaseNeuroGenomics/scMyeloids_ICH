# ============================================================================
# SFigure1_donor_annotations.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Supplementary Figure 1 — donor-level annotation strip (sex, ancestry, age,
# clot location, biopsy location, mRS, diagnosis) plus per-donor QC
# distributions (avg genes/counts) and cell count; companion demographic pie
# charts.
#
# Adapted from a raw uploaded script (SFigure1.R) for this project's
# conventions. Fixes applied:
#   - Input path was a local Mac path (readRDS('/Users/.../Data_backup/...'))
#     — now reads from BASE_DIR, same file every other Figure1_*.R script uses.
#   - library(dplyr) was loaded at the very end (line ~237 of the original)
#     but %>%/group_by() were already used earlier — only worked interactively
#     because dplyr was already attached from prior session state. 00_setup.R
#     loads dplyr globally, so this is fixed for free here.
#   - Dropped dead/superseded code: an early VlnPlot exploration (g1/g1alt)
#     and a standalone CL tile (pl_CL) that were both superseded by later
#     re-definitions before ever being used; an unused GrAge computation;
#     a broken dangling `colors_list <-` line in create_pie_chart() (unused
#     even before being broken — every category already has a hardcoded
#     scale_fill_manual() below it).
#   - Consolidated ~6 scattered, partly-broken PDF-only save attempts
#     (some absolute Mac paths, some fragile "../../../SFigures/..." relative
#     paths) into this project's save_plot_safe() (pdf + png, into FIG_DIR).
#     Dropped one further "combined A+B" save attempt whose plot_layout()
#     width vector (7 values) doesn't match its actual panel count once two
#     multi-panel patchworks are added together (~12 leaf panels) — looked
#     like an unfinished experiment, not a working output.
#
# KNOWN ISSUE, NOT FIXED — CONFIRM INTENDED: the original script builds a
# `sex` tile plot (g2) but then overwrites the g2 variable with the
# avg-nCounts plot before the final composite is built, so the actual final
# figure never includes sex, despite the manuscript's own S.Figure 1
# description saying it shows "biological sex." Preserved as-is (matching
# current actual behavior) rather than guessing where to reinsert it.
#
# Input:  BASE_DIR/2024_11_Final_Immune.rds
# Output: FIG_DIR/SFigure1_A_donor_annotations.(pdf|png)
#         FIG_DIR/SFigure1_B_demographic_pies.(pdf|png)
#         TABLE_DIR/SFigure1_STable1_donor_metadata.csv  (S.Table 1)
# ==============================================================================

source("../../00_setup.R")

mydata <- readRDS(file.path(BASE_DIR, "2024_11_Final_Immune.rds"))

# ---- Data repair for one donor with missing metadata ----
mydata$mRS[mydata$donor == "CK-038"] <- "3"
mydata$mRS_binned[mydata$donor == "CK-038"] <- "G"
mydata$mRS_3class[mydata$donor == "CK-038"] <- "Mid"
mydata$BL[mydata$donor == "CK-038"] <- "Parietal"
mydata$CL[mydata$donor == "CK-038"] <- "Parietal"
mydata$TSH[mydata$donor == "CK-038"] <- 35.7
mydata$dx[mydata$donor == "CK-038"] <- "AD"
mydata$age[mydata$donor == "CK-038"] <- 74
mydata$sex[mydata$donor == "CK-038"] <- "Male"
mydata$race[mydata$donor == "CK-038"] <- "Asian"
mydata$PHV[mydata$donor == "CK-038"] <- 43.4
mydata$PEV[mydata$donor == "CK-038"] <- 27.2

# ---- Anonymize donor IDs ----
Real_donor_ids <- sort(unique(mydata@meta.data$donor))
Crypto_ids <- sprintf("Donor_%02d", 1:length(Real_donor_ids))
names(Crypto_ids) <- Real_donor_ids

df <- mydata@meta.data
df$donor <- unlist(lapply(df$donor, function(x) Crypto_ids[x]))
mydata$donor <- df$donor

# ------------------------------------------------------------------------------
# Panel A: donor annotation strip
# ------------------------------------------------------------------------------
g3 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
  geom_tile(aes(fill = race), colour = "black") + theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("race") +
  scale_fill_manual(values = c("#FFFE99", "#F5C085", "#BDADD4", "#7FC97F"), na.value = NA) +
  theme(axis.title.x = element_text(angle = 90, vjust = .5, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

g4 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
  geom_tile(aes(fill = age), colour = "black") + theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("age") +
  theme(axis.title.x = element_text(angle = 90, vjust = .5, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_gradient(low = "ivory", high = "darkgoldenrod")

g5 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
  geom_tile(aes(fill = CL), colour = "black") + theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("CL") +
  theme(axis.title.x = element_text(angle = 90, vjust = .5, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_manual(values = c("#CAB2D6", "#FDB462", "#8DD3C7", "#FFFFB3",
                                "#BEBADA", "#FB8072", "#CCEBC5", "#FFED6F"))

g6 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
  geom_tile(aes(fill = BL), colour = "black") + theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("BL") +
  theme(axis.title.x = element_text(angle = 90, vjust = .5, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_manual(values = c("#FFB3BA", "#B3E2CD", "#FFFFBA"), na.value = NA)

g7 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
  geom_tile(aes(fill = as.factor(mRS)), colour = "black") + theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("mRS") +
  scale_fill_manual(values = c("#4D4D4D", "#999999", "#E0E0E0", "#F9DAC7", "#EF8962", "#B3242B"), na.value = NA) +
  theme(axis.title.x = element_text(angle = 90, vjust = .5, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

g8 <- ggplot(mydata@meta.data, aes(x = donor, y = 0.2)) +
  geom_tile(aes(fill = dx), colour = "black") + theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("dx") +
  theme(axis.title.x = element_text(angle = 90, vjust = .5, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_manual(values = c("#FDBF6F", "#B2DF8A", "#A6CEE3"), na.value = NA)

g9 <- ggplot(mydata@meta.data, aes(x = donor)) + geom_bar() + theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  coord_flip()

avg_nGenes_df <- mydata@meta.data %>%
  group_by(donor) %>%
  summarise(avg_nGenes = mean(nFeature_RNA, na.rm = TRUE))

g1 <- ggplot(avg_nGenes_df, aes(x = donor, y = 0.2, fill = avg_nGenes)) +
  geom_tile(colour = "black") + theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab("nGenes") +
  theme(axis.title.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_gradient(low = "white", high = "darkgreen", name = "Avg nGenes")

avg_nCounts_df <- mydata@meta.data %>%
  group_by(donor) %>%
  summarise(avg_nCounts = mean(nCount_RNA, na.rm = TRUE))

g2 <- ggplot(avg_nCounts_df, aes(x = donor, y = 0.2, fill = avg_nCounts)) +
  geom_tile(colour = "black") + theme_bw() + coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("nCounts") +
  theme(axis.title.x = element_text(angle = 90, vjust = .5, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_gradient(low = "white", high = "darkblue", name = "Avg nCounts")

sfig1A <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 +
  plot_layout(guides = "collect", widths = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.6)) &
  theme(legend.position = "right")

sfig1A <- sfig1A & theme(legend.position = "right", legend.text = element_text(size = 5))

save_plot_safe(sfig1A, "SFigure1_A_donor_annotations", width = 8, height = 6)

# ------------------------------------------------------------------------------
# Panel B: demographic pie charts
# ------------------------------------------------------------------------------
Categories <- c("sex", "race", "BL", "CL", "dx", "mRS")

# Full per-donor clinical/demographic table (S.Table 1) — superset of
# Categories (used below for the pie charts) plus the remaining fields shown
# across the donor annotation strip and reported elsewhere in the manuscript.
STable1_cols <- c("orig.ident", "donor", "sex", "race", "age", "BL", "CL", "dx",
                   "mRS", "mRS_binned", "mRS_3class", "TSH", "PHV", "PEV")
df_unique <- df[, STable1_cols] %>%
  distinct(donor, .keep_all = TRUE)
df_unique_ord <- df_unique[order(df_unique$donor), ]
rownames(df_unique_ord) <- df_unique_ord$donor

save_table_safe(df_unique_ord, "SFigure1_STable1_donor_metadata.csv")

create_pie_chart <- function(category) {
  category_data <- df_unique_ord[[category]]
  summary_data <- as.data.frame(table(category_data))
  colnames(summary_data) <- c(category, "Count")
  summary_data$Percentage <- round(100 * summary_data$Count / sum(summary_data$Count))

  p <- ggplot(summary_data, aes(x = "", y = Count, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    geom_text(aes(label = paste0(Percentage, "%")), position = position_stack(vjust = 0.5), size = 3) +
    labs(title = paste("Pie Chart:", category), fill = category) +
    theme_void() +
    theme(legend.position = "right", legend.text = element_text(size = 6))

  if (category == "sex")  p <- p + scale_fill_manual(values = c("#E69F00", "#56B4E9"), na.value = NA)
  if (category == "race") p <- p + scale_fill_manual(values = c("#FFFE99", "#F5C085", "#BDADD4", "#7FC97F"), na.value = NA)
  if (category == "dx")   p <- p + scale_fill_manual(values = c("#E78AC3", "#8DA1CB", "#EF8C61"), na.value = NA)
  if (category == "mRS")  p <- p + scale_fill_manual(values = c("#4D4D4D", "#999999", "#E0E0E0", "#F9DAC7", "#EF8962", "#B3242B"), na.value = NA)
  p
}

plots <- lapply(Categories, create_pie_chart)
p12 <- plots[[1]] + plots[[2]]
p34 <- plots[[3]] + plots[[4]]
p56 <- plots[[5]] + plots[[6]]
sfig1B <- p12 / p34 / p56

save_plot_safe(sfig1B, "SFigure1_B_demographic_pies", width = 8, height = 6)
