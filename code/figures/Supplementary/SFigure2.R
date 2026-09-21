#!/usr/bin/env Rscript
# ============================================================================
# SFigure2_FINAL.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================
#
# Extended Data Figure 2, panels c-f, and Supplementary Table 6, regenerated end
# to end from the myeloid Seurat object.
#
#   c  Mean z-scored expression of FreshMG subtype markers, 13 myeloid subtypes
#   d  Subtype dendrogram + crumblr effect size of proportion on six-month mRS
#   e  Subtype-level dreamlet DEGs associated with outcome (dotplot)
#   f  GO Molecular Function terms per subtype, split Good / Poor
#
# Self-contained. Reads the Seurat object and nothing else: no saved model fits,
# no side R libraries, no files from earlier runs. Everything it writes goes
# into ./SFigure2_FINAL/ and nothing outside that folder is touched.
#
# At the end it prints the numbers for the manuscript and a Methods sentence
# with the package versions of the run filled in, ready to paste.
#
#   Rscript SFigure2_FINAL.R
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr);                library(Seurat)
  library(SingleCellExperiment); library(dreamlet)
  library(crumblr);              library(variancePartition)
  library(ggplot2);              library(ggtree)
  library(aplot);                library(cowplot)
  library(reshape2);             library(broom)
  library(metafor);              library(ape)
  library(clusterProfiler);      library(org.Hs.eg.db)
})

# AnnotationDbi exports S4 generics that mask the dplyr verbs.
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename

set.seed(12345)
options(future.globals.maxSize = 100000000 * 1024^2)

# ================================================================== CONFIG ==
OBJ_PATH <- "Myeloid_Pred_ADAM.rds.zstd"
UTILS    <- "workflow/scripts/utils.R"

OUT_ROOT  <- "SFigure2"
N_THREADS <- 20
CONTRAST  <- c(Diff_B_vs_G = "mRS_binnedB-mRS_binnedG")

# Panel d, subtype composition.
MODEL_CRUMBLR  <- ~ (1 | race) + (1 | sex) + (1 | BL) + (1 | dx) + age + TSH + mRS_binned + 0
# Pseudobulk normalisation for panels e/f.
MODEL_PROCESS  <- ~ race + dx + BL + TSH + age + mRS_binned
# Panels e and f, differential expression. Same model as the cell-type analysis.
MODEL_DREAMLET <- ~ (1 | race) + (1 | dx) + (1 | BL) + age + mRS_binned + 0

MIN_COUNT     <- 1
PANEL_E_FDR   <- 0.05    # adj.P.Val cutoff for the dotplot gene set
PANEL_E_LFC   <- 1       # |avg_LFC| cutoff for the dotplot gene set
PANEL_F_FDR   <- 0.05    # adj.P.Val cutoff for the enrichment input
PANEL_F_LFC   <- 0.5     # |avg_LFC| cutoff for the enrichment input
PANEL_F_SHOWN <- 2       # GO terms displayed per cluster

SUBTYPE_ORDER <- c("Adapt_AIF1","Adapt_CCL3","Adapt_HIF1A","Adapt_HIST",
                   "Adapt_HSPA1A","Adapt_IFI44L","Adapt_TMEM163",
                   "Homeo_CECR2","Homeo_PICALM","PVM_CD163","ADAM_GPNMB",
                   "Prolif_MKI67","exAM_ERN1")

MARKERS <- c("AIF1","TMSB10","RPLP2","FTL","S100A9","CCL3","CCL4","CCL3L1",
             "CCL4L2","HIF1A","CREM","TNFRSF1B","SIPA1L1","FAM110B",
             "H1-3","H2AC8","H2BC4","H2BC7",
             "HSPA1A","HSPH1","HSP90AA1","DNAJA4","IFIT1","IFIT2","IFIT3",
             "IFI44","TMEM163","NCK2","SH3RF3","PTPN2","LIMK2","FRMD4A",
             "P2RY12","RASGEF1C","CX3CR1","NAV2","CECR2","PICALM","OXR1",
             "WNT2B","ELMO1","CD163","CSGALNACT1","SELENOP","GPNMB","PTPRG",
             "SAMD4A","MKI67","HELLS","CLSPN","CENPK","ERN1","CSKMT","PLK2")

# Current HGNC symbol -> the older symbol this object carries. Without this,
# AverageExpression() drops the four histone genes silently and the Adapt_HIST
# row of panel c has no defining markers.
ALIASES <- c(`H1-3` = "HIST1H1D", H2AC8 = "HIST1H2AE", H2BC4 = "HIST1H2BC",
             H2BC7 = "HIST1H2BF", SELENOP = "SEPP1")

# =================================================================== SETUP ==
FIG_DIR     <- file.path(OUT_ROOT, "figures")
TABLE_DIR   <- file.path(OUT_ROOT, "tables")
DERIVED_DIR <- file.path(OUT_ROOT, "derived")
for (d in c(FIG_DIR, TABLE_DIR, DERIVED_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

OUT_ABS <- normalizePath(OUT_ROOT, mustWork = TRUE)
out_file <- function(dir, filename) {
  p <- file.path(dir, filename)
  if (!startsWith(normalizePath(dirname(p), mustWork = FALSE), OUT_ABS))
    stop("Refusing to write outside ", OUT_ROOT, ": ", p)
  p
}
# Expensive fits are cached so an interrupted run resumes. The key includes a
# tag from the formula, so changing a model can never reuse an old fit.
model_tag <- function(f) {
  v <- utf8ToInt(paste(deparse(f), collapse = ""))
  sprintf("%06x", sum(v * seq_along(v)) %% 16777216L)
}
cached <- function(name, expr, model = NULL) {
  if (!is.null(model)) name <- paste0(name, "_", model_tag(model))
  f <- file.path(DERIVED_DIR, paste0(name, ".rds"))
  if (file.exists(f)) { message("[cache] reusing ", name); return(readRDS(f)) }
  val <- force(expr); saveRDS(val, f); val
}
save_both <- function(plot, filename, width, height) {
  pdf(out_file(FIG_DIR, paste0(filename, ".pdf")), width = width, height = height)
  print(plot); dev.off()
  png(out_file(FIG_DIR, paste0(filename, ".png")), width = width, height = height,
      units = "in", res = 300)
  print(plot); dev.off()
  message("[saved] ", filename)
}

VERS <- vapply(c("dreamlet","crumblr","variancePartition","clusterProfiler",
                 "org.Hs.eg.db","GO.db","Seurat","limma"),
               function(p) tryCatch(as.character(packageVersion(p)),
                                    error = function(e) "not installed"), character(1))
cat("\n================== ENVIRONMENT ==================\n")
cat(sprintf("  R %s\n", getRversion()))
for (p in names(VERS)) cat(sprintf("  %-18s %s\n", p, VERS[[p]]))
cat("=================================================\n")

source(UTILS)

# ================================================================== OBJECT ==
message("\n[1/6] loading object")
obj <- readCRDS(OBJ_PATH)

# The gene is CECR2 (Cat Eye Syndrome Chromosome Region 2). Earlier versions of
# this analysis labelled the subtype FRMD4A, then CERC2.
obj$Subtype[obj$Subtype == "Homeo_FRMD4A"] <- "Homeo_CECR2"
obj$Subtype[obj$Subtype == "Homeo_CERC2"]  <- "Homeo_CECR2"
stopifnot(setequal(unique(obj$Subtype), SUBTYPE_ORDER))
message("  cells: ", ncol(obj), " | donors: ", length(unique(obj$donor)),
        " | subtypes: ", length(unique(obj$Subtype)))

# ================================================================= PANEL C ==
message("\n[2/6] panel c - marker z-score heatmap")
resolve_symbols <- function(symbols, object) {
  present <- rownames(object); out <- character(0); miss <- character(0)
  for (s in symbols) {
    if (s %in% present) { out <- c(out, s); next }
    alt <- if (s %in% names(ALIASES)) ALIASES[[s]] else NA_character_
    if (!is.na(alt) && alt %in% present) {
      message("  [alias] ", s, " -> ", alt); out <- c(out, alt); next
    }
    miss <- c(miss, s)
  }
  if (length(miss))
    warning("Markers absent from the object, not plotted: ",
            paste(miss, collapse = ", "), call. = FALSE, immediate. = TRUE)
  unique(out)
}
markers_ok <- resolve_symbols(MARKERS, obj)
message("  resolved ", length(markers_ok), "/", length(MARKERS), " markers")
writeLines(setdiff(MARKERS, markers_ok),
           out_file(DERIVED_DIR, "panelC_missing_symbols.txt"))

avg <- AverageExpression(obj, group.by = "Subtype", assays = "RNA",
                         slot = "data", features = markers_ok)$RNA
z <- scale(t(as.matrix(avg)))
rownames(z) <- gsub("-", "_", rownames(z))   # AverageExpression rewrites "_" as "-"
stopifnot(all(SUBTYPE_ORDER %in% rownames(z)))
z <- z[SUBTYPE_ORDER, , drop = FALSE]

zm <- melt(z, varnames = c("subtype", "gene"))
zm$subtype <- factor(zm$subtype, levels = rev(SUBTYPE_ORDER))
zm$gene    <- factor(zm$gene, levels = colnames(z))
p_c <- ggplot(zm, aes(gene, subtype, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3",
                       midpoint = 0, name = "Mean Z-Score",
                       guide = guide_colorbar(frame.linewidth = 0.7,
                                              frame.colour = "black",
                                              barwidth = 0.8, barheight = 10)) +
  labs(title = "FreshMG markers", x = "", y = "") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title  = element_text(hjust = 0.5))
save_both(p_c, "SFigure2_C_marker_zscore_heatmap", 15, 12)
write.csv(as.data.frame(z), out_file(TABLE_DIR, "SFigure2_C_marker_zscores.csv"))

# ============================================================== PSEUDOBULK ==
message("\n[3/6] pseudobulk aggregation")
sce <- as.SingleCellExperiment(obj)
pb  <- cached("pseudobulk_subtype",
              aggregateToPseudoBulk(sce, assay = "counts", cluster_id = "Subtype",
                                    sample_id = "donor",
                                    BPPARAM = SnowParam(5, progressbar = TRUE)))

# ================================================================= PANEL D ==
message("\n[4/6] panel d - crumblr composition vs mRS")
res_d <- cached("crumblr", {
  cobj <- crumblr(cellCounts(pb))
  L    <- makeContrastsDream(MODEL_CRUMBLR, colData(pb),
                             contrasts = c(mRS_contrast = "mRS_binnedB-mRS_binnedG"))
  fit  <- eBayes(dream(cobj, MODEL_CRUMBLR, colData(pb), L = L))
  tt   <- topTable(fit, coef = "mRS_contrast", number = Inf, sort.by = "none")
  tt$assay <- rownames(tt); tt
}, MODEL_CRUMBLR)

meta_d <- res_d %>% dplyr::as_tibble() %>% dplyr::group_by(assay) %>%
  dplyr::do(tidy(rma(yi = logFC, sei = logFC / t, data = ., method = "FE"))) %>%
  dplyr::select(-term, -type) %>%
  dplyr::mutate(FDR = p.adjust(p.value, "fdr")) %>% as.data.frame()
write.csv(meta_d, out_file(TABLE_DIR, "SFigure2_D_crumblr_mRS.csv"), row.names = FALSE)
sig_d <- meta_d$assay[meta_d$FDR < 0.05]
message("  subtypes at FDR < 0.05: ", paste(sig_d, collapse = ", "))

hc        <- buildClusterTreeFromPB(pb)
tree_base <- ggtree(as.phylo(hc), branch.length = "none") +
  geom_tiplab(colour = "black", size = 4, hjust = 0, offset = 0.2)
fig_tree  <- tree_base + xlim(0, layer_scales(tree_base)$x$range$range[2] * 2.2)

meta_d$celltype <- factor(meta_d$assay, rev(ggtree::get_taxa_name(fig_tree)))
p_d <- ggplot(meta_d, aes(celltype, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey", linewidth = 1) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error), width = 0) +
  geom_point(aes(colour = pmin(4, -log10(FDR)), size = pmin(4, -log10(FDR)))) +
  scale_colour_gradient2(name = bquote(-log[10] ~ FDR), limits = c(0, 4),
                         low = "grey90", mid = "red", high = "darkred",
                         midpoint = -log10(0.01)) +
  scale_size_area(name = bquote(-log[10] ~ FDR), limits = c(0, 4)) +
  geom_text(aes(label = ifelse(FDR < 0.05, "+", "")),
            colour = "white", size = 6, vjust = 0.4, hjust = 0.5) +
  theme_classic() + coord_flip() + xlab("") + ylab("mRS") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12), text = element_text(size = 20)) +
  scale_y_continuous(breaks = scales::breaks_pretty(3))
save_both(insert_left(p_d, fig_tree, width = 1.4),
          "SFigure2_D_tree_crumblr_mRS", 8, 4)

# ================================================================= PANEL E ==
message("\n[5/6] panel e - subtype differential expression")
res_proc <- cached("processAssays",
                   processAssays(pb, MODEL_PROCESS, min.count = MIN_COUNT,
                                 BPPARAM = SnowParam(N_THREADS, progressbar = TRUE)),
                   MODEL_PROCESS)
res_dream <- cached("dreamlet",
                    dreamlet(res_proc, MODEL_DREAMLET, contrasts = CONTRAST,
                             BPPARAM = SnowParam(N_THREADS, progressbar = TRUE)),
                    MODEL_DREAMLET)
cat("\n--- model fitted per assay ---\n"); print(details(res_dream)); cat("\n")

DEGs <- as.data.frame(topTable(res_dream, p.value = 1, lfc = 0,
                               number = 1e8, coef = names(CONTRAST)))
DEGs <- DEGs[order(DEGs$adj.P.Val), ]
n_assays <- length(unique(DEGs$assay))
DEGs <- DEGs %>% dplyr::group_by(ID) %>%
  dplyr::mutate(pos_avgLFC = mean(logFC),
                avg_LFC    = sum(logFC) / n_assays) %>%
  dplyr::ungroup() %>% as.data.frame()
write.csv(DEGs, out_file(TABLE_DIR, "SupplementaryTable6_DEGs_myeloid_subtypes.csv"),
          row.names = FALSE)

n_pairs <- sum(DEGs$adj.P.Val < 0.05, na.rm = TRUE)
n_genes <- length(unique(DEGs$ID[DEGs$adj.P.Val < 0.05]))

plot_genes <- unique(DEGs$ID[DEGs$adj.P.Val < PANEL_E_FDR &
                             abs(DEGs$avg_LFC) > PANEL_E_LFC])
dp <- na.omit(DEGs[DEGs$ID %in% plot_genes, ])
dp <- dp[order(dp$avg_LFC), ]
dp$gene  <- factor(as.vector(dp$ID), levels = unique(as.vector(dp$ID)))
dp$assay <- factor(dp$assay, levels = rev(SUBTYPE_ORDER))
message("  panel e genes: ", length(plot_genes))

p_e <- ggplot(dp, aes(x = gene, y = assay, size = abs(logFC), fill = logFC)) +
  geom_point(shape = 21, colour = "black") +
  guides(size = guide_legend(title = "log2FC"),
         fill = guide_legend(title = "-log2FC")) +
  theme_cowplot() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.title = element_blank(),
        axis.text  = element_text(size = 11),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.position = "bottom")
save_both(p_e, "SFigure2_E_subtype_DEG_dotplot", 10, 5)

# ================================================================= PANEL F ==
# Deliberately not cached: a cache key that does not encode the gene set can
# silently return an enrichment computed from different input.
message("\n[6/6] panel f - GO:MF enrichment by subtype and direction")
DEGs$cluster <- ifelse(DEGs$avg_LFC > 0, "Poor", "Good")
DEGs$absLFC  <- abs(DEGs$avg_LFC)
DEGs$gene    <- DEGs$ID
sig_f <- DEGs[DEGs$adj.P.Val <= PANEL_F_FDR & DEGs$absLFC >= PANEL_F_LFC, ]
message("  enrichment input: ", nrow(sig_f), " rows, ",
        length(unique(sig_f$gene)), " genes, ",
        length(unique(sig_f$assay)), " subtypes")

MF_res <- compareCluster(gene ~ cluster + assay, data = sig_f,
                         fun = "enrichGO", ont = "MF", keyType = "SYMBOL",
                         OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05)
write.csv(as.data.frame(MF_res),
          out_file(TABLE_DIR, "SFigure2_F_SourceData_GOMF.csv"), row.names = FALSE)

MF_top <- dotplot(MF_res, showCategory = PANEL_F_SHOWN)$data
MF_top$Significant <- ifelse(MF_top$qvalue < 0.05, "+", "")
MF_top$directed    <- -log10(MF_top$qvalue) * ifelse(MF_top$cluster == "Good", -1, 1)
MF_top$assay       <- gsub("Adapt_|Homeo_|PVM_|Prolif_|ADAM_|exAM_", "", MF_top$assay)
write.csv(MF_top, out_file(TABLE_DIR, "SFigure2_F_plotted_terms.csv"), row.names = FALSE)

p_f <- ggplot(MF_top, aes(x = assay, y = Description, colour = directed,
                          fill = directed, label = Significant)) +
  geom_tile() +
  geom_text(alpha = 1, size = 5, colour = "white") +
  facet_wrap(~ cluster) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3",
                       na.value = "white") +
  scale_colour_gradient2(low = "navy", mid = "white", high = "firebrick3",
                         na.value = "white") +
  theme_cowplot() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 15)) +
  xlab("") + ylab("")
save_both(p_f, "SFigure2_F_GOMF_good_vs_poor", 15, 10)

# ============================================== CHECKS AGAINST THE TEXT =====
good_terms <- sort(unique(as.character(MF_top$Description[MF_top$cluster == "Good"])))
poor_terms <- sort(unique(as.character(MF_top$Description[MF_top$cluster == "Poor"])))
has <- function(v, pat) any(grepl(pat, v, ignore.case = TRUE))

cat("\n############ CHECKS AGAINST THE MANUSCRIPT ############\n")
cat(sprintf("  VEGFB in panel e gene set                 : %s\n",
            "VEGFB" %in% plot_genes))
cat(sprintf("  chemokine / cytokine terms in Good        : %s\n",
            has(good_terms, "chemokine|cytokine")))
cat(sprintf("  chromatin / histone terms in Poor         : %s\n",
            has(poor_terms, "chromatin|histone|methylat")))
cat(sprintf("  CECR2, TMEM163, CD163 significant in d    : %s\n",
            all(c("Homeo_CECR2","Adapt_TMEM163","PVM_CD163") %in% sig_d)))
cat("#######################################################\n")

cat("\nGood terms (", length(good_terms), "):\n", sep = "")
for (t in good_terms) cat("   ", t, "\n", sep = "")
cat("\nPoor terms (", length(poor_terms), "):\n", sep = "")
for (t in poor_terms) cat("   ", t, "\n", sep = "")

# =================================================== FOR THE MANUSCRIPT =====
cat("\n############### NUMBERS FOR THE MANUSCRIPT ###############\n")
cat(sprintf("  unique genes, adjusted p < 0.05 : %d\n", n_genes))
cat(sprintf("  subtype-gene associations       : %d\n", n_pairs))
cat(sprintf("  genes plotted in panel e        : %d\n", length(plot_genes)))
cat("##########################################################\n")

methods <- sprintf(paste0(
  "Gene Ontology Molecular Function over-representation for Extended Data ",
  "Fig. 2f was performed with clusterProfiler v%s against org.Hs.eg.db v%s ",
  "and GO.db v%s, using genes at adjusted p < 0.05 and |average log2 fold ",
  "change| >= 0.5 within each subtype and outcome direction, displaying the ",
  "top two terms per group."),
  VERS[["clusterProfiler"]], VERS[["org.Hs.eg.db"]], VERS[["GO.db"]])

cat("\n-------- METHODS SENTENCE, READY TO PASTE --------\n\n")
cat(methods, "\n\n")
cat("--------------------------------------------------\n")

writeLines(c(sprintf("unique_genes\t%d", n_genes),
             sprintf("pairs\t%d", n_pairs),
             sprintf("panel_e_genes\t%d", length(plot_genes)),
             "", "METHODS SENTENCE:", methods),
           out_file(DERIVED_DIR, "for_the_manuscript.txt"))

message("\nDone. Figures in ", FIG_DIR, " | tables in ", TABLE_DIR)
sessionInfo()
