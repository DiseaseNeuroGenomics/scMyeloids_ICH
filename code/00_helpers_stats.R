# ============================================================================
# 00_helpers_stats.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Shared statistical / plotting helper functions used by both Figure 1 and
# Figure 2 scripts (mRS coefficient plots, meta-analysis, tree plots).
# Source after 00_setup.R.
# ==============================================================================

# ------------------------------------------------------------------------------
# Fixed-effects meta-analysis across a list of topTable-style data frames
# ------------------------------------------------------------------------------
dh_meta_analysis <- function(tabList) {
  message(sprintf("[%s] INFO: Starting fixed-effects meta-analysis across %d datasets", Sys.time(), length(tabList)))

  if (is.null(names(tabList))) {
    names(tabList) <- as.character(seq(length(tabList)))
  }
  for (key in names(tabList)) {
    tabList[[key]]$Dataset <- key
  }
  df <- do.call(rbind, tabList)

  res <- df %>%
    as_tibble() %>%
    group_by(assay) %>%
    do(tidy(rma(yi = logFC, sei = logFC / t, data = ., method = "FE"))) %>%
    select(-term, -type) %>%
    mutate(
      FDR = p.adjust(p.value, "fdr"),
      log10FDR = -log10(FDR)
    )

  message(sprintf("[%s] INFO: Meta-analysis complete.", Sys.time()))
  res
}

# ------------------------------------------------------------------------------
# Hierarchical cluster tree (used as the left-hand dendrogram next to
# coefficient/forest plots)
# ------------------------------------------------------------------------------
dh_plotTree <- function(tree, low = "grey90", mid = "red", high = "darkred", xmax.scale = 1.5) {
  fig <- ggtree(tree, branch.length = "none") +
    geom_tiplab(color = "black", size = 4, hjust = 0, offset = 0.2) +
    theme(legend.position = "top left", plot.title = element_text(hjust = 0.5))

  xmax <- layer_scales(fig)$x$range$range[2]
  fig + xlim(0, xmax * xmax.scale)
}

# ------------------------------------------------------------------------------
# Forest/coefficient plot (effect size +/- 95% CI per cell type), colored and
# sized by -log10(FDR). Set label_axis = FALSE when pairing with a tree that
# already shows cell-type labels (as in Figure 2 panel C).
# ------------------------------------------------------------------------------
dh_plotCoef <- function(tab, fig.tree, low = "grey90", mid = "red", high = "darkred",
                         ylab = "Effect Size (logFC)", label_axis = TRUE) {
  tab$logFC <- tab$estimate
  tab$celltype <- factor(tab$assay, rev(ggtree::get_taxa_name(fig.tree)))
  tab$se <- tab$std.error

  fig.es <- ggplot(tab, aes(x = celltype, y = logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +
    geom_errorbar(aes(ymin = logFC - 1.96 * se, ymax = logFC + 1.96 * se), width = 0) +
    geom_point(aes(color = pmin(4, -log10(FDR)), size = pmin(4, -log10(FDR)))) +
    scale_color_gradient2(
      name = bquote(-log[10]~FDR),
      limits = c(0, 4), low = low, mid = mid, high = high, midpoint = -log10(0.01)
    ) +
    scale_size_area(name = bquote(-log[10]~FDR), limits = c(0, 4)) +
    geom_text(
      aes(label = ifelse(FDR < 0.05, '+', '')),
      color = "white", size = 6, vjust = 0.4, hjust = 0.5
    ) +
    theme_classic() +
    coord_flip() +
    labs(x = "Cell Type", y = ylab) +
    theme(
      axis.text.y = if (label_axis) element_text(size = 12, color = "black") else element_blank(),
      axis.ticks.y = if (label_axis) element_line(color = "black") else element_blank(),
      axis.text.x = element_text(size = 12),
      text = element_text(size = 14)
    ) +
    scale_y_continuous(breaks = scales::breaks_pretty(3))

  fig.es
}

# ------------------------------------------------------------------------------
# Voom mean-variance diagnostic panel (QC for processAssays() output)
# ------------------------------------------------------------------------------
my_VroomPlot <- function(x, ncol = 3, alpha = .5, assays = names(x)) {
  y <- NULL

  assays <- intersect(assays, names(x))
  if (length(assays) == 0) stop("No valid assays selected")

  df_range <- lapply(assays, function(id) {
    if (!is.null(x[[id]]$voom.xy)) {
      with(x[[id]]$voom.xy, data.frame(range(x), range(y), id = id))
    } else {
      NULL
    }
  })
  df_range <- do.call(rbind, df_range)
  if (is.null(df_range)) stop("Voom was not run on this object")

  xlim <- range(df_range$range.x.)
  ylim <- c(0, max(df_range$range.y.))
  xlab <- bquote(log[2](counts + 0.5))
  ylab <- bquote(sqrt(standard~deviation))

  validAssays <- droplevels(factor(unique(df_range)$id, assays))

  df.list <- lapply(validAssays, function(id) with(x[[id]]$voom.xy, data.frame(id, x, y)))
  df_points <- do.call(rbind, df.list)
  df_points$id <- droplevels(factor(df_points$id, assays))
  df_points <- df_points[order(df_points$id), ]

  df.list <- lapply(validAssays, function(id) with(x[[id]]$voom.line, data.frame(id, x, y)))
  df_curve <- do.call(rbind, df.list)

  unique_ids <- unique(df_curve$id)
  plotlist_vroom <- lapply(unique_ids, function(sub_id) {
    df_point_sub <- df_points %>% filter(id == sub_id)
    df_curve_sub <- df_curve %>% filter(id == sub_id)

    ggplot(df_point_sub, aes(x, y)) +
      geom_point(size = 0.1, alpha = alpha) +
      theme_classic() +
      theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~id, ncol = ncol) +
      xlab(xlab) + ylab(ylab) + xlim(xlim) + ylim(ylim) +
      geom_line(data = df_curve_sub, aes(x, y), color = "red")
  })

  ggpubr::ggarrange(plotlist = plotlist_vroom, ncol = ncol, align = "hv")
}

# ------------------------------------------------------------------------------
# Nature Medicine-style ggplot theme, applied to every final panel before saving
# ------------------------------------------------------------------------------
theme_natmed <- function(base_family = "Helvetica", base_size = 11) {
  theme(
    text = element_text(family = base_family),
    plot.title = element_text(family = base_family, size = 12, face = "bold"),
    axis.title = element_text(family = base_family, size = base_size),
    axis.text.x = element_text(family = base_family, size = base_size - 1, color = "black"),
    axis.text.y = element_text(family = base_family, size = base_size - 1, color = "black", hjust = 1),
    legend.title = element_text(family = base_family, size = 10),
    legend.text = element_text(family = base_family, size = 9)
  )
}

init_natmed_fonts <- function() {
  # font_add() lives in sysfonts, not showtext's own namespace — showtext::font_add()
  # fails even though showtext normally re-exposes it on the search path via Depends.
  sysfonts::font_add(family = "Helvetica",
                      regular = "/usr/share/fonts/urw-base35/NimbusSans-Regular.otf",
                      bold    = "/usr/share/fonts/urw-base35/NimbusSans-Bold.otf",
                      italic  = "/usr/share/fonts/urw-base35/NimbusSans-Italic.otf")
  showtext::showtext_auto()
  theme_set(theme_natmed())
}
