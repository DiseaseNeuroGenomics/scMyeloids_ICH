# ============================================================================
# Figure4_FINAL.R
# Panels A-C: Mouse_Behavioral_Dimitris.R modeling, unchanged.
# Panels D-K: 2026_08_ImageAnalysis.R modeling, unchanged, except panel D
# labels now use BH p pooled across all 8 interaction terms.
# Outputs: Figure4_FINAL.pdf
#          Tables/Table18_Figure4_behaviour_FULL.csv
#          Tables/Table19_Figure4DE_morphology_FULL.csv
# ============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(LaplacesDemon)
library(patchwork)
library(ggsignif)

dir.create("Tables", recursive = TRUE, showWarnings = FALSE)

theme_natmed <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line        = element_line(linewidth = 0.3, color = "black"),
      axis.ticks       = element_line(linewidth = 0.3, color = "black"),
      axis.text        = element_text(color = "black", size = base_size),
      axis.title       = element_text(color = "black", size = base_size),
      plot.title       = element_text(size = base_size + 1, hjust = 0, face = "bold"),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size = base_size),
      legend.key.size  = unit(0.3, "cm"),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size, face = "bold"),
      panel.grid       = element_blank()
    )
}

treatment_cols <- c(Control = "#666666", EV = "#D55E00", RA = "#0072B2")
sig_label <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "n.s.")

## ============================================================================
## PART A -- behavioral panels, exact modeling from Mouse_Behavioral_Dimitris.R
## Two-group model per drug, separately. Now also collects n / AIC / convergence
## so we can emit one full supplementary table, not just panel labels.
## ============================================================================
behav_raw <- read.csv("data/Modified_ICH_Microglia_Experimental_Data.csv")
num_cols <- grep("^(Rotarod|MNSS|BW_|Dose)", names(behav_raw), value = TRUE)
for (cc in num_cols) behav_raw[[cc]] <- suppressWarnings(as.numeric(as.character(behav_raw[[cc]])))

measures <- c("RotarodRPMMedians", "RotarodLatencyMedians", "MNSS")
sfx_list <- list(c("BSL","Day3","Day7"), c("BSL","Day3","Day7"), c("2h","Day1","Day3","Day7"))
day_list <- list(c(0,3,7), c(0,3,7), c(0.083,1,3,7))
fam_list <- c("Gamma", "Gamma", "Poisson")
base_tp  <- c("BSL", "BSL", "2h")

behav_results <- data.frame()
behav_n       <- data.frame()
behav_bic     <- data.frame()
behav_conv    <- data.frame()

for (i in seq_along(measures)) {
  meas <- measures[i]; sfx <- sfx_list[[i]]; dys <- day_list[[i]]; fm <- fam_list[i]
  cols <- paste0(meas, "_", sfx)
  
  for (drug in c("RA", "Everolimus")) {
    sub <- behav_raw[behav_raw$Condition %in% c("Control", drug), ]
    sc <- c(); for (cc in cols) sc <- c(sc, sub[[cc]])
    
    long <- data.frame(
      Animal    = rep(sub$Animal.Code, length(sfx)),
      Condition = rep(sub$Condition,   length(sfx)),
      Dose      = rep(sub$Dose,        length(sfx)),
      Timepoint = rep(sfx, each = nrow(sub)),
      Days      = rep(dys, each = nrow(sub)),
      Score     = sc
    )
    long <- long[!is.na(long$Score), ]
    long$Condition <- factor(long$Condition, levels = c("Control", drug))
    long$Timepoint <- factor(long$Timepoint, levels = sfx, ordered = TRUE)
    long$Dose_s    <- as.numeric(scale(long$Dose))
    
    nt <- as.data.frame(table(long$Condition, long$Timepoint))
    names(nt) <- c("Condition", "Timepoint", "n")
    behav_n <- rbind(behav_n, data.frame(Measure = meas, Drug = drug, nt))
    
    ff <- if (fm == "Gamma") Gamma(link = "log") else poisson(link = "log")
    
    # ---- per-timepoint model ----
    m1 <- glmer(Score ~ Timepoint * Condition + Dose_s + (1 | Animal), family = ff, data = long)
    cm <- m1@optinfo$conv$lme4$messages; if (is.null(cm)) cm <- "OK"
    behav_conv <- rbind(behav_conv, data.frame(Measure = meas, Drug = drug, Model = "per_timepoint",
                                               convergence = paste(cm, collapse = " | ")))
    
    if (fm == "Gamma") {
      m_alt <- lmer(Score ~ Timepoint * Condition + Dose_s + (1 | Animal), data = long, REML = FALSE)
      behav_bic <- rbind(behav_bic, data.frame(Measure = meas, Drug = drug,
                                               family = c("Gamma(log)", "Gaussian"),
                                               AIC = round(c(AIC(m1), AIC(m_alt)), 1)))
    }
    
    ph <- emmeans(m1, pairwise ~ Condition | Timepoint)
    ct <- as.data.frame(summary(ph)$contrasts)
    zz <- if ("z.ratio" %in% names(ct)) ct$z.ratio else ct$t.ratio
    behav_results <- rbind(behav_results, data.frame(
      Measure = meas, Drug = drug, Test = "per_timepoint",
      Timepoint = as.character(ct$Timepoint),
      estimate = ct$estimate, SE = ct$SE, stat = zz, p = ct$p.value))
    
    # ---- trajectory (slope) model ----
    m2 <- glmer(Score ~ Days * Condition + Dose_s + (1 | Animal), family = ff, data = long)
    cm2 <- m2@optinfo$conv$lme4$messages; if (is.null(cm2)) cm2 <- "OK"
    behav_conv <- rbind(behav_conv, data.frame(Measure = meas, Drug = drug, Model = "trajectory",
                                               convergence = paste(cm2, collapse = " | ")))
    
    sl <- emtrends(m2, pairwise ~ Condition, var = "Days")
    st <- as.data.frame(summary(sl)$contrasts)
    zs <- if ("z.ratio" %in% names(st)) st$z.ratio else st$t.ratio
    behav_results <- rbind(behav_results, data.frame(
      Measure = meas, Drug = drug, Test = "trajectory", Timepoint = "slope_per_day",
      estimate = st$estimate, SE = st$SE, stat = zs, p = st$p.value))
  }
}

behav_results$role <- "efficacy"
for (i in seq_along(measures)) {
  behav_results$role[behav_results$Measure == measures[i] & behav_results$Timepoint == base_tp[i]] <- "balance check"
}
behav_results$role[behav_results$Test == "trajectory"] <- "trajectory"

behav_results$fold_drug_vs_vehicle <- round(exp(-behav_results$estimate), 3)

# BH within measure, per-timepoint efficacy rows only, both drugs pooled -- matches Mouse_Behavioral_Dimitris.R exactly
behav_results$p_BH <- NA
for (i in seq_along(measures)) {
  idx <- behav_results$Measure == measures[i] & behav_results$Test == "per_timepoint" & behav_results$role == "efficacy"
  behav_results$p_BH[idx] <- p.adjust(behav_results$p[idx], "BH")
}
behav_results$sig <- ifelse(!is.na(behav_results$p_BH) & behav_results$p_BH < 0.05, "*", "")

behav_results$Model <- ifelse(behav_results$Measure == "MNSS", "Poisson (log)", "Gamma (log)")
behav_results$n_vehicle <- NA; behav_results$n_drug <- NA
behav_results$AIC_Gamma <- NA; behav_results$AIC_Gaussian <- NA; behav_results$convergence <- NA

for (r in 1:nrow(behav_results)) {
  tp <- behav_results$Timepoint[r]
  if (behav_results$Test[r] == "trajectory") tp <- base_tp[match(behav_results$Measure[r], measures)]
  
  sel <- behav_n$Measure == behav_results$Measure[r] & behav_n$Drug == behav_results$Drug[r] & behav_n$Timepoint == tp
  nv <- behav_n$n[sel & behav_n$Condition == "Control"]
  nd <- behav_n$n[sel & behav_n$Condition == behav_results$Drug[r]]
  if (length(nv) == 1) behav_results$n_vehicle[r] <- nv
  if (length(nd) == 1) behav_results$n_drug[r]    <- nd
  
  selb <- behav_bic$Measure == behav_results$Measure[r] & behav_bic$Drug == behav_results$Drug[r]
  if (any(selb)) {
    behav_results$AIC_Gamma[r]    <- behav_bic$AIC[selb & behav_bic$family == "Gamma(log)"][1]
    behav_results$AIC_Gaussian[r] <- behav_bic$AIC[selb & behav_bic$family == "Gaussian"][1]
  }
  
  selc <- behav_conv$Measure == behav_results$Measure[r] & behav_conv$Drug == behav_results$Drug[r] &
    behav_conv$Model == behav_results$Test[r]
  if (any(selc)) behav_results$convergence[r] <- behav_conv$convergence[selc][1]
}

behav_results$Comparison <- paste(behav_results$Drug, "vs Vehicle")
behav_results$Measure <- recode(behav_results$Measure,
                                RotarodRPMMedians     = "Rotarod RPM (median)",
                                RotarodLatencyMedians = "Rotarod latency to fall (median)",
                                MNSS                  = "mNSS")

behav_table <- behav_results[, c("Measure","Comparison","Test","Timepoint","role",
                                 "n_vehicle","n_drug","estimate","SE","stat",
                                 "fold_drug_vs_vehicle","p","p_BH","sig",
                                 "Model","AIC_Gamma","AIC_Gaussian","convergence")]
names(behav_table) <- c("Measure","Comparison","Test","Timepoint","Role",
                        "n_vehicle","n_drug","Estimate_log","SE","Statistic",
                        "Fold_change_vs_vehicle","P_raw","P_BH","Significant",
                        "Model","AIC_Gamma","AIC_Gaussian","Convergence")

write.csv(behav_table, "Tables/Table18_Figure4_behaviour_FULL.csv", row.names = FALSE)

## ---- plotting: panel labels now show Day3, Day7, AND slope, per drug ----

make_behavior_panel <- function(score_prefix, timepoint_levels, timepoint_recode, y_label, panel_label) {
  plot_df <- behav_raw %>%
    mutate(Treatment = recode(Condition, "Everolimus" = "EV"),
           Treatment = factor(Treatment, levels = c("Control","EV","RA"))) %>%
    filter(Treatment %in% c("Control","EV","RA")) %>%
    select(Animal.Code, Treatment, starts_with(score_prefix)) %>%
    pivot_longer(cols = starts_with(score_prefix), names_to = "Timepoint", values_to = "Score") %>%
    mutate(Timepoint = recode(Timepoint, !!!timepoint_recode),
           Timepoint = factor(Timepoint, levels = timepoint_levels, ordered = TRUE)) %>%
    filter(!is.na(Timepoint) & !is.na(Score))
  
  meas_full <- behav_results$Measure[1]  # placeholder, overwritten below by exact match
  meas_key <- recode(score_prefix,
                     RotarodRPMMedians = "Rotarod RPM (median)",
                     RotarodLatencyMedians = "Rotarod latency to fall (median)",
                     MNSS = "mNSS")
  eff_tps <- setdiff(timepoint_levels, timepoint_levels[1])  # everything except baseline
  
  fmt_line <- function(drug) {
    parts <- c()
    for (tp in eff_tps) {
      row <- behav_results[behav_results$Measure == meas_key & behav_results$Drug == drug &
                             behav_results$Test == "per_timepoint" & behav_results$Timepoint == tp, ]
      if (nrow(row) == 1) parts <- c(parts, paste0(tp, ": p=", format(row$p_BH, digits = 2), " ", sig_label(row$p_BH)))
    }
    slope_row <- behav_results[behav_results$Measure == meas_key & behav_results$Drug == drug &
                                 behav_results$Test == "trajectory", ]
    if (nrow(slope_row) == 1) parts <- c(parts, paste0("slope: p=", format(slope_row$p, digits = 2), " ", sig_label(slope_row$p)))
    paste(parts, collapse = ", ")
  }
  
  label_text <- paste0(
    "EV vs Control (BH): ", fmt_line("Everolimus"), "\n",
    "RA vs Control (BH): ", fmt_line("RA")
  )
  
  ggplot(plot_df, aes(Timepoint, Score, color = Treatment, group = Treatment)) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, linewidth = 0.4) +
    stat_summary(fun = mean, geom = "point", size = 1.4) +
    annotate("text", x = -Inf, y = Inf, label = label_text,
             hjust = -0.05, vjust = 1.2, size = 1.9, color = "black", lineheight = 1.1) +
    scale_color_manual(values = treatment_cols) +
    labs(x = NULL, y = y_label, title = panel_label) +
    theme_natmed()
}

p_rpm  <- make_behavior_panel("RotarodRPMMedians", c("BSL","Day3","Day7"),
                              c(RotarodRPMMedians_BSL="BSL", RotarodRPMMedians_Day3="Day3", RotarodRPMMedians_Day7="Day7"),
                              "Rotarod RPM (median)", "a   Rotarod RPM")

p_lat  <- make_behavior_panel("RotarodLatencyMedians", c("BSL","Day3","Day7"),
                              c(RotarodLatencyMedians_BSL="BSL", RotarodLatencyMedians_Day3="Day3", RotarodLatencyMedians_Day7="Day7"),
                              "Rotarod latency (s, median)", "b   Rotarod latency")

p_mnss <- make_behavior_panel("MNSS", c("2h","Day1","Day3","Day7"),
                              c(MNSS_2h="2h", MNSS_Day1="Day1", MNSS_Day3="Day3", MNSS_Day7="Day7"),
                              "mNSS score", "c   mNSS")


## ============================================================================
## PART B -- morphology panels, exact modeling from 2026_08_ImageAnalysis.R
## Now exports every coefficient per morphotype (not just the 2 interaction
## terms) plus the windowed boxplot table, into one combined supplementary CSV.
## ============================================================================

morph_file <- "data/all_samples_clustered.csv.gz"
df_full <- read_csv(morph_file, show_col_types = FALSE) %>%
  separate_wider_delim(cols = sample, delim = "_", names = c("Animal","Treatment","Section","Other")) %>%
  mutate(cluster_name = factor(cluster_name, c("Ameboid","Ramified","Hypertrophic","Rod-like"))) %>%
  mutate(dist = dist_um) %>%
  filter(dist < 1500)

celltypes <- c("Ameboid", "Hypertrophic", "Ramified", "Rod-like")
cutoffs   <- c(500, 1000, 1500)

fit_ct <- function(CT) {
  form <- (cluster_name == CT) ~ dist * Treatment + (1 | Animal/Section)
  fit  <- glmmTMB(form, df_full, family = binomial())
  co <- coef(summary(fit))$cond %>% data.frame() %>%
    rownames_to_column("coef") %>% rename(se = Std..Error, pvalue = Pr...z..) %>%
    mutate(CellType = CT)
  list(fit = fit, co = co)
}

fits <- setNames(lapply(celltypes, fit_ct), celltypes)

# full coefficient table, every term, every morphotype
morph_trend_full <- bind_rows(lapply(fits, function(x) x$co))

# BH across the 8 interaction terms only (this is the family Results text describes)
is_interaction <- morph_trend_full$coef %in% c("dist:TreatmentEV", "dist:TreatmentRA")
morph_trend_full$p_BH_interaction_family <- NA
morph_trend_full$p_BH_interaction_family[is_interaction] <- p.adjust(morph_trend_full$pvalue[is_interaction], "BH")
morph_trend_full$PanelType <- "distance_trend (panel D)"
morph_trend_full$Cutoff <- NA

interaction_p <- morph_trend_full %>%
  filter(coef %in% c("dist:TreatmentEV","dist:TreatmentRA")) %>%
  mutate(Drug = ifelse(coef == "dist:TreatmentEV", "EV", "RA")) %>%
  select(CellType, Drug, p_BH = p_BH_interaction_family)

make_trend_panel <- function(CT, panel_label) {
  fit <- fits[[CT]]$fit
  p_ev <- interaction_p$p_BH[interaction_p$CellType == CT & interaction_p$Drug == "EV"]
  p_ra <- interaction_p$p_BH[interaction_p$CellType == CT & interaction_p$Drug == "RA"]
  
  label_text <- paste0(
    "EV x dist (BH): p=", format(p_ev, digits = 2), " ", sig_label(p_ev), "\n",
    "RA x dist (BH): p=", format(p_ra, digits = 2), " ", sig_label(p_ra)
  )
  
  pred <- expand.grid(dist = seq(10, 1500, length.out = 300), Treatment = c("Control","EV","RA"))
  X <- model.matrix(~ dist * Treatment, pred)
  b <- fixef(fit)$cond
  S <- vcov(fit)$cond
  pred$eta    <- as.numeric(X %*% b)
  pred$eta.se <- sqrt(rowSums((X %*% S) * X))
  pred$mu     <- invlogit(pred$eta)
  pred$lo     <- invlogit(pred$eta - 1.96 * pred$eta.se)
  pred$hi     <- invlogit(pred$eta + 1.96 * pred$eta.se)
  
  ggplot(pred, aes(dist, mu, color = Treatment, fill = Treatment)) +
    geom_vline(xintercept = c(500, 1000), linetype = "dashed", linewidth = 0.3, color = "grey40") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.6) +
    annotate("text", x = -Inf, y = Inf, label = label_text,
             hjust = -0.05, vjust = 1.15, size = 2.1, color = "black", lineheight = 1.1) +
    scale_color_manual(values = treatment_cols) +
    scale_fill_manual(values = treatment_cols) +
    labs(x = "Distance from hematoma edge (um)", y = paste0("P(", CT, ")"), title = panel_label) +
    theme_natmed()
}

fit_window_treatment <- function(data, CT) {
  form <- (cluster_name == CT) ~ Treatment + (1 | Animal/Section)
  fit  <- glmmTMB(form, data, family = binomial())
  coef(summary(fit))$cond %>% data.frame() %>%
    rownames_to_column("coef") %>% rename(se = Std..Error, pvalue = Pr...z..) %>%
    filter(coef %in% c("TreatmentEV", "TreatmentRA")) %>% mutate(CellType = CT)
}

window_stats <- map_dfr(cutoffs, function(cut) {
  data <- df_full %>% filter(dist < cut)
  map_dfr(celltypes, function(CT) fit_window_treatment(data, CT) %>% mutate(Cutoff = cut))
}) %>% mutate(FDR = p.adjust(pvalue, method = "BH"))

morph_window_full <- window_stats %>%
  mutate(PanelType = "windowed_boxplot (panel E)",
         p_BH_interaction_family = NA) %>%
  rename(p_BH_window_family = FDR)

morph_trend_full$p_BH_window_family <- NA

morphology_table <- bind_rows(
  morph_trend_full %>% select(CellType, PanelType, coef, Estimate, se, pvalue,
                              p_BH_interaction_family, p_BH_window_family, Cutoff),
  morph_window_full %>% select(CellType, PanelType, coef, Estimate, se, pvalue,
                               p_BH_interaction_family, p_BH_window_family, Cutoff)
)

write.csv(morphology_table, "Tables/Table19_Figure4DE_morphology_FULL.csv", row.names = FALSE)

make_box_panel <- function(CT, panel_label) {
  window_data <- map_dfr(cutoffs, function(cut) {
    df_full %>% filter(dist < cut) %>% group_by(Animal, Treatment) %>%
      summarise(prop = mean(cluster_name == CT), .groups = "drop") %>%
      mutate(Window = paste0("<", cut, " um"))
  }) %>% mutate(Window = factor(Window, paste0("<", cutoffs, " um")))
  
  p <- ggplot(window_data, aes(Window, prop, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.7,
                 position = position_dodge(0.75), width = 0.6) +
    geom_point(aes(color = Treatment), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
               size = 0.8, alpha = 0.6) +
    scale_fill_manual(values = treatment_cols) +
    scale_color_manual(values = treatment_cols) +
    labs(x = NULL, y = paste0("Mean proportion ", CT), title = panel_label) +
    theme_natmed()
  
  stats_here <- window_stats %>% filter(CellType == CT)
  y_max <- max(window_data$prop, na.rm = TRUE)
  step  <- (max(window_data$prop) - min(window_data$prop)) * 0.12 + 0.02
  
  for (i in seq_along(cutoffs)) {
    cut <- cutoffs[i]; x0 <- i - 0.25; x_ev <- i; x_ra <- i + 0.25
    ev_row <- stats_here %>% filter(Cutoff == cut, coef == "TreatmentEV")
    ra_row <- stats_here %>% filter(Cutoff == cut, coef == "TreatmentRA")
    if (nrow(ev_row) == 1) p <- p + geom_signif(annotations = sig_label(ev_row$FDR),
                                                y_position = y_max + step, xmin = x0, xmax = x_ev, tip_length = 0.01, textsize = 2.5, vjust = 0.3)
    if (nrow(ra_row) == 1) p <- p + geom_signif(annotations = sig_label(ra_row$FDR),
                                                y_position = y_max + step * 2, xmin = x0, xmax = x_ra, tip_length = 0.01, textsize = 2.5, vjust = 0.3)
  }
  p + coord_cartesian(ylim = c(NA, y_max + step * 3))
}

## ============================================================================
## PART C -- assemble
## ============================================================================

top_row <- p_rpm | p_lat | p_mnss

bottom_half <- wrap_plots(
  make_trend_panel("Ameboid", "d   Ameboid"), make_trend_panel("Hypertrophic", "e   Hypertrophic"),
  make_trend_panel("Ramified", "f   Ramified"), make_trend_panel("Rod-like", "g   Rod-like"),
  make_box_panel("Ameboid", "h   Ameboid"), make_box_panel("Hypertrophic", "i   Hypertrophic"),
  make_box_panel("Ramified", "j   Ramified"), make_box_panel("Rod-like", "k   Rod-like"),
  ncol = 4
)

fig_full <- (top_row / bottom_half) +
  plot_layout(heights = c(1, 2.5), guides = "collect") &
  theme(legend.position = "top")
fig_full
fig_full
ggsave("Figure4_FINAL.pdf", fig_full, width = 15, height = 8, units = "in", device = grDevices::cairo_pdf)
ggsave("Figure4_FINAL.png", fig_full, width = 15, height = 8, units = "in", device = grDevices::cairo_pdf)

cat("Behavioral table:", nrow(behav_table), "rows ->", "Tables/Table18_Figure4_behaviour_FULL.csv\n")
cat("Morphology table:", nrow(morphology_table), "rows ->", "Tables/Table19_Figure4DE_morphology_FULL.csv\n")
