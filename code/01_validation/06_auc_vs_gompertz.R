## ============================================================
## 06_auc_vs_gompertz.R — B.4 of the analytical plan
## Comparison of empirical AUC (pipeline statistic) vs
## logistic-model parameters (growthcurver::SummarizeGrowth)
##
## Questions:
##   (1) What fraction of OD600 curves fit well to the model
##       logistic (R² ≥ 0.95 / 0.90)?
##   (2) Are the compounds that fit poorly hits or non-hits?
##   (3) How does k (normalised carrying capacity) correlate
##       with AUC_rel of the pipeline?
##   (4) Does the carrying capacity k detects the same 18 core hits
##       that AUC? (compared by AUROC)
##
## Note: growthcurver uses raw OD600 data (without baseline-correction)
##       because the logistic model requires N(0) = N0 > 0.
##       The pipeline AUC (B-score + median) is loaded directly
##       from final_results_robust.txt for the final comparison.
##
## Outputs:
##   Table_growthcurver_bywell.txt   — fit by well × plate-run
##   Table_growthcurver_summary.txt  — global metrics by class
##   Fig_growthcurver_rsq.pdf        — R² distribution by class
##   Fig_growthcurver_k_vs_auc.pdf   — scatter k_norm vs AUC_rel
##   Fig_growthcurver_bad_fits.pdf   — profiles of curves with worst fit
## ============================================================

## ── Anchor working directory ──────────────────────────────────────────────────
.find_root <- function() {
  d <- normalizePath(getwd())
  for (i in seq_len(6)) {
    if (dir.exists(file.path(d, "data/raw"))) return(d)
    d <- dirname(d)
  }
  stop("No se found the repository root.")
}
setwd(.find_root())
cat(sprintf("Working directory: %s\n", getwd()))

library(gcplyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)

for (pkg in c("growthcurver", "pROC")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, quiet = TRUE, repos = "https://cloud.r-project.org")
}
library(growthcurver)

OUT_DIR  <- "./data/intermediate/validation"
DEFS_DIR <- file.path(OUT_DIR, "definitive")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DEFS_DIR, showWarnings = FALSE, recursive = TRUE)

PLATES_EXCLUDE <- c(1, 18)
TIME_MAX_H     <- 60
RSQ_GOOD       <- 0.95

## ── Core hits and AUC of the pipeline ─────────────────────────────────────────────
res_lmm   <- read.table("./data/intermediate/primary_pipeline/statistics_LMM.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
res_anova <- read.table("./data/intermediate/primary_pipeline/statistics_robust.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
core_hits <- setdiff(
  intersect(res_lmm[res_lmm$FDR_BH  < 0.05, "Chemical"],
            res_anova[res_anova$FDR_BH < 0.05, "Chemical"]),
  "Chicago sky blue 6B"
)
cat(sprintf("Core hits: %d\n", length(core_hits)))

## AUC of the pipeline (B-score + median, already corregido)
res_pipeline <- read.table(
  "./data/intermediate/primary_pipeline/final_results_robust.txt",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
## Mapa run + well_pos → Chemical (for controls and compounds)
chemical_map <- unique(res_pipeline[, c("run", "well_pos", "Chemical")])

## AUC_rel medio per compound (referencia of the pipeline)
auc_rel_mean <- res_pipeline %>%
  group_by(Chemical) %>%
  summarise(AUC_rel_pipeline = mean(AUC_relative, na.rm = TRUE), .groups = "drop")

## ── Loop principal: fit growthcurver by well × plate-run ────────────────
files_all <- list.files("./data/raw/")
files     <- files_all[!sapply(files_all,
               function(f) as.integer(sub("plate_(\\d+)_.*", "\\1", f)) %in% PLATES_EXCLUDE)]
cat(sprintf("Plate-runs procthatr: %d\n", length(files)))

fit_list <- vector("list", length(files))

for (i in seq_along(files)) {
  fname    <- files[i]
  plate_id <- as.integer(sub("plate_(\\d+)_.*", "\\1", fname))
  run_id   <- sub("\\.csv$", "", fname)

  wide <- tryCatch(
    read_wides(files = paste0("./data/raw/", fname)),
    error = function(and) NULL
  )
  if (is.null(wide)) next

  wide     <- wide[as.numeric(wide$horas) <= TIME_MAX_H, ]
  t_obs    <- as.numeric(wide$horas)
  non_well <- c("file", "horas", "tipo", "Filename")
  well_cols <- setdiff(names(wide), non_well)
  if (length(well_cols) == 0) next

  well_pos    <- as.integer(sub(".*\\.(\\d+)$", "\\1", well_cols))
  bad         <- is.na(well_pos)
  well_cols   <- well_cols[!bad]
  well_pos    <- well_pos[!bad]
  col_idx_all <- ((well_pos - 1L) %% 12L) + 1L
  is_ctrl     <- col_idx_all %in% c(1L, 12L)

  ## AUC empirical baseline-corrected (identical to the pipeline, without B-score)
  t0_row <- which.min(t_obs)
  bl <- as.numeric(wide[t0_row, well_cols])
  auc_raw <- sapply(seq_along(well_cols), function(j) {
    and <- as.numeric(wide[, well_cols[j]]) - bl[j]
    and[and < 0] <- 0
    gcplyr::auc(x = t_obs, and = and)
  })
  med_ctrl_auc <- median(auc_raw[is_ctrl], na.rm = TRUE)
  if (is.na(med_ctrl_auc) || med_ctrl_auc < 1e-6) next

  ## growthcurver (OD raw, without baseline correction)
  gc_fits <- lapply(seq_along(well_cols), function(j) {
    od_raw <- as.numeric(wide[, well_cols[j]])
    if (all(is.na(od_raw))) return(NULL)
    ## growthcurver requires N > 0; reemplazar NA and negative by minimum observado
    od_raw[is.na(od_raw)] <- min(od_raw[!is.na(od_raw)], na.rm = TRUE)
    od_raw[od_raw <= 0]   <- min(od_raw[od_raw > 0], na.rm = TRUE)
    tryCatch({
      fit <- growthcurver::SummarizeGrowth(t_obs, od_raw)
      ## R² = 1 - SS_res/SS_tot computed from the nls model
      r2 <- tryCatch({
        1 - sum(residuals(fit$model)^2) / sum((od_raw - mean(od_raw))^2)
      }, error = function(and) NA_real_)
      data.frame(
        r.squared = r2,
        k         = fit$vals$k,
        r_fit     = fit$vals$r,
        t_mid     = fit$vals$t_mid,
        sigma     = fit$vals$sigma,
        auc_l     = fit$vals$auc_l,
        note      = if (!is.null(fit$vals$note)) fit$vals$note else ""
      )
    }, error = function(and) {
      data.frame(r.squared=NA_real_, k=NA_real_, r_fit=NA_real_,
                 t_mid=NA_real_, sigma=NA_real_, auc_l=NA_real_, note="error")
    })
  })

  gc_df <- do.call(rbind, lapply(gc_fits, function(x)
    if (is.null(x)) data.frame(r.squared=NA_real_, k=NA_real_, r_fit=NA_real_,
                                t_mid=NA_real_, sigma=NA_real_, auc_l=NA_real_, note="NULL")
    else x))
  rownames(gc_df) <- NULL

  ## Mediank controls
  k_ctrl  <- gc_df$k[is_ctrl]
  med_k   <- median(k_ctrl, na.rm = TRUE)

  df_plate <- data.frame(
    run       = run_id,
    plate     = plate_id,
    well_pos  = well_pos,
    col_idx   = col_idx_all,
    is_ctrl   = is_ctrl,
    auc_raw   = auc_raw,
    auc_norm  = auc_raw / med_ctrl_auc,
    r.squared = gc_df$r.squared,
    k         = gc_df$k,
    r_fit     = gc_df$r_fit,
    t_mid     = gc_df$t_mid,
    sigma     = gc_df$sigma,
    k_norm    = if (!is.na(med_k) && med_k > 1e-6) gc_df$k / med_k else NA_real_,
    auc_l     = gc_df$auc_l,
    note      = gc_df$note,
    stringsAsFactors = FALSE
  )
  fit_list[[i]] <- df_plate

  if (i %% 10 == 0) cat(sprintf("  %d/%d plate-runs...\n", i, length(files)))
}

fit_df <- do.call(rbind, Filter(Negate(is.null), fit_list))
rownames(fit_df) <- NULL
cat(sprintf("Total wells analizados: %d\n", nrow(fit_df)))

## ── Join Chemical from map of the pipeline ────────────────────────────────────
fit_df <- merge(fit_df, chemical_map,
                by = c("run", "well_pos"), all.x = TRUE)
fit_df$Chemical[is.na(fit_df$Chemical) & fit_df$is_ctrl] <- "CONTROL"

## Clase: CONTROL / HIT / COMPOUND
fit_df$class <- "COMPOUND"
fit_df$class[fit_df$is_ctrl] <- "CONTROL"
fit_df$class[!is.na(fit_df$Chemical) & fit_df$Chemical %in% core_hits] <- "HIT"

## ── Statistics fit by class ─────────────────────────────────────────
cat("\n=== Calidad fit logistic (R²) ===\n")

rsq_stats <- fit_df %>%
  filter(!is.na(r.squared)) %>%
  group_by(class) %>%
  summarise(
    n             = n(),
    rsq_median    = round(median(r.squared, na.rm=TRUE), 4),
    rsq_mean      = round(mean(r.squared,   na.rm=TRUE), 4),
    pct_ge_095    = round(mean(r.squared >= 0.95, na.rm=TRUE) * 100, 1),
    pct_ge_090    = round(mean(r.squared >= 0.90, na.rm=TRUE) * 100, 1),
    pct_lt_080    = round(mean(r.squared <  0.80, na.rm=TRUE) * 100, 1),
    .groups = "drop"
  )
print(rsq_stats, row.names=FALSE)

## Global
rsq_global <- fit_df %>%
  filter(!is.na(r.squared)) %>%
  summarise(
    n          = n(),
    rsq_median = round(median(r.squared, na.rm=TRUE), 4),
    pct_ge_095 = round(mean(r.squared >= 0.95, na.rm=TRUE) * 100, 1),
    pct_ge_090 = round(mean(r.squared >= 0.90, na.rm=TRUE) * 100, 1),
    pct_lt_080 = round(mean(r.squared <  0.80, na.rm=TRUE) * 100, 1)
  )
cat("\nGlobal (all the wells):\n")
print(rsq_global, row.names=FALSE)

## ── Correlacion k_norm vs AUC_norm (non-hit compounds) ─────────────────────
## Excluir wells without Chemical asignado
cmpd_df <- fit_df %>%
  filter(!is_ctrl, !is.na(Chemical), !is.na(k_norm), !is.na(auc_norm)) %>%
  group_by(Chemical) %>%
  summarise(
    k_norm_mean   = mean(k_norm,   na.rm=TRUE),
    auc_norm_mean = mean(auc_norm, na.rm=TRUE),
    rsq_mean      = mean(r.squared, na.rm=TRUE),
    is_hit        = first(class) == "HIT",
    .groups = "drop"
  ) %>%
  left_join(auc_rel_mean, by="Chemical")

ok <- !is.na(cmpd_df$k_norm_mean) & !is.na(cmpd_df$auc_norm_mean)
r_pearson  <- cor(cmpd_df$k_norm_mean[ok], cmpd_df$auc_norm_mean[ok])
r_spearman <- cor(cmpd_df$k_norm_mean[ok], cmpd_df$auc_norm_mean[ok], method="spearman")
cat(sprintf("\nCorrelacion k_norm vs AUC_norm (per compound): Pearson r=%.3f, Spearman ρ=%.3f\n",
            r_pearson, r_spearman))

## with AUC_rel of the pipeline (if available)
ok2 <- ok & !is.na(cmpd_df$AUC_rel_pipeline)
if (sum(ok2) > 0) {
  r_pipeline <- cor(cmpd_df$k_norm_mean[ok2], cmpd_df$AUC_rel_pipeline[ok2])
  cat(sprintf("Correlacion k_norm vs AUC_rel (pipeline B-score): Pearson r=%.3f\n", r_pipeline))
}

## ── AUROC: k_norm vs AUC_norm as score for detect hits ──────────────────
auroc_k <- tryCatch({
  ro <- pROC::roc(cmpd_df$is_hit[ok], cmpd_df$k_norm_mean[ok],
                  direction = ">", quiet = TRUE)
  as.numeric(pROC::auc(ro))
}, error = function(and) NA_real_)

auroc_auc <- tryCatch({
  ro <- pROC::roc(cmpd_df$is_hit[ok], cmpd_df$auc_norm_mean[ok],
                  direction = ">", quiet = TRUE)
  as.numeric(pROC::auc(ro))
}, error = function(and) NA_real_)

cat(sprintf("AUROC k_norm   (detection hits): %.4f\n", auroc_k))
cat(sprintf("AUROC AUC_norm (detection hits): %.4f\n", auroc_auc))

## ── Tabla global metrics ──────────────────────────────────────────────────
metrics_df <- data.frame(
  statistic = c("N wells total", "N wells fit attempted", "N wells fit failed (NA r²)",
                "R² median (all wells)", "R² median (compounds)", "R² median (controls)",
                "R² median (hits)", "% wells R²≥0.95", "% wells R²≥0.90", "% wells R²<0.80",
                "Corr k_norm vs AUC_norm (Pearson)", "Corr k_norm vs AUC_norm (Spearman)",
                "AUROC k_norm (hits)", "AUROC AUC_norm (hits)"),
  value = c(nrow(fit_df),
            sum(!is.na(fit_df$r.squared)),
            sum(is.na(fit_df$r.squared)),
            round(median(fit_df$r.squared, na.rm=TRUE), 4),
            round(median(fit_df$r.squared[!fit_df$is_ctrl], na.rm=TRUE), 4),
            round(median(fit_df$r.squared[fit_df$is_ctrl],  na.rm=TRUE), 4),
            round(median(fit_df$r.squared[fit_df$class=="HIT"], na.rm=TRUE), 4),
            round(mean(fit_df$r.squared >= 0.95, na.rm=TRUE)*100, 1),
            round(mean(fit_df$r.squared >= 0.90, na.rm=TRUE)*100, 1),
            round(mean(fit_df$r.squared <  0.80, na.rm=TRUE)*100, 1),
            round(r_pearson, 4),
            round(r_spearman, 4),
            round(auroc_k, 4),
            round(auroc_auc, 4)),
  stringsAsFactors = FALSE
)
cat("\n=== Tabla metrics B.4 ===\n")
print(metrics_df, row.names=FALSE)

write.table(fit_df[, c("run","plate","well_pos","class","Chemical",
                        "auc_raw","auc_norm","r.squared","k","k_norm",
                        "r_fit","t_mid","sigma","auc_l","note")],
            file.path(OUT_DIR, "Table_growthcurver_bywell.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(metrics_df, file.path(OUT_DIR, "Table_growthcurver_summary.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)
cat("Tablas saved.\n")

## ── Figuras ───────────────────────────────────────────────────────────────────
CLASE_COLS <- c(CONTROL="#2CA02C", COMPOUND="#AAAAAA", HIT="#D62728")
CLASE_LABS <- c(CONTROL="Controls (n=16/plate)",
                COMPOUND="Compounds (non-hits)",
                HIT=sprintf("Core hits (n=%d)", length(core_hits)))

## ── Fig 1: Distribution R² by class ──────────────────────────────────────
rsq_plot_df <- fit_df %>%
  filter(!is.na(r.squared)) %>%
  mutate(class = factor(class, levels=c("CONTROL","COMPOUND","HIT")))

p_rsq_dens <- ggplot(rsq_plot_df, aes(x=r.squared, fill=class, color=class)) +
  geom_density(alpha=0.35, linewidth=0.7) +
  geom_vline(xintercept=0.95, linetype="dashed", color="black", linewidth=0.6) +
  geom_vline(xintercept=0.90, linetype="dotted", color="black", linewidth=0.5) +
  scale_fill_manual(values=CLASE_COLS, labels=CLASE_LABS, name=NULL) +
  scale_color_manual(values=CLASE_COLS, labels=CLASE_LABS, name=NULL) +
  theme_bw(base_size=13) +
  labs(title="Distribution of logistic-fit R² by well class",
       subtitle="Dashed: R²=0.95  Dotted: R²=0.90",
       x="Goodness of fit (R²)", and="Density")

p_rsq_box <- ggplot(rsq_plot_df, aes(x=class, and=r.squared, fill=class)) +
  geom_boxplot(notch=TRUE, outlier.size=0.5, outlier.alpha=0.3) +
  geom_hline(yintercept=0.95, linetype="dashed", color="black", linewidth=0.6) +
  scale_fill_manual(values=CLASE_COLS, guide="none") +
  scale_x_discrete(labels=CLASE_LABS) +
  theme_bw(base_size=13) +
  theme(axis.text.x=element_text(angle=20, hjust=1)) +
  labs(title="R² boxplot by well class",
       x=NULL, and="Goodness of fit (R²)")

## Barchart: % below thresholds by class
thresh_df <- rsq_plot_df %>%
  group_by(class) %>%
  summarise(
    pct_ge095 = mean(r.squared >= 0.95)*100,
    pct_090_095 = mean(r.squared >= 0.90 & r.squared < 0.95)*100,
    pct_lt090   = mean(r.squared <  0.90)*100,
    .groups="drop"
  ) %>%
  pivot_longer(cols=starts_with("pct"), names_to="range", values_to="pct") %>%
  mutate(range=factor(range,
           levels=c("pct_lt090","pct_090_095","pct_ge095"),
           labels=c("R²<0.90","0.90≤R²<0.95","R²≥0.95")))

p_thresh <- ggplot(thresh_df, aes(x=class, and=pct, fill=range)) +
  geom_col(width=0.6) +
  scale_fill_manual(values=c("R²<0.90"="#D62728","0.90≤R²<0.95"="#FF7F0E","R²≥0.95"="#2CA02C"),
                    name="Fit quality") +
  scale_x_discrete(labels=CLASE_LABS) +
  theme_bw(base_size=13) +
  theme(axis.text.x=element_text(angle=20, hjust=1)) +
  labs(title="Fit quality distribution by class (%)",
       x=NULL, and="% of wells")

## -- Draft figure (review layout, original dimensions) ----------------
pdf(file.path(OUT_DIR, "Fig_growthcurver_rsq_draft.pdf"), width = 14, height = 5)
grid.arrange(p_rsq_dens, p_rsq_box, p_thresh, ncol = 3)
dev.off()

## -- Publication figure (PLOS format + TIFF definitive) --------------------
rsq_pub <- cowplot::plot_grid(
  p_rsq_dens, p_rsq_box, p_thresh,
  ncol = 3,
  labels = c("A","B","C"),
  label_size = 18, label_fontface = "bold"
)
pdf(file.path(OUT_DIR, "Fig_growthcurver_rsq.pdf"), width = 6.83, height = 3.5)
print(rsq_pub)
dev.off()
cat("Figura saved: Fig_growthcurver_rsq.pdf\n")

## ── Fig 2: k_norm vs AUC_norm scatter ────────────────────────────────────────
cmpd_plot <- cmpd_df %>%
  filter(!is.na(k_norm_mean), !is.na(auc_norm_mean)) %>%
  mutate(clase2 = ifelse(is_hit, "HIT", "COMPOUND"))

p_scatter <- ggplot(cmpd_plot, aes(x=auc_norm_mean, and=k_norm_mean, color=clase2)) +
  geom_point(size=1.2, alpha=0.6) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  scale_color_manual(values=c(COMPOUND="#AAAAAA", HIT="#D62728"),
                     labels=c(COMPOUND="Non-hit compounds", HIT=sprintf("Core hits (n=%d)",length(core_hits))),
                     name=NULL) +
  annotate("text", x=Inf, and=-Inf, hjust=1.05, vjust=-0.3, size=3.5,
           label=sprintf("Pearson r = %.3f\nSpearman rho = %.3f", r_pearson, r_spearman)) +
  theme_bw(base_size=13) +
  labs(title="Logistic k (carrying capacity) vs AUC",
       subtitle="Both normalised by plate-run control median",
       x="AUC_norm (empirical, baseline-corrected)",
       and="k_norm (logistic carrying capacity)")

## k_norm vs AUC_rel of the pipeline (B-score corrected)
if (sum(!is.na(cmpd_plot$AUC_rel_pipeline)) > 0) {
  p_scatter2 <- ggplot(cmpd_plot[!is.na(cmpd_plot$AUC_rel_pipeline),],
                       aes(x=AUC_rel_pipeline, and=k_norm_mean, color=clase2)) +
    geom_point(size=1.2, alpha=0.6) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
    scale_color_manual(values=c(COMPOUND="#AAAAAA", HIT="#D62728"),
                       name=NULL, guide="none") +
    annotate("text", x=Inf, and=-Inf, hjust=1.05, vjust=-0.3, size=3.5,
             label=sprintf("Pearson r = %.3f", if(exists("r_pipeline")) round(r_pipeline,3) else NA)) +
    theme_bw(base_size=13) +
    labs(title="k_norm vs AUC_rel (pipeline, B-score corrected)",
         x="AUC_rel (pipeline, B-score + median norm.)",
         and="k_norm (logistic carrying capacity)")
} else {
  p_scatter2 <- p_scatter
}

## AUROC comparison panel
auroc_comp <- data.frame(
  score  = c("k (carrying capacity)", "AUC (empirical)"),
  auroc  = c(auroc_k, auroc_auc)
)
p_auroc <- ggplot(auroc_comp, aes(x=score, and=auroc, fill=score)) +
  geom_col(width=0.5) +
  geom_text(aes(label=sprintf("%.4f", auroc)), vjust=-0.4, size=4) +
  scale_fill_manual(values=c("k (carrying capacity)"="#1F77B4", "AUC (empirical)"="#2CA02C"),
                    guide="none") +
  coord_cartesian(ylim=c(0.9, 1.01)) +
  theme_bw(base_size=13) +
  labs(title="AUROC for core-hit recovery",
       subtitle="18 core hits vs all other compounds",
       x=NULL, and="AUROC")

## -- Draft figure (review layout, original dimensions) ----------------
pdf(file.path(OUT_DIR, "Fig_growthcurver_k_vs_auc_draft.pdf"), width = 14, height = 5)
grid.arrange(p_scatter, p_scatter2, p_auroc, ncol = 3)
dev.off()

## -- Publication figure (PLOS format + TIFF definitive) --------------------
kauc_pub <- cowplot::plot_grid(
  p_scatter, p_scatter2, p_auroc,
  ncol = 3,
  labels = c("A","B","C"),
  label_size = 18, label_fontface = "bold"
)
pdf(file.path(OUT_DIR, "Fig_growthcurver_k_vs_auc.pdf"), width = 6.83, height = 3.5)
print(kauc_pub)
dev.off()
cat("Figura saved: Fig_growthcurver_k_vs_auc.pdf\n")

## ── Fig 3: profiles the 9 curves with worst fit (compounds) ─────────────
cat("Generating profiles of curves with worst fit...\n")

worst_wells <- fit_df %>%
  filter(!is_ctrl, !is.na(r.squared), !is.na(Chemical)) %>%
  arrange(r.squared) %>%
  slice_head(n = 9) %>%
  select(run, well_pos, class, Chemical, r.squared, k_norm, auc_norm)

if (nrow(worst_wells) > 0) {
  bad_curves <- vector("list", nrow(worst_wells))
  for (wi in seq_len(nrow(worst_wells))) {
    fname_b <- paste0(worst_wells$run[wi], ".csv")
    wide_b <- tryCatch(read_wides(files=paste0("./data/raw/",fname_b)), error=function(and) NULL)
    if (is.null(wide_b)) next
    wide_b <- wide_b[as.numeric(wide_b$horas) <= TIME_MAX_H, ]
    t_b    <- as.numeric(wide_b$horas)
    non_well <- c("file","horas","tipo","Filename")
    wc_b <- setdiff(names(wide_b), non_well)
    wp_b <- as.integer(sub(".*\\.(\\d+)$","\\1",wc_b))
    bad_b <- is.na(wp_b); wc_b <- wc_b[!bad_b]; wp_b <- wp_b[!bad_b]
    j <- which(wp_b == worst_wells$well_pos[wi])
    if (length(j) == 0) next
    od_raw <- as.numeric(wide_b[, wc_b[j[1]]])
    od_raw_clean <- od_raw
    od_raw_clean[is.na(od_raw_clean)] <- min(od_raw_clean, na.rm=TRUE)
    od_raw_clean[od_raw_clean <= 0] <- min(od_raw_clean[od_raw_clean > 0], na.rm=TRUE)
    ## Logistic fitted
    fit_b <- tryCatch(growthcurver::SummarizeGrowth(t_b, od_raw_clean), error=function(and) NULL)
    chem_short <- substr(worst_wells$Chemical[wi], 1, 22)
    label_b <- sprintf("%s\nR²=%.3f  k=%.2f", chem_short,
                       worst_wells$r.squared[wi], worst_wells$k_norm[wi])
    df_b <- data.frame(t=t_b, od=od_raw, type="observed", label=label_b)
    if (!is.null(fit_b)) {
      t_seq <- seq(0, max(t_b), length.out=200)
      od_fit <- tryCatch(predict(fit_b$model, newdata=list(t=t_seq)), error=function(and) NULL)
      if (!is.null(od_fit)) {
        df_b <- rbind(df_b, data.frame(t=t_seq, od=od_fit, type="logistic", label=label_b))
      }
    }
    bad_curves[[wi]] <- df_b
  }

  bad_df <- do.call(rbind, Filter(Negate(is.null), bad_curves))
  if (!is.null(bad_df) && nrow(bad_df) > 0) {
    p_bad <- ggplot(bad_df, aes(x=t, and=od, color=type, linetype=type)) +
      geom_line(linewidth=0.7) +
      scale_color_manual(values=c(observed="#333333", logistic="#D62728"), name=NULL) +
      scale_linetype_manual(values=c(observed="solid", logistic="dashed"), name=NULL) +
      facet_wrap(~ label, ncol=3, scales="free_y") +
      theme_bw(base_size=11) +
      theme(strip.text=element_text(size=9), legend.position="top") +
      labs(title="Growth curves with poorest logistic fit (worst 9 compounds)",
           x="Time (h)", and="OD600 (raw)")
    pdf(file.path(OUT_DIR, "Fig_growthcurver_bad_fits_draft.pdf"), width = 12, height = 9)
    print(p_bad)
    dev.off()
    pdf(file.path(OUT_DIR, "Fig_growthcurver_bad_fits.pdf"), width = 6.83, height = 7)
    print(p_bad)
    dev.off()
    cat("Figura saved: Fig_growthcurver_bad_fits.pdf\n")
  }
}

## ── Summary ───────────────────────────────────────────────────────────────────
file.copy(file.path(OUT_DIR, "Table_growthcurver_summary.txt"),
          file.path(DEFS_DIR, "S7_Table_Growthcurver.txt"), overwrite = TRUE)

cat("\n=== Summary B.4 — AUC vs logistic model (growthcurver) ===\n")
cat(sprintf("  R² median global:       %.4f\n", median(fit_df$r.squared, na.rm=TRUE)))
cat(sprintf("  %% wells R² ≥ 0.95:      %.1f%%\n", mean(fit_df$r.squared >= 0.95, na.rm=TRUE)*100))
cat(sprintf("  %% wells R² ≥ 0.90:      %.1f%%\n", mean(fit_df$r.squared >= 0.90, na.rm=TRUE)*100))
cat(sprintf("  Correlacion k_norm ~ AUC_norm: r=%.3f (Pearson), ρ=%.3f (Spearman)\n",
            r_pearson, r_spearman))
cat(sprintf("  AUROC k_norm (hits):     %.4f\n", auroc_k))
cat(sprintf("  AUROC AUC_norm (hits):   %.4f\n", auroc_auc))
cat("\n=== 06_auc_vs_gompertz.R completedo ===\n")
