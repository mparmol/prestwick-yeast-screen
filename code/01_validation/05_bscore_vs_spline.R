## ============================================================
## 05_bscore_vs_spline.R — B.3 of the analytical plan
## Comparison of positional-effect correction methods:
##   (0) Without correction (baseline)
##   (1) B-score: median polish 8×10 (method used in the pipeline)
##   (2) Thin-plate spline: gam(auc ~ s(row, col), mgcv)
##   (3) Loess 2D: loess(auc ~ row * col)  [only referencia]
##
## Primary metric: SD residual controls after correction
## (controls are NOT corrected → measures how much positional effect
##  remains in the compound space after each method)
##
## Metrics secondary:
##   · ICC between biological replicates (irr::icc)
##   · AUROC recovery core hits (pROC)
##   · SD AUC controls the map 8×12 (diagnostic visual)
##
## Outputs:
##   Table_spatial_correction.txt
##   Fig_spatial_correction.pdf      — Fig new for the paper
##   Fig_spatial_heatmaps.pdf        — heatmaps of positional effect by method
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
library(reshape2)
library(cowplot)

for (pkg in c("mgcv", "pROC", "irr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, quiet = TRUE, repos = "https://cloud.r-project.org")
}

OUT_DIR  <- "./data/intermediate/validation"
DEFS_DIR <- file.path(OUT_DIR, "definitive")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DEFS_DIR, showWarnings = FALSE, recursive = TRUE)

PLATES_EXCLUDE <- c(1, 18)
TIME_MAX_H     <- 60
FDR_THRESH     <- 0.05

## ── Core hits (gold standard for AUROC) ────────────────────────────────────
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

## ── Functions positional correction ───────────────────────────────────────

## Mapeo well_pos (1-96) → (row_idx, col_idx)
well_to_rowcol <- function(wp) {
  wp <- as.integer(wp)
  data.frame(row_idx = ceiling(wp / 12L),
             col_idx = ((wp - 1L) %% 12L) + 1L)
}

## B-score: median polish on matrix 8×10 compounds (cols 2-11)
apply_bscore <- function(auc_df) {
  rc <- well_to_rowcol(auc_df$well_pos)
  auc_df$row_idx <- rc$row_idx
  auc_df$col_idx <- rc$col_idx
  cmpd <- auc_df[auc_df$col_idx %in% 2:11, ]
  if (nrow(cmpd) < 4) return(auc_df)
  mat <- matrix(NA_real_, 8, 10)
  for (k in seq_len(nrow(cmpd))) {
    r <- cmpd$row_idx[k]; c <- cmpd$col_idx[k] - 1L
    if (r >= 1L && r <= 8L && c >= 1L && c <= 10L) mat[r, c] <- cmpd$auc[k]
  }
  mp <- tryCatch(medpolish(mat, na.rm = TRUE, trace.iter = FALSE),
                 error = function(and) NULL)
  if (is.null(mp)) return(auc_df)
  for (k in seq_len(nrow(cmpd))) {
    r <- cmpd$row_idx[k]; c <- cmpd$col_idx[k] - 1L
    if (r >= 1L && r <= 8L && c >= 1L && c <= 10L) {
      idx <- which(auc_df$well_pos == cmpd$well_pos[k])
      re <- if (!is.na(mp$row[r])) mp$row[r] else 0
      ce <- if (!is.na(mp$col[c])) mp$col[c] else 0
      auc_df$auc[idx] <- auc_df$auc[idx] - re - ce
    }
  }
  auc_df
}

## Thin-plate spline: gam(auc ~ s(row_idx, col_idx)) on TODOS the wells
## (compounds + controls; spline fit with all the placa)
apply_spline <- function(auc_df) {
  rc  <- well_to_rowcol(auc_df$well_pos)
  auc_df$row_idx <- rc$row_idx
  auc_df$col_idx <- rc$col_idx
  ok  <- !is.na(auc_df$auc)
  if (sum(ok) < 10) return(auc_df)
  fit <- tryCatch(
    mgcv::gam(auc ~ s(row_idx, col_idx, k = 9), data = auc_df[ok, ]),
    error = function(and) NULL
  )
  if (is.null(fit)) return(auc_df)
  pred <- predict(fit, newdata = auc_df)
  ## Residuo + media global (mantener level)
  grand_mean <- mean(auc_df$auc[ok], na.rm = TRUE)
  auc_df$auc <- auc_df$auc - pred + grand_mean
  auc_df
}

## ── Loop principal by plate-run ────────────────────────────────────────────
cat("Procthatndo plate-runs...\n")

my_design  <- read.csv("./Placas_test/metadata.txt", sep = "\t", encoding = "latin1")
files_all  <- list.files("./data/raw/")
files      <- files_all[!sapply(files_all,
               function(f) as.integer(sub("plate_(\\d+)_.*", "\\1", f)) %in% PLATES_EXCLUDE)]

methods <- c("none", "bscore", "spline")

results_list <- vector("list", length(files))

for (i in seq_along(files)) {
  fname    <- files[i]
  plate_id <- as.integer(sub("plate_(\\d+)_.*", "\\1", fname))
  run_id   <- sub("\\.csv$", "", fname)

  wide <- tryCatch(
    read_wides(files = paste0("./data/raw/", fname)),
    error = function(and) NULL
  )
  if (is.null(wide)) next

  wide <- wide[as.numeric(wide$horas) <= TIME_MAX_H, ]
  t_obs <- as.numeric(wide$horas)

  non_well <- c("file", "horas", "tipo", "Filename")
  well_cols <- setdiff(names(wide), non_well)
  if (length(well_cols) == 0) next

  ## Extract position numeric columns tipo "10.A.86" → 86
  well_pos <- as.integer(sub(".*\\.(\\d+)$", "\\1", well_cols))
  bad <- is.na(well_pos); well_cols <- well_cols[!bad]; well_pos <- well_pos[!bad]

  ## Baseline: rthisr OD t=0
  t0_row <- which.min(t_obs)
  bl <- as.numeric(wide[t0_row, well_cols])

  ## col_idx for identify controls (col 1 and col 12)
  col_idx_all <- ((well_pos - 1L) %% 12L) + 1L
  is_ctrl     <- col_idx_all %in% c(1L, 12L)

  ## AUC by well
  auc_vec <- sapply(seq_along(well_cols), function(j) {
    and <- as.numeric(wide[, well_cols[j]]) - bl[j]
    and[and < 0] <- 0
    gcplyr::auc(x = t_obs, and = and)
  })

  auc_df <- data.frame(
    run      = run_id,
    plate    = plate_id,
    well_pos = well_pos,
    col_idx  = col_idx_all,
    row_idx  = ceiling(well_pos / 12L),
    is_ctrl  = is_ctrl,
    auc      = auc_vec,
    stringsAsFactors = FALSE
  )

  ## Median of controls (used in the normalisation of all methods)
  med_ctrl <- median(auc_df$auc[is_ctrl], na.rm = TRUE)
  if (is.na(med_ctrl) || med_ctrl < 1e-6) next

  ## Apply each method correction and compute metrics
  run_results <- lapply(methods, function(meth) {
    df_corr <- auc_df
    if (meth == "bscore") df_corr <- apply_bscore(df_corr)
    if (meth == "spline") df_corr <- apply_spline(df_corr)

    ## Normalizar by median controls (always the same referencia)
    df_corr$auc_norm <- df_corr$auc / med_ctrl

    ## SD of normalised controls (mide how much noise posicional remains)
    ctrl_norms <- df_corr$auc_norm[is_ctrl]
    sd_ctrl    <- sd(ctrl_norms, na.rm = TRUE)
    cv_ctrl    <- sd_ctrl / mean(ctrl_norms, na.rm = TRUE) * 100

    cmpd_df <- df_corr[!is_ctrl, ]
    if (nrow(cmpd_df) == 0) return(NULL)
    pos_resid <- tryCatch({
      fit_diag <- mgcv::gam(auc_norm ~ s(row_idx, col_idx, k = 9), data = cmpd_df)
      sd(fitted(fit_diag) - mean(fitted(fit_diag)))
    }, error = function(and) NA_real_)

    data.frame(
      run           = run_id,
      plate         = plate_id,
      method        = meth,
      sd_ctrl_norm  = round(sd_ctrl,    4),
      cv_ctrl_pct   = round(cv_ctrl,    2),
      pos_resid_sd  = round(pos_resid,  4),
      n_wells       = nrow(df_corr),
      stringsAsFactors = FALSE
    )
  })

  results_list[[i]] <- do.call(rbind, Filter(Negate(is.null), run_results))

  if (i %% 10 == 0) cat(sprintf("  %d/%d plate-runs...\n", i, length(files)))
}

res_df <- do.call(rbind, Filter(Negate(is.null), results_list))
rownames(res_df) <- NULL

## ── ICC and AUROC per method (requires re-running everything to have all AUCs) ─

res_pipeline <- read.table(
  "./data/intermediate/primary_pipeline/final_results_robust.txt",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

## ICC and AUROC we recompute the normalised AUC cada method
## on the whole library and then compute global metrics
cat("Computing ICC and AUROC per method...\n")

## Tabla run+well_pos → Chemical (from the pipeline results)
chemical_map <- unique(res_pipeline[, c("run", "well_pos", "Chemical")])

global_metrics <- lapply(methods, function(meth) {

  all_norms <- vector("list", length(files))

  for (i in seq_along(files)) {
    fname    <- files[i]
    plate_id <- as.integer(sub("plate_(\\d+)_.*", "\\1", fname))
    run_id   <- sub("\\.csv$", "", fname)

    wide <- tryCatch(
      read_wides(files = paste0("./data/raw/", fname)),
      error = function(and) NULL
    )
    if (is.null(wide)) next
    wide <- wide[as.numeric(wide$horas) <= TIME_MAX_H, ]
    t_obs <- as.numeric(wide$horas)

    non_well  <- c("file", "horas", "tipo", "Filename")
    well_cols <- setdiff(names(wide), non_well)
    if (length(well_cols) == 0) next
    well_pos  <- as.integer(sub(".*\\.(\\d+)$", "\\1", well_cols))
    bad <- is.na(well_pos); well_cols <- well_cols[!bad]; well_pos <- well_pos[!bad]

    t0_row <- which.min(t_obs)
    bl <- as.numeric(wide[t0_row, well_cols])
    col_idx_all <- ((well_pos - 1L) %% 12L) + 1L
    is_ctrl     <- col_idx_all %in% c(1L, 12L)

    auc_vec <- sapply(seq_along(well_cols), function(j) {
      and <- as.numeric(wide[, well_cols[j]]) - bl[j]
      and[and < 0] <- 0
      gcplyr::auc(x = t_obs, and = and)
    })

    auc_df <- data.frame(
      run = run_id, plate = plate_id, well_pos = well_pos,
      col_idx = col_idx_all, row_idx = ceiling(well_pos / 12L),
      is_ctrl = is_ctrl, auc = auc_vec, stringsAsFactors = FALSE
    )
    med_ctrl <- median(auc_df$auc[is_ctrl], na.rm = TRUE)
    if (is.na(med_ctrl) || med_ctrl < 1e-6) next

    if (meth == "bscore") auc_df <- apply_bscore(auc_df)
    if (meth == "spline") auc_df <- apply_spline(auc_df)
    auc_df$auc_norm <- auc_df$auc / med_ctrl

    ## Join compound name from res_pipeline (run + well_pos → Chemical)
    auc_df_merge <- merge(
      auc_df[!is_ctrl, ],
      chemical_map,
      by = c("run", "well_pos"),
      all.x = TRUE
    )
    auc_df_merge <- auc_df_merge[!is.na(auc_df_merge$Chemical), ]
    if (nrow(auc_df_merge) == 0) next

    all_norms[[i]] <- auc_df_merge[, c("run", "plate", "well_pos", "Chemical", "auc_norm")]
  }

  df_all <- do.call(rbind, Filter(Negate(is.null), all_norms))
  if (is.null(df_all) || nrow(df_all) == 0) return(NULL)

  ## ICC: pivot wide by replicate; necesitamos run ~ Chemical
  ## Extract rep identifier from run ("plate_10_rep_1" → rep1)
  df_all$rep <- sub(".*_rep_?", "rep", df_all$run)

  ## Mean across replicates per compound per rep (in case of multiple wells per replicate)
  df_agg <- df_all %>%
    group_by(Chemical, rep) %>%
    summarise(auc_norm = mean(auc_norm, na.rm = TRUE), .groups = "drop")

  ## Pivot for ICC
  df_wide <- df_agg %>%
    pivot_wider(names_from = rep, values_from = auc_norm, values_fn = mean)

  ## Numeric columns only (replicates), drop Chemical
  num_cols <- names(df_wide)[sapply(df_wide, is.numeric)]
  keep_cols <- num_cols[colMeans(is.na(df_wide[, num_cols, drop=FALSE])) < 0.5]
  df_icc <- as.data.frame(df_wide[, keep_cols, drop=FALSE])
  complete_rows <- complete.cases(df_icc)
  icc_val <- tryCatch({
    irr::icc(df_icc[complete_rows, ], model = "twoway",
             type = "agreement", unit = "single")$value
  }, error = function(and) NA_real_)

  ## AUROC: score = mean AUC_norm; label = core_hit
  score_by_cmpd <- df_all %>%
    group_by(Chemical) %>%
    summarise(score = mean(auc_norm, na.rm = TRUE), .groups = "drop")
  score_by_cmpd$is_hit <- score_by_cmpd$Chemical %in% core_hits

  auroc_val <- tryCatch({
    ro <- pROC::roc(score_by_cmpd$is_hit, score_by_cmpd$score,
                    direction = ">", quiet = TRUE)
    as.numeric(pROC::auc(ro))
  }, error = function(and) NA_real_)

  data.frame(
    method = meth,
    icc    = round(icc_val,   4),
    auroc  = round(auroc_val, 4),
    stringsAsFactors = FALSE
  )
})

gm_df <- do.call(rbind, Filter(Negate(is.null), global_metrics))

## ── Tabla summary ─────────────────────────────────────────────────────────
sum_df <- res_df %>%
  group_by(method) %>%
  summarise(
    sd_ctrl_norm_median = round(median(sd_ctrl_norm, na.rm = TRUE), 4),
    cv_ctrl_pct_median  = round(median(cv_ctrl_pct,  na.rm = TRUE), 2),
    pos_resid_sd_median = round(median(pos_resid_sd, na.rm = TRUE), 4),
    n_runs              = n(),
    .groups = "drop"
  )
sum_df <- merge(sum_df, gm_df, by = "method", all.x = TRUE)
sum_df$label <- c(none   = "(0) No correction",
                  bscore = "(1) B-score (median polish)",
                  spline = "(2) Thin-plate spline")[sum_df$method]

cat("\n=== Comparison method positional correction ===\n")
print(sum_df[, c("label","sd_ctrl_norm_median","cv_ctrl_pct_median",
                  "pos_resid_sd_median","icc","auroc")], row.names = FALSE)

write.table(res_df,  file.path(OUT_DIR, "Table_spatial_correction_byrun.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sum_df, file.path(OUT_DIR, "Table_spatial_correction.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Tablas saved.\n")

METHOD_LABELS <- c(
  none   = "(0) No correction",
  bscore = "(1) B-score",
  spline = "(2) Thin-plate spline"
)
res_df$label <- factor(METHOD_LABELS[res_df$method],
                        levels = METHOD_LABELS[c("none","bscore","spline")])

PAL3 <- c("(0) No correction"       = "#AAAAAA",
          "(1) B-score"              = "#2CA02C",
          "(2) Thin-plate spline"    = "#1F77B4")

## Panel A: Normalised control SD per method (boxplot by plate-run)
p_sd <- ggplot(res_df, aes(x = label, and = sd_ctrl_norm, fill = label)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.7, outlier.alpha = 0.4) +
  scale_fill_manual(values = PAL3, guide = "none") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "Within-plate control SD after normalisation",
       subtitle = "Lower = better residual positional correction",
       x = NULL, and = "SD of normalised control AUC")

## Panel B: residual posicional compounds
p_resid <- ggplot(res_df[!is.na(res_df$pos_resid_sd), ],
                   aes(x = label, and = pos_resid_sd, fill = label)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.7, outlier.alpha = 0.4) +
  scale_fill_manual(values = PAL3, guide = "none") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "Residual positional effect in compound AUC",
       subtitle = "SD of fitted spatial trend after correction",
       x = NULL, and = "SD of spatial trend (normalised AUC)")

## Panel C: ICC
p_icc <- ggplot(sum_df, aes(x = label, and = icc, fill = label)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%.3f", icc)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = PAL3, guide = "none") +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "Inter-replicate ICC",
       subtitle = "Higher = better reproducibility",
       x = NULL, and = "ICC (two-way mixed, agreement)")

## Panel D: AUROC
p_auroc <- ggplot(sum_df, aes(x = label, and = auroc, fill = label)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%.4f", auroc)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = PAL3, guide = "none") +
  coord_cartesian(ylim = c(0.9, 1.01)) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "AUROC for core-hit recovery",
       subtitle = "18 core hits vs all other compounds",
       x = NULL, and = "AUROC")

## -- Draft figure (review layout, original dimensions) ----------------
pdf(file.path(OUT_DIR, "Fig_spatial_correction_draft.pdf"), width = 12, height = 10)
grid.arrange(p_sd, p_resid, p_icc, p_auroc, ncol = 2)
dev.off()

## -- Publication figure (PLOS format + TIFF definitive) --------------------
spatial_pub <- cowplot::plot_grid(
  p_sd, p_resid, p_icc, p_auroc,
  ncol = 2,
  labels = c("A","B","C","D"),
  label_size = 18, label_fontface = "bold"
)
pdf(file.path(OUT_DIR, "Fig_spatial_correction.pdf"), width = 6.83, height = 7)
print(spatial_pub)
dev.off()
cat("Fig_spatial_correction.pdf\n")

## ── Heatmaps effect posicional: plate representativa ──────────────────
## Use the plate with the highest control SD as the most illustrative example
cat("Generating heatmaps positional correction...\n")

worst_run <- res_df %>%
  filter(method == "none") %>%
  arrange(desc(pos_resid_sd)) %>%
  slice(1) %>%
  pull(run)

cat(sprintf("Plate-run for heatmap: %s\n", worst_run))

## Recompute AUC for that plate-run with the three methods
fname_w <- paste0(worst_run, ".csv")
wide_w  <- read_wides(files = paste0("./data/raw/", fname_w))
wide_w  <- wide_w[as.numeric(wide_w$horas) <= TIME_MAX_H, ]
t_obs_w <- as.numeric(wide_w$horas)
non_well  <- c("file", "horas", "tipo", "Filename")
wc_w    <- setdiff(names(wide_w), non_well)
wp_w    <- as.integer(sub(".*\\.(\\d+)$", "\\1", wc_w))
bad_w   <- is.na(wp_w); wc_w <- wc_w[!bad_w]; wp_w <- wp_w[!bad_w]
t0_w    <- which.min(t_obs_w)
bl_w    <- as.numeric(wide_w[t0_w, wc_w])

auc_w <- sapply(seq_along(wc_w), function(j) {
  and <- as.numeric(wide_w[, wc_w[j]]) - bl_w[j]
  and[and < 0] <- 0
  gcplyr::auc(x = t_obs_w, and = and)
})
ci_w  <- ((wp_w - 1L) %% 12L) + 1L
ri_w  <- ceiling(wp_w / 12L)
is_c_w <- ci_w %in% c(1L, 12L)
med_w  <- median(auc_w[is_c_w], na.rm = TRUE)

hm_list <- lapply(methods, function(meth) {
  adf <- data.frame(run="x", plate=0, well_pos=wp_w, col_idx=ci_w,
                    row_idx=ri_w, is_ctrl=is_c_w, auc=auc_w)
  if (meth == "bscore") adf <- apply_bscore(adf)
  if (meth == "spline") adf <- apply_spline(adf)
  adf$auc_norm <- adf$auc / med_w
  adf$method_label <- METHOD_LABELS[meth]
  adf
})
hm_df <- do.call(rbind, hm_list)
hm_df$method_label <- factor(hm_df$method_label,
                               levels = METHOD_LABELS[c("none","bscore","spline")])

## Colour scale: deviation from the median value
hm_df$dev <- hm_df$auc_norm - 1

p_hm <- ggplot(hm_df, aes(x = col_idx, and = rev(row_idx), fill = dev)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(data = hm_df[hm_df$is_ctrl, ],
            aes(label = "C"), size = 2.5, color = "black") +
  scale_fill_gradient2(low = "#D73027", mid = "white", high = "#1A9850",
                       midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish,
                       name = "AUC_norm − 1") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(breaks = 1:8,
                     labels = rev(c("A","B","C","D","E","F","G","H"))) +
  facet_wrap(~ method_label, ncol = 3) +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(size = 11, face = "bold"),
        panel.grid = element_blank()) +
  labs(title = sprintf("Spatial heatmap of normalised AUC — plate-run: %s", worst_run),
       subtitle = "Green = higher than control median; red = lower; C = control wells",
       x = "Column", and = "Row")

pdf(file.path(OUT_DIR, "Fig_spatial_heatmaps_draft.pdf"), width = 14, height = 5)
print(p_hm)
dev.off()
pdf(file.path(OUT_DIR, "Fig_spatial_heatmaps.pdf"), width = 6.83, height = 3.5)
print(p_hm)
dev.off()
file.copy(file.path(OUT_DIR, "Table_spatial_correction.txt"),
          file.path(DEFS_DIR, "S6_Table_SpatialCorrection.txt"), overwrite = TRUE)
cat("Figura saved: Fig_spatial_heatmaps.pdf\n")

## ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=== Summary B.3 — B-score vs Thin-plate spline ===\n")
for (k in seq_len(nrow(sum_df))) {
  cat(sprintf("  %-30s  SD_ctrl=%.4f  CV%%=%.1f  posResid=%.4f  ICC=%.3f  AUROC=%.4f\n",
              sum_df$label[k], sum_df$sd_ctrl_norm_median[k],
              sum_df$cv_ctrl_pct_median[k], sum_df$pos_resid_sd_median[k],
              sum_df$icc[k], sum_df$auroc[k]))
}
cat("\n=== 05_bscore_vs_spline.R completedo ===\n")
