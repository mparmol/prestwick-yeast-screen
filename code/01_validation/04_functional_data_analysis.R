## ============================================================
## 04_functional_data_analysis.R — A.5 of the analytical plan
## Functional Data Analysis (FDA) on curves OD600
##
## Analysis:
##   1. fPCA (fpca.face, refund) on all the curves compound
##      → components interpretables (lag, tasa, plateau)
##   2. ANOVA on scores fPC per compound vs controls
##      → lista hits functionales; comparison with hits of the AUC-ANOVA
##   3. Proyeccion core hits space fPC
##      → characterisation of the mode of action temporal the hits
##
## Outputs:
##   Table_fPCA_summary.txt       — explained variance by PC
##   Table_fPCA_hits.txt          — compounds with effect functional significant
##   Fig_fPCA_components.pdf      — media + effects the 4 primeras PCs
##   Fig_fPCA_scores.pdf          — scatter PC1 vs PC2 (hits highlighted)
##   Fig_fPCA_hit_profiles.pdf    — curves individuales core hits
## ============================================================

## ── Anchor working directory ──────────────────────────────────────────────────
.find_root <- function() {
  d <- normalizePath(getwd())
  for (i in seq_len(6)) {
    if (dir.exists(file.path(d, "data/raw"))) return(d)
    d <- dirname(d)
  }
  stop("Could not find the project root (data/raw not found).")
}
setwd(.find_root())
cat(sprintf("Working directory: %s\n", getwd()))

library(gcplyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)

## refund for fPCA
for (pkg in c("refund", "reshape2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, quiet = TRUE, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("No could instalar ", pkg, ". Install it manually.")
}

OUT_DIR  <- "./data/intermediate/validation"
DEFS_DIR <- file.path(OUT_DIR, "definitive")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DEFS_DIR, showWarnings = FALSE, recursive = TRUE)

PLATES_EXCLUDE <- c(1, 18)
TIME_MAX_H     <- 60
## Grid temporal common: cada 30 min 0 60 h
TGRID          <- seq(0, TIME_MAX_H, by = 0.5)
N_FPCS         <- 6        # PCs retener
FDR_THRESH     <- 0.05

## ── Functions auxiliares ─────────────────────────────────────────────────────
plate_num_from_file <- function(fname)
  as.integer(sub("plate_(\\d+)_.*", "\\1", fname))

## Extract well_pos numeric of the name column "10.A.86" → 86
well_pos_from_colname <- function(x)
  as.integer(sub(".*\\.(\\d+)$", "\\1", x))

## Interpolation linear a curve the grid temporal common
interp_to_grid <- function(t_obs, y_obs, tgrid) {
  ok <- !is.na(y_obs) & !is.na(t_obs)
  if (sum(ok) < 2) return(rep(NA_real_, length(tgrid)))
  approx(t_obs[ok], y_obs[ok], xout = tgrid, rule = 2)$and
}

## ── Cargar pipeline results for metadata ──────────────────────────────
res_lmm   <- read.table("./data/intermediate/primary_pipeline/statistics_LMM.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
res_anova <- read.table("./data/intermediate/primary_pipeline/statistics_robust.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
EXCLUDE_ARTEFACTS <- "Chicago sky blue 6B"
core_hits <- setdiff(
  intersect(res_lmm[res_lmm$FDR_BH  < 0.05, "Chemical"],
            res_anova[res_anova$FDR_BH < 0.05, "Chemical"]),
  EXCLUDE_ARTEFACTS
)
cat(sprintf("Core hits (ANOVA ∩ LMM, FDR<0.05): %d\n", length(core_hits)))

## ── Cargar design and files plate ────────────────────────────────────────
my_design  <- read.csv(".", sep = "\t", encoding = "latin1")
files_all  <- list.files("./data/raw/")
files      <- files_all[!sapply(files_all,
               function(f) plate_num_from_file(f) %in% PLATES_EXCLUDE)]
cat(sprintf("Plate-runs procthatr: %d\n", length(files)))

## ── Loop: load and normalizar curves ─────────────────────────────────────────
cat("\n")

## Lists accumulation: row = curve (compound or control)
curve_mat_list  <- vector("list", length(files))   # matrices effect [wells × tgrid]
curve_norm_list <- vector("list", length(files))   # matrices normalised (for figures)
meta_list       <- vector("list", length(files))

for (i in seq_along(files)) {
  fname    <- files[i]
  plate_id <- plate_num_from_file(fname)
  run_id   <- sub("\\.csv$", "", fname)

  wide <- tryCatch(
    read_wides(files = paste0("./data/raw/", fname)),
    error = function(and) NULL
  )
  if (is.null(wide)) next

  ## Filter by time maximum and convert numeric
  wide <- wide[as.numeric(wide$horas) <= TIME_MAX_H, ]
  t_obs <- as.numeric(wide$horas)

  ## wells: we exclude "file" and "horas" (gcplyr adds "file")
  non_well_cols <- c("file", "horas", "tipo", "Filename")
  well_cols <- setdiff(names(wide), non_well_cols)
  if (length(well_cols) == 0) next

  ## Extract position numeric: "10.A.86" → 86
  well_pos <- well_pos_from_colname(well_cols)
  bad <- is.na(well_pos)
  well_cols <- well_cols[!bad]
  well_pos  <- well_pos[!bad]

  ## Baseline: rthisr OD t=0 well well (identical to the pipeline)
  t0_row <- which.min(t_obs)
  baseline_vals <- as.numeric(wide[t0_row, well_cols])

  ## col_idx 1 and 12 (first and last column of the layout 8×12)
  col_idx_all <- ((well_pos - 1L) %% 12L) + 1L
  is_ctrl_pos <- col_idx_all %in% c(1L, 12L)
  ctrl_cols   <- well_cols[is_ctrl_pos]

  if (length(ctrl_cols) == 0) next

  ## AUC controls for median normalisation
  ctrl_auc <- sapply(ctrl_cols, function(wc) {
    and <- as.numeric(wide[, wc]) - baseline_vals[match(wc, well_cols)]
    and[and < 0] <- 0
    gcplyr::auc(x = t_obs, and = and)
  })
  med_ctrl <- median(ctrl_auc, na.rm = TRUE)
  if (is.na(med_ctrl) || med_ctrl < 1e-6) next

  n_wells  <- length(well_cols)
  mat_norm <- matrix(NA_real_, nrow = n_wells, ncol = length(TGRID))

  for (j in seq_along(well_cols)) {
    wc  <- well_cols[j]
    and   <- as.numeric(wide[, wc]) - baseline_vals[j]
    and[and < 0] <- 0
    y_interp <- interp_to_grid(t_obs, and, TGRID)
    mat_norm[j, ] <- y_interp / max(med_ctrl, 1e-6)
  }

  ctrl_idx_in_well <- is_ctrl_pos
  ctrl_mean_curve  <- colMeans(mat_norm[ctrl_idx_in_well, , drop = FALSE], na.rm = TRUE)


  mat_effect <- sweep(mat_norm, 2, ctrl_mean_curve)   # cada row less the media ctrl

  curve_mat_list[[i]]  <- mat_effect
  curve_norm_list[[i]] <- mat_norm
  meta_list[[i]] <- data.frame(
    run      = run_id,
    plate    = plate_id,
    well_col = well_cols,
    well_pos = well_pos,
    col_idx  = col_idx_all,
    is_ctrl  = is_ctrl_pos,
    stringsAsFactors = FALSE
  )

  if (i %% 10 == 0) cat(sprintf("  %d/%d plate-runs...\n", i, length(files)))
}

valid_idx      <- !sapply(curve_mat_list, is.null)
curve_mat_all  <- do.call(rbind, curve_mat_list[valid_idx])   # curves effect
curve_norm_all <- do.call(rbind, curve_norm_list[valid_idx])  # curves normalised
meta_all       <- do.call(rbind, meta_list[valid_idx])
rownames(meta_all) <- NULL

cat(sprintf(": %d  ×  %d puntos temporal\n",
            nrow(curve_mat_all), ncol(curve_mat_all)))

## ── Filter NA: remove rows with more of the 20% puntos missing ───────────
na_frac <- rowMeans(is.na(curve_mat_all))
keep    <- na_frac < 0.20
curve_mat_all  <- curve_mat_all[keep, ]
curve_norm_all <- curve_norm_all[keep, ]
meta_all       <- meta_all[keep, ]
cat(sprintf("Curvas after filtered NA: %d\n", nrow(curve_mat_all)))

## Impute NA remaining by interpolation linear simple row row
impute_row <- function(and) {
  if (all(!is.na(and))) return(and)
  idx <- seq_along(and)
  ok  <- !is.na(and)
  if (sum(ok) < 2) return(and)
  approx(idx[ok], and[ok], xout = idx, rule = 2)$and
}
curve_mat_all  <- t(apply(curve_mat_all,  1, impute_row))
curve_norm_all <- t(apply(curve_norm_all, 1, impute_row))

## ── fPCA on curves effect compound ─────────────────────────────────
## Seforr compounds controls
is_ctrl_vec    <- meta_all$is_ctrl
cmpd_mat       <- curve_mat_all[!is_ctrl_vec, ]   # curves effect
ctrl_mat       <- curve_mat_all[is_ctrl_vec,  ]
cmpd_norm_mat  <- curve_norm_all[!is_ctrl_vec, ]  # for figures profiles
ctrl_norm_mat  <- curve_norm_all[is_ctrl_vec,  ]
meta_cmpd      <- meta_all[!is_ctrl_vec, ]
meta_ctrl      <- meta_all[is_ctrl_vec,  ]

cat(sprintf("Curvas compound: %d | Curvas control: %d\n",
            nrow(cmpd_mat), nrow(ctrl_mat)))

## fPCA on the compounds (excluding controls for no contaminate the decomposition)
cat("Ajustando fPCA...\n")
set.seed(2026)
fpca_fit <- refund::fpca.face(
  Y        = cmpd_mat,
  center   = TRUE,
  argvals  = TGRID,
  npc      = N_FPCS,
  pve      = 0.99,      # retener hasta explicar 99% variance
  var      = TRUE
)

## Varianza explained
pve_vec  <- fpca_fit$pve
cum_pve  <- cumsum(pve_vec)
n_pc_current <- length(pve_vec)
cat(sprintf("PCs retenidas: %d\n", n_pc_current))
cat("Varianza explained by PC:", round(pve_vec * 100, 1), "\n")
cat("Acumulada:                 ", round(cum_pve * 100, 1), "\n")

## Guardar summary PCs
pve_df <- data.frame(
  PC           = seq_len(n_pc_current),
  pve_pct      = round(pve_vec * 100, 2),
  cum_pve_pct  = round(cum_pve * 100, 2)
)
write.table(pve_df, file.path(OUT_DIR, "Table_fPCA_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## ── Puntuaciones fPC for all the wells (projection) ────────────────────────
## Proyectar curves control the space fPC computed on compounds
mu_hat    <- fpca_fit$mu
efuncs    <- fpca_fit$efunctions   # [n_tgrid × n_pc]
dt        <- TGRID[2] - TGRID[1]

## Scores for compounds (already computed by fpca.face)
scores_cmpd <- fpca_fit$scores   # [n_cmpd × n_pc]

## Proyeccion manual controls the same space
center_curve <- function(mat, mu) sweep(mat, 2, mu)
ctrl_centered <- center_curve(ctrl_mat, mu_hat)
scores_ctrl <- ctrl_centered %*% efuncs * dt   # integral discreta

## Compilar all the scores with metadata
scores_all <- rbind(
  cbind(meta_cmpd, as.data.frame(scores_cmpd)),
  cbind(meta_ctrl, as.data.frame(scores_ctrl))
)
colnames(scores_all)[(ncol(meta_cmpd)+1):ncol(scores_all)] <-
  paste0("PC", seq_len(ncol(scores_cmpd)))

## ── Join metadata compound (name of the compound) ─────────────────────────
## Res_final_robusto: Chemical, Well, auc, well_pos, run, plate
res_cmpd <- read.table("./data/intermediate/primary_pipeline/final_results_robust.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
res_cmpd$well_pos <- as.integer(res_cmpd$well_pos)
scores_all$well_pos <- as.integer(scores_all$well_pos)

scores_all <- merge(
  scores_all,
  res_cmpd[, c("run", "well_pos", "Chemical")],
  by = c("run", "well_pos"),
  all.x = TRUE
)
scores_all$Chemical[is.na(scores_all$Chemical)] <- "_control_"
scores_all$is_core_hit <- scores_all$Chemical %in% core_hits


cat("Computing statistics per compound (scores fPC medios + AUC_rel)...\n")

## Aggregate per compound: score medio PC1 (and PC2 if it exists)
pc_cols <- paste0("PC", seq_len(n_pc_current))

scores_by_cmpd <- scores_all %>%
  filter(Chemical != "_control_") %>%
  group_by(Chemical) %>%
  summarise(
    across(all_of(pc_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
    n_reps = n(),
    .groups = "drop"
  )

## Join AUC_rel of the pipeline
res_anova_stats <- read.table(
  "./data/intermediate/primary_pipeline/statistics_robust.txt",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
## AUC_rel per compound (media replicates)
auc_by_cmpd <- res_cmpd %>%
  group_by(Chemical) %>%
  summarise(AUC_rel_mean = mean(AUC_relative, na.rm = TRUE), .groups = "drop")

scores_by_cmpd <- merge(scores_by_cmpd, auc_by_cmpd, by = "Chemical", all.x = TRUE)
scores_by_cmpd$is_core_hit <- scores_by_cmpd$Chemical %in% core_hits

## PC1 ~ AUC_rel (Pearson and Spearman)
ok <- !is.na(scores_by_cmpd$PC1_mean) & !is.na(scores_by_cmpd$AUC_rel_mean)
r_pearson  <- cor(scores_by_cmpd$PC1_mean[ok], scores_by_cmpd$AUC_rel_mean[ok],
                  method = "pearson")
r_spearman <- cor(scores_by_cmpd$PC1_mean[ok], scores_by_cmpd$AUC_rel_mean[ok],
                  method = "spearman")
cat(sprintf("\nfPC1 ~ AUC_rel: Pearson r = %.4f  Spearman rho = %.4f\n",
            r_pearson, r_spearman))
cat(sprintf("R² (Pearson) = %.4f — proportion AUC_rel variance explained by fPC1\n",
            r_pearson^2))

## Compilar table summary
fpc_summary <- merge(
  scores_by_cmpd,
  res_anova_stats[, c("Chemical", "FDR_BH")],
  by = "Chemical", all.x = TRUE
)
names(fpc_summary)[names(fpc_summary) == "FDR_BH"] <- "FDR_AUC_ANOVA"
fpc_summary <- fpc_summary[order(fpc_summary$PC1_mean), ]

write.table(fpc_summary, file.path(OUT_DIR, "Table_fPCA_hits.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Tabla saved: Table_fPCA_hits.txt\n")

cat(sprintf("\n=== FDA ===\n"))
cat(sprintf("fPC1 explains %.1f%% of the variance functional\n", pve_vec[1]*100))
cat(sprintf("fPC1 ~ AUC_rel: R² = %.3f (Pearson)\n", r_pearson^2))
cat(sprintf("→ AUC captures ~%.0f%% of the information pharmacological este screen\n",
            r_pearson^2 * pve_vec[1] * 100))


n_show <- min(n_pc_current, 4)
pc_scale <- apply(scores_cmpd[, seq_len(n_show), drop=FALSE], 2, sd, na.rm=TRUE)

fig1_list <- vector("list", n_show + 1)

## Panel A: curve media global
df_mean <- data.frame(t = TGRID, and = mu_hat)
fig1_list[[1]] <- ggplot(df_mean, aes(x = t, and = and)) +
  geom_line(linewidth = 1.1, color = "black") +
  theme_bw(base_size = 13) +
  labs(title = "Global mean curve (normalised OD600)",
       x = "Time (h)", and = "OD600 / median control")

## Paneles B–E: PC1 PC_n_show
for (k in seq_len(n_show)) {
  ef    <- efuncs[, k]
  sc    <- pc_scale[k]
  df_pc <- data.frame(
    t    = TGRID,
    mean = mu_hat,
    plus = mu_hat + sc * ef,
    minus= mu_hat - sc * ef
  )
  pve_k <- round(pve_vec[k] * 100, 1)

  fig1_list[[k + 1]] <- ggplot(df_pc, aes(x = t)) +
    geom_ribbon(aes(ymin = minus, ymax = plus), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(and = mean), linewidth = 0.9) +
    geom_line(aes(and = plus),  color = "steelblue", linetype = "dashed", linewidth = 0.7) +
    geom_line(aes(and = minus), color = "tomato",    linetype = "dashed", linewidth = 0.7) +
    theme_bw(base_size = 13) +
    labs(title = sprintf("fPC%d (%.1f%% variance)", k, pve_k),
         subtitle = "Solid = mean; dashed = ±1 SD of PC scores",
         x = "Time (h)", and = "OD600 / median control")
}

## -- Draft figure (review layout, original dimensions) ----------------
pdf(file.path(OUT_DIR, "Fig_fPCA_components_draft.pdf"), width = 14, height = 10)
do.call(grid.arrange, c(fig1_list, list(ncol = 3)))
dev.off()

## -- Publication figure (PLOS format + TIFF definitive) --------------------
fig1_pub <- cowplot::plot_grid(
  plotlist = fig1_list,
  ncol = 3,
  labels = c("A","B","C","D","E"),
  label_size = 18, label_fontface = "bold"
)
pdf(file.path(OUT_DIR, "Fig_fPCA_components.pdf"), width = 6.83, height = 5)
print(fig1_pub)
dev.off()
cat("Fig_fPCA_components.pdf\n")

## ── Fig 2A: Correlacion fPC1 vs AUC_rel ──────────────────────────────────────
p_corr <- ggplot(fpc_summary[!is.na(fpc_summary$AUC_rel_mean), ],
                  aes(x = PC1_mean, and = AUC_rel_mean,
                      color = is_core_hit)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", formula = and ~ x, se = TRUE,
              color = "steelblue", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick"),
                     labels = c("Non-hit", "Core hit (n=18)"),
                     name = NULL) +
  annotate("text", x = -Inf, and = Inf, hjust = -0.1, vjust = 1.3,
           label = sprintf("R² = %.3f (Pearson)", r_pearson^2),
           size = 3.5, color = "steelblue") +
  theme_bw(base_size = 13) +
  labs(title = sprintf("fPC1 ~ AUC_rel: R² = %.3f", r_pearson^2),
       subtitle = sprintf("fPC1 explains %.1f%% of functional variance; AUC captures ~%.0f%% of pharmacological information",
                          pve_vec[1]*100, r_pearson^2 * pve_vec[1]*100),
       x = sprintf("Mean fPC1 score (%.1f%% of functional variance)", pve_vec[1]*100),
       and = "Mean AUC_rel (B-score + median pipeline)")

p_hist <- ggplot(fpc_summary, aes(x = PC1_mean, fill = is_core_hit)) +
  geom_histogram(binwidth = diff(range(fpc_summary$PC1_mean, na.rm=TRUE))/50,
                 color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                    labels = c("Non-hit", "Core hit"), name = NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  theme_bw(base_size = 13) +
  labs(title = "Distribution of fPC1 scores",
       subtitle = "fPC1 < 0 = growth inhibition relative to plate control mean",
       x = "Mean fPC1 score", and = "Count")

## ── Fig 2C: Scatter fPC1 vs AUC_rel with etiquetas hits ────────────────────
hits_df  <- fpc_summary[fpc_summary$is_core_hit & !is.na(fpc_summary$AUC_rel_mean), ]

p_scatter <- ggplot(fpc_summary[!is.na(fpc_summary$AUC_rel_mean), ],
                     aes(x = PC1_mean, and = AUC_rel_mean)) +
  geom_point(color = "grey50", size = 0.9, alpha = 0.5) +
  geom_point(data = hits_df, color = "firebrick", size = 2.5) +
  geom_text(data = hits_df,
            aes(label = Chemical),
            size = 2.2, vjust = -0.8, color = "firebrick", check_overlap = TRUE) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  theme_bw(base_size = 13) +
  labs(title = "Core hits in functional space",
       subtitle = "Red = 18 core hits; dotted lines = no-effect reference",
       x = "Mean fPC1 score", and = "Mean AUC_rel")

## -- Draft figure (review layout, original dimensions) ----------------
pdf(file.path(OUT_DIR, "Fig_fPCA_scores_draft.pdf"), width = 14, height = 10)
grid.arrange(p_corr, p_hist, p_scatter, ncol = 2)
dev.off()

## -- Publication figure (PLOS format + TIFF definitive) --------------------
scores_pub <- cowplot::plot_grid(
  p_corr, p_hist, p_scatter,
  ncol = 2,
  labels = c("A","B","C"),
  label_size = 18, label_fontface = "bold"
)
pdf(file.path(OUT_DIR, "Fig_fPCA_scores.pdf"), width = 6.83, height = 7)
print(scores_pub)
dev.off()
cat("Fig_fPCA_scores.pdf\n")

## ── Fig 3: Curvas growth the 18 core hits (means between replicates) ───────
cat("Generating profiles of curves core hits...\n")

## Para cada core hit: curve media normalised between replicates
hit_curve_list <- vector("list", length(core_hits))

for (hi in seq_along(core_hits)) {
  chem <- core_hits[hi]
  ## scores_all tiene run + well_pos; meta_all tiene the same columns and same order
  idx_s <- which(scores_all$Chemical == chem & !scores_all$is_ctrl)
  if (length(idx_s) == 0) next
  key_all  <- paste(meta_all$run, meta_all$well_pos)
  key_hit  <- paste(scores_all$run[idx_s], scores_all$well_pos[idx_s])
  rows_in_mat <- match(key_hit, key_all)
  rows_in_mat <- rows_in_mat[!is.na(rows_in_mat)]
  if (length(rows_in_mat) == 0) next

  ## Use curves normalised (no effect) for visualisation
  chem_curves <- curve_norm_all[rows_in_mat, , drop = FALSE]
  mean_curve  <- colMeans(chem_curves, na.rm = TRUE)
  sd_curve    <- apply(chem_curves, 2, sd, na.rm = TRUE)

  hit_curve_list[[hi]] <- data.frame(
    Chemical = chem,
    t        = TGRID,
    y_mean   = mean_curve,
    y_sd     = sd_curve,
    n_reps   = nrow(chem_curves)
  )
}

ctrl_rows_mat    <- which(meta_all$is_ctrl)
ctrl_mean_global <- colMeans(curve_norm_all[ctrl_rows_mat, ], na.rm = TRUE)
ctrl_sd_global   <- apply(curve_norm_all[ctrl_rows_mat, ], 2, sd, na.rm = TRUE)
ctrl_df <- data.frame(
  Chemical = "Control",
  t = TGRID, y_mean = ctrl_mean_global, y_sd = ctrl_sd_global, n_reps = length(ctrl_rows_mat)
)

all_hit_curves <- rbind(do.call(rbind, Filter(Negate(is.null), hit_curve_list)), ctrl_df)

PAL_HIT <- scales::hue_pal()(length(core_hits))
names(PAL_HIT) <- core_hits
PAL_ALL <- c(PAL_HIT, "Control" = "black")

p_profiles <- ggplot(all_hit_curves, aes(x = t, and = y_mean, color = Chemical)) +
  geom_line(data = all_hit_curves[all_hit_curves$Chemical == "Control", ],
            linewidth = 1.2, linetype = "dashed") +
  geom_line(data = all_hit_curves[all_hit_curves$Chemical != "Control", ],
            linewidth = 0.7, alpha = 0.85) +
  scale_color_manual(values = PAL_ALL, guide = "none") +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 9) +
  labs(title = "Growth curve profiles: 18 core hits vs controls (mean across replicates)",
       subtitle = "All curves normalised by plate-run control median; dashed = control mean",
       x = "Time (h)", and = "OD600 / median control")

## Faceted by hit for greater clarity
p_profiles_facet <- ggplot(all_hit_curves[all_hit_curves$Chemical != "Control", ],
                             aes(x = t, and = y_mean)) +
  geom_ribbon(aes(ymin = pmax(y_mean - y_sd, 0),
                  ymax = y_mean + y_sd),
              alpha = 0.2, fill = "tomato") +
  geom_line(color = "tomato", linewidth = 0.8) +
  geom_line(data = ctrl_df, aes(x = t, and = y_mean),
            color = "black", linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ Chemical, ncol = 4, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9)) +
  labs(title = "Core hit growth profiles (red = hit ± SD; dashed = control mean)",
       x = "Time (h)", and = "OD600 / median control")

pdf(file.path(OUT_DIR, "Fig_fPCA_hit_profiles_draft.pdf"), width = 14, height = 12)
print(p_profiles_facet)
dev.off()
pdf(file.path(OUT_DIR, "Fig_fPCA_hit_profiles.pdf"), width = 6.83, height = 8)
print(p_profiles_facet)
dev.off()
file.copy(file.path(OUT_DIR, "Table_fPCA_summary.txt"),
          file.path(DEFS_DIR, "S4_Table_fPCA_pve.txt"), overwrite = TRUE)
file.copy(file.path(OUT_DIR, "Table_fPCA_hits.txt"),
          file.path(DEFS_DIR, "S5_Table_fPCA_hits.txt"), overwrite = TRUE)
cat("saved: Fig_fPCA_hit_profiles.pdf\n")

## ── Summary final ─────────────────────────────────────────────────────────────
cat("\n=== Summary FDA ===\n")
cat(sprintf("%d compound + %d control\n",
            nrow(cmpd_mat), nrow(ctrl_mat)))
cat(sprintf("PCs %d (explican %.1f%% variance)\n",
            n_pc_current, max(cum_pve)*100))
cat(sprintf("fPC1 ~ AUC_rel: Pearson R = %.4f, R² = %.4f, Spearman rho = %.4f\n",
            r_pearson, r_pearson^2, r_spearman))
cat(sprintf("Conclusion: fPC1 (%.1f%% var. functional) ≡ inhibition amplitud (≡ AUC)\n",
            pve_vec[1]*100))
            sum(fpc_summary$is_core_hit)))

cat("\n=== 04_functional_data_analysis.R completedo ===\n")
