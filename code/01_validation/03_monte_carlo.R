## ============================================================
## 03_monte_carlo.R — A.3 of the analytical plan
## ¡Monte Carlo for validar FDR and sensitivity of the pipeline
##
## Varibefore evaluadas:
##   V1: single-centroid  (without B-score, referencia = 1 well DTW)
##   V2: median-ctrl      (without B-score, referencia = median controls)
##   V3: bscore + median  (positional correction + median controls)
## ── Anchor working directory at the repository root ─────────────────────────
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
##   V4: bscore + median + doble-criterio (ANOVA < 0.05 & F_robust < 0.05)
##
##
## Outputs: Table_MC_summary.txt + Fig_MC_results.pdf + Fig_MC_power.pdf
## ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)

## pROC without library() for evitar that pROC::auc() enmascare oafter functions
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC", quiet = TRUE)

set.seed(2026)

OUT_DIR  <- "./data/intermediate/validation"
DEFS_DIR <- file.path(OUT_DIR, "definitive")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DEFS_DIR, showWarnings = FALSE, recursive = TRUE)

N_ITER      <- 1000     # iterations Monte Carlo
N_PLATES    <- 17       # layouts placa
N_REPS      <- 3        # replicates biological
N_CMPDS     <- 80       # compounds by placa
N_CTRL      <- 16       # controls by plate-run
MEAN_PLATE  <- 44       # median AUC controls (unidades reales of the dataset)
SD_PLATE    <- 8        # variabilidad between plates
SD_WITHIN   <- 5.3      # noise intra-plate (SD controls ≈ 5.3 observado)
SD_POS      <- 1.5      # magnitude effects position row/col
N_HIT_PLATE <- 5        # hits true by placa
TRUE_EFFECTS <- c(0.10, 0.30, 0.50, 0.70, 0.90)  # fraccion AUC retenida
FDR_THRESH  <- 0.05

bscore_vec <- function(vals, n_rows = 8, n_cols = 10) {
  mat <- matrix(vals, nrow = n_rows, ncol = n_cols)
  mp  <- tryCatch(medpolish(mat, trace.iter = FALSE, na.rm = TRUE),
                  error = function(and) NULL)
  if (is.null(mp)) return(vals)
  as.vector(mp$overall + mp$residuals)
}

simulate_screen <- function(true_eff = TRUE_EFFECTS) {
  cmpd_rows <- vector("list", N_PLATES * N_REPS)
  ctrl_rows <- vector("list", N_PLATES * N_REPS)
  idx <- 1L
  hit_registry <- character(0)

  for (p in seq_len(N_PLATES)) {
    plate_mu <- rnorm(1, MEAN_PLATE, SD_PLATE)
    hit_pos  <- sample(N_CMPDS, N_HIT_PLATE)
    eff_asgn <- sample(true_eff, N_HIT_PLATE, replace = TRUE)
    cmpd_ids <- paste0("P", p, "_C", seq_len(N_CMPDS))
    hit_registry <- c(hit_registry, cmpd_ids[hit_pos])

    row_eff <- rnorm(8,  0, SD_POS)
    col_eff <- rnorm(10, 0, SD_POS)
    pos_eff <- as.vector(outer(row_eff, col_eff))   # length N_CMPDS

    for (r in seq_len(N_REPS)) {
      ctrl_auc  <- rnorm(N_CTRL, plate_mu, SD_WITHIN)
      cmpd_base <- plate_mu + pos_eff
      cmpd_auc  <- cmpd_base + rnorm(N_CMPDS, 0, SD_WITHIN)
      cmpd_auc[hit_pos] <- plate_mu * eff_asgn + rnorm(N_HIT_PLATE, 0, SD_WITHIN * 0.5)

      cmpd_rows[[idx]] <- data.frame(
        plate    = p, rep = r,
        compound = cmpd_ids,
        auc_raw  = cmpd_auc,
        true_hit = seq_len(N_CMPDS) %in% hit_pos,
        well_pos = seq_len(N_CMPDS),
        stringsAsFactors = FALSE
      )
      ctrl_rows[[idx]] <- data.frame(
        plate   = p, rep = r,
        auc_raw = ctrl_auc,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  list(cmpd = bind_rows(cmpd_rows),
       ctrl = bind_rows(ctrl_rows),
       hit_ids = unique(hit_registry))
}


run_pipeline <- function(sim, variant = c("V1","V2","V3","V4")) {
  variant <- match.arg(variant)
  cmpd <- sim$cmpd
  ctrl <- sim$ctrl

  norm_cmpd <- vector("list", N_PLATES * N_REPS)
  norm_ctrl <- vector("list", N_PLATES * N_REPS)
  k <- 1L

  for (p in seq_len(N_PLATES)) {
    for (r in seq_len(N_REPS)) {
      cpr <- cmpd[cmpd$plate == p & cmpd$rep == r, ]
      ctr <- ctrl[ctrl$plate == p & ctrl$rep == r, ]
      neg <- ctr$auc_raw
      med <- median(neg)
      ref <- if (variant == "V1") neg[1L] else med   # centroide = 1er well

      cmpd_norm_v <- switch(variant,
        V1 = cpr$auc_raw / max(ref, 1e-6),
        V2 = cpr$auc_raw / max(med, 1e-6),
        V3 = bscore_vec(cpr$auc_raw) / max(med, 1e-6),
        V4 = bscore_vec(cpr$auc_raw) / max(med, 1e-6)
      )
      ctrl_norm_v <- neg / max(ref, 1e-6)   # controls: same that compounds
      if (variant %in% c("V2","V3","V4")) ctrl_norm_v <- neg / max(med, 1e-6)

      norm_cmpd[[k]] <- data.frame(
        compound = cpr$compound, plate = p, rep = r,
        true_hit = cpr$true_hit, auc_norm = cmpd_norm_v,
        stringsAsFactors = FALSE
      )
      norm_ctrl[[k]] <- data.frame(
        plate = p, rep = r, auc_norm = ctrl_norm_v,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  nc_df  <- bind_rows(norm_cmpd)
  nct_df <- bind_rows(norm_ctrl)

  compounds <- unique(nc_df$compound)
  n_cmpds   <- length(compounds)
  pval_main <- numeric(n_cmpds)
  pval_withf <- numeric(n_cmpds)   # only is used V4
  score_v   <- numeric(n_cmpds)
  is_true   <- logical(n_cmpds)

  for (i in seq_along(compounds)) {
    ci <- nc_df[nc_df$compound == compounds[i], ]
    pl <- ci$plate[1L]

    ct_all  <- nct_df[nct_df$plate == pl, "auc_norm"]
    ct_12   <- nct_df[nct_df$plate == pl & nct_df$rep %in% 1:2, "auc_norm"]
    ct_3    <- nct_df[nct_df$plate == pl & nct_df$rep == 3L,    "auc_norm"]

    cmpd_all <- ci$auc_norm             # 3 valuees
    cmpd_12  <- ci[ci$rep %in% 1:2, "auc_norm"]
    cmpd_3   <- ci[ci$rep == 3L,    "auc_norm"]

    pval_main[i] <- tryCatch(
      anova(lm(c(cmpd_all, ct_all) ~
                c(rep("cmpd", length(cmpd_all)),
                  rep("ctrl", length(ct_all)))))$Pr[1L],
      error = function(and) 1
    )

    if (variant == "V4") {
      pval_withf[i] <- tryCatch(
        anova(lm(c(cmpd_12, ct_12) ~
                  c(rep("cmpd", length(cmpd_12)),
                    rep("ctrl", length(ct_12)))))$Pr[1L],
        error = function(and) 1
      )
    }

    score_v[i] <- -log10(max(pval_main[i], 1e-300))
    is_true[i] <- any(ci$true_hit)
  }

  ## BH-FDR and hit calling
  fdr_main <- p.adjust(pval_main, method = "BH")

  called <- if (variant == "V4") {
    fdr_withf <- p.adjust(pval_withf, method = "BH")
    fdr_main < FDR_THRESH & fdr_withf < FDR_THRESH
  } else {
    fdr_main < FDR_THRESH
  }

  data.frame(compound = compounds, true_hit = is_true,
              called = called, score = score_v,
              fdr_main = fdr_main,
              variant = variant, stringsAsFactors = FALSE)
}

## ── Monte Carlo principal ─────────────────────────────────────────────────────
cat(sprintf("Iniciando %d iterations × 4 variants...\n", N_ITER))
mc_rows <- vector("list", N_ITER * 4)
k <- 1L

for (it in seq_len(N_ITER)) {
  sim <- simulate_screen()
  for (v in c("V1","V2","V3","V4")) {
    res <- tryCatch(run_pipeline(sim, variant = v), error = function(and) NULL)
    if (is.null(res)) { k <- k + 1L; next }

    n_true   <- sum(res$true_hit)
    n_called <- sum(res$called, na.rm = TRUE)
    n_tp     <- sum(res$called & res$true_hit, na.rm = TRUE)
    n_fp     <- sum(res$called & !res$true_hit, na.rm = TRUE)

    auroc <- tryCatch({
      ro <- pROC::roc(res$true_hit, res$score, quiet = TRUE)
      as.numeric(pROC::auc(ro))
    }, error = function(and) NA_real_)

    mc_rows[[k]] <- data.frame(
      iter     = it, variant = v,
      sens     = if (n_true   > 0) n_tp / n_true else NA_real_,
      fdr_obs  = if (n_called > 0) n_fp / n_called else 0,
      auroc    = auroc,
      n_tp = n_tp, n_fp = n_fp, n_called = n_called,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  if (it %% 200 == 0) cat(sprintf("  iter %d/%d\n", it, N_ITER))
}

mc_df <- bind_rows(mc_rows)

## ── Summary ──────────────────────────────────────────────────────────────────
VARIANT_LABELS <- c(
  V1 = "V1: Single centroid",
  V2 = "V2: Median ctrl",
  V3 = "V3: B-score + Median",
  V4 = "V4: B-score + Med + dual-criterion"
)

summary_mc <- mc_df %>%
  group_by(variant) %>%
  summarise(
    sens_mean  = round(mean(sens,    na.rm=TRUE), 3),
    sens_sd    = round(sd(sens,      na.rm=TRUE), 3),
    fdr_mean   = round(mean(fdr_obs, na.rm=TRUE), 3),
    fdr_sd     = round(sd(fdr_obs,   na.rm=TRUE), 3),
    auroc_mean = round(mean(auroc,   na.rm=TRUE), 3),
    auroc_sd   = round(sd(auroc,     na.rm=TRUE), 3),
    .groups = "drop"
  )
summary_mc$label <- VARIANT_LABELS[summary_mc$variant]

cat("\n=== Summary Monte Carlo ===\n")
print(summary_mc[, c("label","sens_mean","sens_sd","fdr_mean","fdr_sd",
                       "auroc_mean","auroc_sd")], row.names = FALSE)

write.table(mc_df,      file.path(OUT_DIR, "MC_raw_results.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(summary_mc, file.path(OUT_DIR, "Table_MC_summary.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

## ── Figuras ──────────────────────────────────────────────────────────────────
mc_df$label <- factor(VARIANT_LABELS[mc_df$variant],
                       levels = VARIANT_LABELS[c("V1","V2","V3","V4")])

PAL4 <- c("V1: Single centroid"               = "#AAAAAA",
           "V2: Median ctrl"                   = "#6699CC",
           "V3: B-score + Median"              = "#FF9933",
           "V4: B-score + Med + dual-criterion"= "#2CA02C")

p_sens <- ggplot(mc_df, aes(x=label, and=sens, fill=label)) +
  geom_boxplot(outlier.size=0.4, outlier.alpha=0.3, notch=TRUE) +
  scale_fill_manual(values=PAL4, guide="none") +
  geom_hline(yintercept=0.8, linetype="dashed", color="darkgreen", linewidth=0.8) +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=25, hjust=1, size=11)) +
  labs(title="Sensitivity (true hit recovery)",
       subtitle=sprintf("%d MC iterations | FDR threshold = %.2f | dashed = 80%%",
                         N_ITER, FDR_THRESH),
       x=NULL, and="Sensitivity")

p_fdr <- ggplot(mc_df, aes(x=label, and=fdr_obs, fill=label)) +
  geom_boxplot(outlier.size=0.4, outlier.alpha=0.3, notch=TRUE) +
  scale_fill_manual(values=PAL4, guide="none") +
  geom_hline(yintercept=FDR_THRESH, linetype="dashed", color="red", linewidth=0.8) +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=25, hjust=1, size=11)) +
  labs(title="Empirical FDR",
       subtitle=sprintf("Red dashed = nominal FDR = %.2f", FDR_THRESH),
       x=NULL, and="Observed FDR")

p_auroc <- ggplot(mc_df, aes(x=label, and=auroc, fill=label)) +
  geom_boxplot(outlier.size=0.4, outlier.alpha=0.3, notch=TRUE) +
  scale_fill_manual(values=PAL4, guide="none") +
  coord_cartesian(ylim=c(0.4, 1)) +
  geom_hline(yintercept=0.5, linetype="dotted", color="gray50", linewidth=0.6) +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=25, hjust=1, size=11)) +
  labs(title="AUROC",
       subtitle="Dotted = random classifier",
       x=NULL, and="AUROC")

## Tabla summary as panel D (without tableGrob)
summary_text <- paste0(
  sprintf("%-42s  Sens: %.3f±%.3f  FDR: %.3f±%.3f  AUROC: %.3f±%.3f\n",
           summary_mc$label,
           summary_mc$sens_mean, summary_mc$sens_sd,
           summary_mc$fdr_mean,  summary_mc$fdr_sd,
           summary_mc$auroc_mean,summary_mc$auroc_sd),
  collapse=""
)
p_tbl <- ggplot() + theme_void() +
  annotate("text", x=0.05, and=0.6, label=summary_text, hjust=0, vjust=1,
           family="mono", size=3.5) +
  annotate("text", x=0.05, and=0.9, label="Summary statistics",
           hjust=0, fontface="bold", size=4.5) +
  xlim(0,1) + ylim(0,1)

## -- Draft figure (review layout, original dimensions) ----------------
pdf(file.path(OUT_DIR, "Fig_MC_results_draft.pdf"), width = 14, height = 12)
grid.arrange(p_sens, p_fdr, p_auroc, p_tbl, ncol = 2)
dev.off()

## -- Publication figure (PLOS format + TIFF definitive) --------------------
mc_pub <- cowplot::plot_grid(
  p_sens, p_fdr, p_auroc, p_tbl,
  ncol = 2,
  labels = c("A","B","C","D"),
  label_size = 18, label_fontface = "bold"
)
pdf(file.path(OUT_DIR, "Fig_MC_results.pdf"), width = 6.83, height = 7)
print(mc_pub)
dev.off()
## ── Curvas poder vs magnitude effect ────────────────────────────────────
cat("\ncurves poder (V1 vs V4)...\n")
N_ITER_PWR  <- 300
effect_grid <- seq(0.1, 0.9, by=0.1)

pow_rows <- vector("list", N_ITER_PWR * length(effect_grid) * 2)
k <- 1L

for (eff in effect_grid) {
  for (it in seq_len(N_ITER_PWR)) {
    sim <- simulate_screen(true_eff = rep(eff, N_HIT_PLATE))
    for (v in c("V1","V4")) {
      res <- tryCatch(run_pipeline(sim, variant=v), error=function(and) NULL)
      if (is.null(res)) { k <- k + 1L; next }
      n_true <- sum(res$true_hit)
      sens_v <- if (n_true > 0)
        sum(res$called & res$true_hit, na.rm=TRUE) / n_true else NA_real_
      pow_rows[[k]] <- data.frame(effect=eff, iter=it, variant=v, sens=sens_v,
                                    stringsAsFactors=FALSE)
      k <- k + 1L
    }
  }
  cat(sprintf("  effect=%.1f completedo\n", eff))
}

pow_df <- bind_rows(pow_rows)
pow_df$label <- factor(VARIANT_LABELS[pow_df$variant],
                         levels = VARIANT_LABELS[c("V1","V4")])

pow_summary <- pow_df %>%
  group_by(effect, label) %>%
  summarise(mean_sens=mean(sens, na.rm=TRUE), sd_sens=sd(sens, na.rm=TRUE),
             .groups="drop")

p_power <- ggplot(pow_summary,
                   aes(x=effect, and=mean_sens, color=label, fill=label, group=label)) +
  geom_ribbon(aes(ymin=pmax(mean_sens-sd_sens,0), ymax=pmin(mean_sens+sd_sens,1)),
               alpha=0.15, color=NA) +
  geom_line(linewidth=1.1) +
  geom_point(size=2.2) +
  scale_color_manual(values=PAL4[c("V1: Single centroid",
                                     "V4: B-score + Med + dual-criterion")]) +
  scale_fill_manual(values=PAL4[c("V1: Single centroid",
                                    "V4: B-score + Med + dual-criterion")]) +
  geom_hline(yintercept=0.8, linetype="dashed", color="gray40") +
  scale_x_continuous(breaks=effect_grid,
                      labels=function(x) sprintf("%.1f", x)) +
  coord_cartesian(ylim=c(0,1)) +
  theme_bw(base_size=13) +
  theme(legend.position="bottom") +
  labs(title="Power curves: sensitivity vs. effect size",
       subtitle=sprintf("%d MC iterations per effect level; ribbon = ±1 SD", N_ITER_PWR),
       x="True hit AUC (fraction of control retained; 0.10 = strong inhibitor)",
       and="Mean sensitivity", color=NULL, fill=NULL)

pdf(file.path(OUT_DIR, "Fig_MC_power_curves_draft.pdf"), width = 10, height = 6)
print(p_power)
dev.off()
pdf(file.path(OUT_DIR, "Fig_MC_power_curves.pdf"), width = 6.83, height = 4.5)
print(p_power)
dev.off()
write.table(pow_df, file.path(OUT_DIR, "MC_power_raw.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)
file.copy(file.path(OUT_DIR, "Table_MC_summary.txt"),
          file.path(DEFS_DIR, "S3_Table_MonteCarlo.txt"), overwrite = TRUE)

cat("\nOutputs saved en:", OUT_DIR, "\n")
cat("  MC_raw_results.txt\n")
cat("  Table_MC_summary.txt\n")
cat("  Fig_MC_results.pdf\n")
cat("  Fig_MC_power_curves.pdf\n")
cat("  MC_power_raw.txt\n")
cat("\n=== 03_monte_carlo.R completedo ===\n")
