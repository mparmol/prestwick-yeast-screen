## ============================================================
## 01_zprime_ssmd.R — A.4 of the analytical plan
## Plate quality metrics: Z'-factor, SSMD, RZ'-factor, CV
## Inputs : pipeline results (data/intermediate/primary_pipeline/)
## Outputs: Table_QC_zprime_ssmd.txt + Fig_QC_zprime_ssmd.pdf
## ============================================================

## ── Anchor working directory at the repository root ─────────────────────────
## Works whether launched from the repo root or any subfolder
.find_root <- function() {
  d <- normalizePath(getwd())
  for (i in seq_len(6)) {
    if (dir.exists(file.path(d, "data/raw"))) return(d)
    d <- dirname(d)
  }
  stop("Could not find the project root (data/raw not found). ",
       "Run setwd() to the repo root before source().")
}
setwd(.find_root())
cat(sprintf("Working directory: %s\n", getwd()))

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)

OUT_DIR  <- "./data/intermediate/validation"
DEFS_DIR <- file.path(OUT_DIR, "definitive")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DEFS_DIR, showWarnings = FALSE, recursive = TRUE)

## ── Load pipeline results ────────────────────────────────────────
res_cmpd  <- read.table("./data/intermediate/primary_pipeline/final_results_robust.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
res_ctrl  <- read.table("./data/intermediate/primary_pipeline/control_results_robust.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
res_lmm   <- read.table("./data/intermediate/primary_pipeline/statistics_LMM.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
res_anova <- read.table("./data/intermediate/primary_pipeline/statistics_robust.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## ── Core hits: ANOVA ∩ LMM convergence (FDR < 0.05 in both) ────────────
## Chicago sky blue 6B excluded priori: colorimetric artefact (abs. ~620 nm)
EXCLUDE_ARTEFACTS <- "Chicago sky blue 6B"

hits_lmm   <- res_lmm[res_lmm$FDR_BH  < 0.05, "Chemical"]
hits_anova <- res_anova[res_anova$FDR_BH < 0.05, "Chemical"]
core_hits  <- setdiff(intersect(hits_lmm, hits_anova), EXCLUDE_ARTEFACTS)
cat(sprintf("Core hits (ANOVA ∩ LMM, FDR<0.05, artefacts excluded): %d\n",
            length(core_hits)))
cat(paste0("  ", core_hits, collapse = "\n"), "\n\n")

## ── Functions metrics ────────────────────────────────────────────────────
z_prime <- function(pos, neg) {
  if (length(pos) < 2 || length(neg) < 2) return(NA_real_)
  denom <- abs(mean(pos, na.rm = TRUE) - mean(neg, na.rm = TRUE))
  if (denom < 1e-9) return(NA_real_)
  1 - 3 * (sd(pos, na.rm = TRUE) + sd(neg, na.rm = TRUE)) / denom
}

ssmd_fn <- function(pos, neg) {
  if (length(pos) < 2 || length(neg) < 2) return(NA_real_)
  denom <- sqrt(var(pos, na.rm = TRUE) + var(neg, na.rm = TRUE))
  if (denom < 1e-9) return(NA_real_)
  (mean(neg, na.rm = TRUE) - mean(pos, na.rm = TRUE)) / denom
}

rz_prime <- function(pos, neg) {
  if (length(pos) < 2 || length(neg) < 2) return(NA_real_)
  denom <- abs(median(pos, na.rm = TRUE) - median(neg, na.rm = TRUE))
  if (denom < 1e-9) return(NA_real_)
  1 - 3 * (mad(pos, na.rm = TRUE) + mad(neg, na.rm = TRUE)) / denom
}

## ── Loop by plate-run ──────────────────────────────────────────────────────
plate_runs <- sort(unique(res_ctrl$run))
qc_list    <- vector("list", length(plate_runs))

for (i in seq_along(plate_runs)) {
  pr       <- plate_runs[i]
  pl       <- unique(res_ctrl[res_ctrl$run == pr, "plate"])

  neg_aucs <- res_ctrl[res_ctrl$run == pr, "auc"]
  hit_rows <- res_cmpd[res_cmpd$Chemical %in% core_hits & res_cmpd$run == pr, ]
  pos_aucs <- hit_rows$auc   # raw AUC (post B-score) the hits this run

  n_hits    <- nrow(hit_rows)
  hit_names <- if (n_hits > 0) paste(unique(hit_rows$Chemical), collapse = "; ") else "—"

  qc_list[[i]] <- data.frame(
    run          = pr,
    plate        = pl,
    n_ctrl       = length(neg_aucs),
    mean_ctrl    = round(mean(neg_aucs),  3),
    sd_ctrl      = round(sd(neg_aucs),    3),
    cv_ctrl_pct  = round(100 * sd(neg_aucs) / mean(neg_aucs), 2),
    n_hits       = n_hits,
    zprime       = round(z_prime(pos_aucs, neg_aucs),  3),
    ssmd         = round(ssmd_fn(pos_aucs, neg_aucs),  3),
    rzprime      = round(rz_prime(pos_aucs, neg_aucs), 3),
    hits_in_run  = hit_names,
    stringsAsFactors = FALSE
  )
}

qc_df <- do.call(rbind, qc_list)
qc_df <- qc_df[order(qc_df$plate), ]
rownames(qc_df) <- NULL

cat("=== QC metrics by plate-run ===\n")
print(qc_df[, c("run","cv_ctrl_pct","n_hits","zprime","ssmd","rzprime")],
      row.names = FALSE)
cat(sprintf("\nMediana Z' (runs with hits): %.3f  [rango %.3f–%.3f]\n",
            median(qc_df$zprime,  na.rm=TRUE),
            min(qc_df$zprime,   na.rm=TRUE),
            max(qc_df$zprime,   na.rm=TRUE)))
cat(sprintf("Mediana SSMD (runs with hits): %.3f\n",
            median(qc_df$ssmd, na.rm=TRUE)))
cat(sprintf("Runs with Z' > 0.5: %d de %d with hits\n",
            sum(qc_df$zprime > 0.5, na.rm=TRUE),
            sum(!is.na(qc_df$zprime))))

write.table(qc_df, file.path(OUT_DIR, "Table_QC_zprime_ssmd.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nTable_QC_zprime_ssmd.txt\n")

qc_df$run_label <- factor(qc_df$run,
                           levels = qc_df$run[order(qc_df$plate)])

p_cv <- ggplot(qc_df, aes(x = run_label, and = cv_ctrl_pct,
                            fill = cv_ctrl_pct > 20)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                    guide = "none") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 0.7) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  labs(title = "Control CV per plate-run", x = NULL, and = "CV (%)")

p_zprime <- ggplot(qc_df[!is.na(qc_df$zprime), ],
                    aes(x = run_label, and = zprime, fill = zprime > 0.5)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "orange", "TRUE" = "steelblue"),
                    guide = "none") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgreen", linewidth = 0.9) +
  geom_hline(yintercept = 0,   linetype = "dotted", color = "red",       linewidth = 0.7) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  labs(title = "Z'-factor per plate-run (core hits as positive controls)",
       subtitle = "Green dashed = Z' 0.5 (excellent); plate-runs without hits omitted",
       x = NULL, and = "Z'-factor")

p_ssmd <- ggplot(qc_df[!is.na(qc_df$ssmd), ],
                  aes(x = run_label, and = ssmd, fill = ssmd > 3)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "orange", "TRUE" = "steelblue"),
                    guide = "none") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "darkgreen", linewidth = 0.9) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  labs(title = "SSMD per plate-run",
       subtitle = "Green dashed = SSMD 3 (excellent; Zhang 2007)",
       x = NULL, and = "SSMD")

p_rzprime <- ggplot(qc_df[!is.na(qc_df$rzprime), ],
                     aes(x = run_label, and = rzprime, fill = rzprime > 0.5)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "orange", "TRUE" = "steelblue"),
                    guide = "none") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgreen", linewidth = 0.9) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  labs(title = "Robust Z'-factor (RZ') per plate-run",
       subtitle = "Uses median/MAD instead of mean/SD",
       x = NULL, and = "RZ'-factor")

## -- Draft figure (review layout, original dimensions) ----------------
pdf(file.path(OUT_DIR, "Fig_QC_zprime_ssmd_draft.pdf"), width = 14, height = 12)
grid.arrange(p_cv, p_zprime, p_ssmd, p_rzprime, ncol = 2)
dev.off()

## -- Publication figure (PLOS format + TIFF definitive) --------------------
qc_pub <- cowplot::plot_grid(
  p_cv, p_zprime, p_ssmd, p_rzprime,
  ncol = 2,
  labels = c("A","B","C","D"),
  label_size = 18, label_fontface = "bold"
)
pdf(file.path(OUT_DIR, "Fig_QC_zprime_ssmd.pdf"), width = 6.83, height = 7)
print(qc_pub)
dev.off()
file.copy(file.path(OUT_DIR, "Table_QC_zprime_ssmd.txt"),
          file.path(DEFS_DIR, "S1_Table_QC.txt"), overwrite = TRUE)
cat("Fig_QC_zprime_ssmd.pdf\n")
cat("\n=== 01_zprime_ssmd.R completedo ===\n")
