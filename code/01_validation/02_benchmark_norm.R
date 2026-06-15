## ============================================================
## 02_benchmark_norm.R — A.2 of the analytical plan
## Benchmark cabeza-a-cabeza 6 estrategias normalisation
##   (1) raw         — AUC without norm
##   (2) poc         — AUC / mean(ctrl)  [percent-of-control]
##   (3) rzscore     — (AUC - median(ctrl)) / MAD(ctrl)  [robust Z]
##   (4) dtw_centroid— AUC / AUC_centroid  [centroide DTW, baseline anterior]

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
##   (5) median      — AUC / median(ctrl)
##   (6) bscore_med  — B-score(AUC) / median(ctrl)  [pipeline final]
## Metrics QC: SD_ctrl_norm, Z'-factor, ICC between triplicados, AUROC
## Outputs: Table_benchmark_norm.txt + Fig_benchmark_norm.pdf (new Fig 5)
## ============================================================

library(gcplyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)

## pROC for AUROC; irr for ICC — NO are called with library() for evitar
## that pROC::auc() enmascare gcplyr::auc()
for (pkg in c("pROC", "irr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, quiet = TRUE)
}

OUT_DIR  <- "./data/intermediate/validation"
DEFS_DIR <- file.path(OUT_DIR, "definitive")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DEFS_DIR, showWarnings = FALSE, recursive = TRUE)

PLATES_EXCLUDE <- c(1, 18)
TIME_MAX_H     <- 60

## ── Functions auxiliares ─────────────────────────────────────────────────────
plate_num_from_file <- function(fname)
  as.integer(sub("plate_(\\d+)_.*", "\\1", fname))

well_to_rowcol <- function(well_pos) {
  wp  <- as.integer(well_pos)
  data.frame(row_idx = ceiling(wp / 12),
             col_idx = ((wp - 1L) %% 12L) + 1L)
}

bscore_plate <- function(auc_df) {
  rc  <- well_to_rowcol(auc_df$well_pos)
  auc_df$row_idx <- rc$row_idx
  auc_df$col_idx <- rc$col_idx
  cmpd <- auc_df[auc_df$col_idx %in% 2:11, ]
  if (nrow(cmpd) < 4) return(auc_df)
  mat  <- matrix(NA, 8, 10)
  for (k in seq_len(nrow(cmpd))) {
    r <- cmpd$row_idx[k]
    c <- cmpd$col_idx[k] - 1L
    if (r >= 1 && r <= 8 && c >= 1 && c <= 10) mat[r, c] <- cmpd$auc[k]
  }
  mp <- tryCatch(medpolish(mat, na.rm = TRUE, trace.iter = FALSE),
                 error = function(and) NULL)
  if (is.null(mp)) return(auc_df)
  for (k in seq_len(nrow(cmpd))) {
    r <- cmpd$row_idx[k]; c <- cmpd$col_idx[k] - 1L
    if (r >= 1 && r <= 8 && c >= 1 && c <= 10) {
      idx <- which(auc_df$well_pos == cmpd$well_pos[k])
      auc_df$auc[idx] <- auc_df$auc[idx] -
        ifelse(is.na(mp$row[r]), 0, mp$row[r]) -
        ifelse(is.na(mp$col[c]), 0, mp$col[c])
    }
  }
  auc_df
}

z_prime_fn <- function(pos, neg) {
  if (length(pos) < 2 || length(neg) < 2) return(NA_real_)
  d <- abs(mean(pos, na.rm=TRUE) - mean(neg, na.rm=TRUE))
  if (d < 1e-9) return(NA_real_)
  1 - 3*(sd(pos, na.rm=TRUE) + sd(neg, na.rm=TRUE)) / d
}

my_design <- read.csv(".", sep="\t", encoding="latin1")
files_all  <- list.files("./data/raw/")
files      <- files_all[!sapply(files_all, function(f)
                plate_num_from_file(f) %in% PLATES_EXCLUDE)]

controls_dtw <- read.table(
  "./data/",
  sep="\t", header=TRUE, stringsAsFactors=FALSE)

## Mapear "X10.A.37" → well_pos = 37, for cada run
parse_centroid_pos <- function(repr_str) {
  parts <- strsplit(repr_str, "\\.")[[1]]
  as.integer(parts[3])
}
controls_dtw$centroid_pos <- sapply(controls_dtw$Representante,
                                      parse_centroid_pos)
## Normalizar name run for doesr join
controls_dtw$run_key <- gsub("^placa_", "plate_",
                               tolower(controls_dtw$Placa))

## Core hits (gold standard for AUROC)
res_lmm   <- read.table("./data/intermediate/primary_pipeline/statistics_LMM.txt",
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)
res_anova <- read.table("./data/intermediate/primary_pipeline/statistics_robust.txt",
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)
## Chicago sky blue 6B excluded: artefacto colorimetrico (abs. ~620 nm)
EXCLUDE_ARTEFACTS <- "Chicago sky blue 6B"
core_hits <- setdiff(
  intersect(res_lmm[res_lmm$FDR_BH  < 0.05, "Chemical"],
            res_anova[res_anova$FDR_BH < 0.05, "Chemical"]),
  EXCLUDE_ARTEFACTS
)
cat(sprintf("Gold-standard core hits (artefacts excluded): %d\n", length(core_hits)))

## ── Loop principal: compute AUC and 6 normalizaciones by plate-run ──────────
cat("Procthatndo", length(files), "plate-runs...\n")

all_norms <- NULL   # long-format final

for (i in seq_along(files)) {
  fname    <- files[i]
  plate_id <- plate_num_from_file(fname)
  run_id   <- sub("\\.csv$", "", fname)

  ## ── Carga raw ────────────────────────────────────────────────────────────
  wide <- read_wides(files = paste0("./data/raw/", fname))
  wide <- wide[as.numeric(wide$horas) <= TIME_MAX_H, ]
  for (or in 3:ncol(wide)) wide[, or] <- as.numeric(wide[, or])

  baseline <- as.matrix(wide[1, 3:ncol(wide)])
  wide[, 3:ncol(wide)] <- sweep(as.matrix(wide[, 3:ncol(wide)]), 2, baseline)
  dat_m <- wide[, 3:ncol(wide)]; dat_m[dat_m < 0] <- 0
  wide[, 3:ncol(wide)] <- dat_m

  long  <- melt(wide, id.vars = c("file","horas"))
  colnames(long)[3] <- "Well"
  sub_design <- my_design[my_design$Well %in% long$Well, ]
  merged     <- merge_dfs(long, sub_design)
  merged$horas <- as.numeric(merged$horas)
  merged$value <- as.numeric(merged$value)

  ## ── AUC by well ─────────────────────────────────────────────────────────
  ## gcplyr::auc() with pROC::auc()
  auc_sum <- summarize(group_by(merged, Chemical, Well),
                        auc = gcplyr::auc(x=horas, and=value),
                        .groups = "drop")
  auc_sum$well_pos <- as.integer(
    sapply(strsplit(as.character(auc_sum$Well), "\\."), `[`, 3))

  ctrl_rows <- auc_sum[auc_sum$Chemical == "CONTROL", ]
  cmpd_rows <- auc_sum[auc_sum$Chemical != "CONTROL", ]

  neg_aucs   <- ctrl_rows$auc
  mean_ctrl  <- mean(neg_aucs, na.rm=TRUE)
  med_ctrl   <- median(neg_aucs, na.rm=TRUE)
  mad_ctrl   <- mad(neg_aucs, na.rm=TRUE)

  dtw_row <- controls_dtw[controls_dtw$Placa == run_id, ]
  if (nrow(dtw_row) == 1) {
    cent_pos  <- dtw_row$centroid_pos[1]
    cent_auc  <- auc_sum[auc_sum$well_pos == cent_pos, "auc"]
    cent_val  <- if (length(cent_auc) == 1 && !is.na(cent_auc)) cent_auc else mean_ctrl
  } else {
    cent_val <- mean_ctrl   # fallback
  }

  ## ── B-score for compounds ──────────────────────────────────────────────
  cmpd_bs <- bscore_plate(as.data.frame(cmpd_rows))

  norm_df <- data.frame(
    Chemical    = cmpd_rows$Chemical,
    Well        = cmpd_rows$Well,
    auc_raw     = cmpd_rows$auc,
    auc_bs      = cmpd_bs$auc,
    run         = run_id,
    plate       = plate_id,
    stringsAsFactors = FALSE
  )

  norm_df$raw          <- norm_df$auc_raw
  norm_df$poc          <- norm_df$auc_raw / mean_ctrl
  norm_df$rzscore      <- (norm_df$auc_raw - med_ctrl) / max(mad_ctrl, 1e-6)
  norm_df$dtw_centroid <- norm_df$auc_raw / max(cent_val, 1e-6)
  norm_df$median       <- norm_df$auc_raw / max(med_ctrl, 1e-6)
  norm_df$bscore_med   <- norm_df$auc_bs  / max(med_ctrl, 1e-6)

  ## normalised (for SD_ctrl_norm)
  ## auc_bs = auc_raw controls (B-score no is applied on wells control)
  ctrl_norm_df <- data.frame(
    Chemical     = ctrl_rows$Chemical,
    Well         = ctrl_rows$Well,
    auc_raw      = ctrl_rows$auc,
    auc_bs       = ctrl_rows$auc,
    run          = run_id,
    plate        = plate_id,
    raw          = ctrl_rows$auc,
    poc          = ctrl_rows$auc / mean_ctrl,
    rzscore      = (ctrl_rows$auc - med_ctrl) / max(mad_ctrl, 1e-6),
    dtw_centroid = ctrl_rows$auc / max(cent_val, 1e-6),
    median       = ctrl_rows$auc / max(med_ctrl, 1e-6),
    bscore_med   = ctrl_rows$auc / max(med_ctrl, 1e-6),
    stringsAsFactors = FALSE
  )

  all_norms <- dplyr::bind_rows(
    all_norms,
    dplyr::bind_rows(
      cbind(norm_df,      type = "compound"),
      cbind(ctrl_norm_df, type = "control")
    )
  )
  if (i %% 10 == 0) cat(sprintf("  %d/%d procthatdos\n", i, length(files)))
}

cat("\n")

strategies <- c("raw","poc","rzscore","dtw_centroid","median","bscore_med")

## 1. SD of normalised controls (by plate-run, luego median)
ctrl_long <- all_norms[all_norms$type == "control", ]
ctrl_long <- pivot_longer(ctrl_long,
                           cols = all_of(strategies),
                           names_to = "strategy",
                           values_to = "norm_val")

sd_ctrl <- ctrl_long %>%
  group_by(strategy, run) %>%
  summarise(sd_run = sd(norm_val, na.rm=TRUE), .groups="drop") %>%
  group_by(strategy) %>%
  summarise(sd_ctrl_norm_median = median(sd_run, na.rm=TRUE),
            sd_ctrl_norm_mean   = mean(sd_run, na.rm=TRUE),
            .groups="drop")

cmpd_long <- all_norms[all_norms$type == "compound", ]

compute_icc <- function(df, strat) {
  df_s <- df[, c("Chemical","run","plate", strat)]
  colnames(df_s)[4] <- "val"
  wide_icc <- df_s %>%
    group_by(Chemical, plate) %>%
    mutate(rep = row_number()) %>%
    filter(rep <= 3) %>%
    ungroup() %>%
    pivot_wider(id_cols = c("Chemical","plate"),
                names_from  = "rep",
                values_from = "val",
                names_prefix = "rep")
  mat <- as.matrix(wide_icc[, grep("^rep", colnames(wide_icc))])
  mat <- mat[complete.cases(mat), ]
  if (nrow(mat) < 5) return(NA_real_)
  tryCatch(irr::icc(mat, model="twoway", type="agreement")$value,
           error=function(and) NA_real_)
}

icc_vals <- sapply(strategies, function(s) compute_icc(cmpd_long, s))
icc_df   <- data.frame(strategy = strategies,
                        icc = round(icc_vals, 4),
                        stringsAsFactors = FALSE)


compute_auroc <- function(df, strat) {
  df_s <- df[, c("Chemical", "plate", "run", strat)]
  colnames(df_s)[4] <- "val"
  df_mean <- df_s %>%
    group_by(Chemical, plate) %>%
    summarise(mean_val = mean(val, na.rm=TRUE), .groups="drop")
  ## Score: deviation of the null
  if (strat == "rzscore") {
    df_mean$score <- -df_mean$mean_val  
  } else {
    df_mean$score <- 1 - df_mean$mean_val  
  }
  df_mean$is_hit <- df_mean$Chemical %in% core_hits
  if (sum(df_mean$is_hit) < 2) return(NA_real_)
  tryCatch({
    roc_obj <- pROC::roc(df_mean$is_hit, df_mean$score,
                          quiet=TRUE, direction="<")
    as.numeric(pROC::auc(roc_obj))
  }, error=function(and) NA_real_)
}

auroc_vals <- sapply(strategies, function(s) compute_auroc(cmpd_long, s))
auroc_df   <- data.frame(strategy = strategies,
                           auroc = round(auroc_vals, 4),
                           stringsAsFactors = FALSE)

## 4. Z'-factor (for plate-runs with hits)
cmpd_long_ext <- all_norms[all_norms$type == "compound", ]
ctrl_long_ext <- all_norms[all_norms$type == "control", ]

compute_zprime_strat <- function(s) {
  zvals <- c()
  for (pr in unique(cmpd_long_ext$run)) {
    cmpd_pr <- cmpd_long_ext[cmpd_long_ext$run == pr, ]
    ctrl_pr <- ctrl_long_ext[ctrl_long_ext$run == pr, ]
    hits_pr <- cmpd_pr[cmpd_pr$Chemical %in% core_hits, s]
    neg_pr  <- ctrl_pr[, s]
    if (length(hits_pr) >= 2)
      zvals <- c(zvals, z_prime_fn(hits_pr, neg_pr))
  }
  median(zvals, na.rm=TRUE)
}

zprime_vals <- sapply(strategies, compute_zprime_strat)
zprime_df   <- data.frame(strategy = strategies,
                            zprime_median = round(zprime_vals, 4),
                            stringsAsFactors = FALSE)

## ── Tabla summary of the benchmark ──────────────────────────────────────────────
benchmark_tbl <- sd_ctrl %>%
  left_join(icc_df,    by="strategy") %>%
  left_join(auroc_df,  by="strategy") %>%
  left_join(zprime_df, by="strategy") %>%
  select(strategy, sd_ctrl_norm_median, zprime_median, icc, auroc)

strategy_order <- c("raw","poc","rzscore","dtw_centroid","median","bscore_med")
strategy_labels <- c("(1) Raw AUC",
                      "(2) Percent-of-control",
                      "(3) Robust Z-score",
                      "(4) DTW centroid (baseline)",
                      "(5) Median of controls",
                      "(6) Median + B-score (pipeline)")

benchmark_tbl$strategy <- factor(benchmark_tbl$strategy,
                                   levels = strategy_order,
                                   labels = strategy_labels)
benchmark_tbl <- benchmark_tbl[order(benchmark_tbl$strategy), ]

cat("\n=== Benchmark normalizaciones ===\n")
print(benchmark_tbl, row.names=FALSE)

write.table(benchmark_tbl,
            file.path(OUT_DIR, "Table_benchmark_norm.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)


p_sd <- ggplot(benchmark_tbl,
                aes(x=strategy, and=sd_ctrl_norm_median,
                    fill=strategy == "(6) Median + B-score (pipeline)")) +
  geom_col() +
  scale_fill_manual(values=c("FALSE"="steelblue","TRUE"="darkgreen"),
                    guide="none") +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=35, hjust=1, size=11)) +
  labs(title="SD of normalised controls (lower = better)",
       subtitle="Median across plate-runs",
       x=NULL, and="Median SD of normalised ctrl AUC")

## Panel B: Z'-factor (more alto = mejor)
p_zp <- ggplot(benchmark_tbl[!is.na(benchmark_tbl$zprime_median), ],
                aes(x=strategy, and=zprime_median,
                    fill=strategy == "(6) Median + B-score (pipeline)")) +
  geom_col() +
  geom_hline(yintercept=0.5, linetype="dashed", color="red", linewidth=0.8) +
  scale_fill_manual(values=c("FALSE"="steelblue","TRUE"="darkgreen"),
                    guide="none") +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=35, hjust=1, size=11)) +
  labs(title="Median Z'-factor (higher = better)",
       subtitle="Only plate-runs containing core hits; red dashed = 0.5",
       x=NULL, and="Median Z'-factor")

p_icc <- ggplot(benchmark_tbl[!is.na(benchmark_tbl$icc), ],
                 aes(x=strategy, and=icc,
                     fill=strategy == "(6) Median + B-score (pipeline)")) +
  geom_col() +
  scale_fill_manual(values=c("FALSE"="steelblue","TRUE"="darkgreen"),
                    guide="none") +
  coord_cartesian(ylim=c(0,1)) +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=35, hjust=1, size=11)) +
  labs(title="ICC between triplicates (higher = better)",
       x=NULL, and="ICC (two-way, agreement)")

p_auroc <- ggplot(benchmark_tbl[!is.na(benchmark_tbl$auroc), ],
                   aes(x=strategy, and=auroc,
                       fill=strategy == "(6) Median + B-score (pipeline)")) +
  geom_col() +
  geom_hline(yintercept=0.5, linetype="dashed", color="gray", linewidth=0.6) +
  scale_fill_manual(values=c("FALSE"="steelblue","TRUE"="darkgreen"),
                    guide="none") +
  coord_cartesian(ylim=c(0,1)) +
  theme_bw(base_size=13) +
  theme(axis.text.x = element_text(angle=35, hjust=1, size=11)) +
  labs(title="AUROC for core hits recovery (higher = better)",
       subtitle="Core hits (ANOVA ∩ LMM, FDR<0.05) as gold standard",
       x=NULL, and="AUROC")


cmpd_long_plot <- pivot_longer(cmpd_long[cmpd_long$type=="compound" |
                                           is.null(cmpd_long$type), ],
                                cols=all_of(strategies),
                                names_to="strategy",
                                values_to="norm_val")
cmpd_long_plot <- cmpd_long_plot[!is.na(cmpd_long_plot$norm_val), ]
cmpd_long_plot$strategy <- factor(cmpd_long_plot$strategy,
                                    levels=strategy_order,
                                    labels=strategy_labels)
cmpd_long_plot$norm_val[cmpd_long_plot$strategy=="(3) Robust Z-score" &
                           abs(cmpd_long_plot$norm_val) > 5] <- NA

p_dist <- ggplot(cmpd_long_plot, aes(x=norm_val, fill=strategy)) +
  geom_histogram(bins=60, color="white", linewidth=0.1) +
  facet_wrap(~strategy, scales="free", ncol=2) +
  theme_bw(base_size=13) +
  theme(legend.position="none",
        strip.text = element_text(size=11)) +
  labs(title="Distribution of normalised AUC per strategy",
       x="Normalised AUC", and="Count")

## -- Draft figure (review layout, original dimensions) ----------------
pdf(file.path(OUT_DIR, "Fig_benchmark_norm_draft.pdf"), width = 14, height = 18)
grid.arrange(p_sd, p_zp, p_icc, p_auroc, p_dist,
             layout_matrix = rbind(c(1,2), c(3,4), c(5,5)))
dev.off()

## -- Publication figure (PLOS format + TIFF definitive) --------------------
top_row <- cowplot::plot_grid(p_sd, p_zp, ncol = 2,
                               labels = c("A","B"), label_size = 18, label_fontface = "bold")
mid_row <- cowplot::plot_grid(p_icc, p_auroc, ncol = 2,
                               labels = c("C","D"), label_size = 18, label_fontface = "bold")
bot_row <- cowplot::plot_grid(p_dist, labels = "E", label_size = 18, label_fontface = "bold")
bench_pub <- cowplot::plot_grid(top_row, mid_row, bot_row, ncol = 1, rel_heights = c(1,1,1.2))

pdf(file.path(OUT_DIR, "Fig_benchmark_norm.pdf"), width = 6.83, height = 10)
print(bench_pub)
dev.off()
file.copy(file.path(OUT_DIR, "Table_benchmark_norm.txt"),
          file.path(DEFS_DIR, "S2_Table_NormBenchmark.txt"), overwrite = TRUE)

## ── Guardar also data long-format for uso posterior ────────────────────
saveRDS(cmpd_long, file.path(OUT_DIR, "benchmark_cmpd_norms.rds"))

cat("\nOutputs saved en:", OUT_DIR, "\n")
cat("  Table_benchmark_norm.txt\n")
cat("  Fig_benchmark_norm.pdf  (new Fig 5 of the manuscrito)\n")
cat("  benchmark_cmpd_norms.rds\n")
cat("\n=== 02_benchmark_norm.R completedo ===\n")
