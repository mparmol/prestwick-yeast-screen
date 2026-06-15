
library(gcplyr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(gplots)
library(gridExtra)

out_dir <- "."

PLATES_EXCLUDE <- c(1, 18)          # outlier plates with control CV > 45%
TIME_MAX_H     <- 60                # truncate at 60 h
CTRL_CV_WARN   <- 20                # warn if control CV exceeds this %

plate_num_from_file <- function(fname) {
  as.integer(sub("plate_(\\d+)_.*", "\\1", fname))
}


well_to_rowcol <- function(well_pos) {
  well_pos <- as.integer(well_pos)
  row_idx  <- ceiling(well_pos / 12)   # 1=A, 2=B, ..., 8=H
  col_idx  <- ((well_pos - 1) %% 12) + 1  # 1-12
  data.frame(row_idx = row_idx, col_idx = col_idx)
}


bscore_plate <- function(auc_df) {
  rc <- well_to_rowcol(auc_df$well_pos)
  auc_df$row_idx <- rc$row_idx
  auc_df$col_idx <- rc$col_idx

  cmpd <- auc_df[auc_df$col_idx %in% 2:11, ]

  if (nrow(cmpd) < 4) return(auc_df)  # too few samples

  mat <- matrix(NA, nrow = 8, ncol = 10,
                dimnames = list(1:8, 2:11))
  for (k in seq_len(nrow(cmpd))) {
    r <- cmpd$row_idx[k]
    c <- cmpd$col_idx[k] - 1  # offset: column 2 → index 1
    if (r >= 1 && r <= 8 && c >= 1 && c <= 10)
      mat[r, c] <- cmpd$auc[k]
  }

  mp <- tryCatch(medpolish(mat, na.rm = TRUE, trace.iter = FALSE),
                 error = function(and) NULL)

  if (is.null(mp)) return(auc_df)

  for (k in seq_len(nrow(cmpd))) {
    r <- cmpd$row_idx[k]
    c <- cmpd$col_idx[k] - 1
    if (r >= 1 && r <= 8 && c >= 1 && c <= 10) {
      idx <- which(auc_df$well_pos == cmpd$well_pos[k])

      row_eff <- if (!is.na(mp$row[r])) mp$row[r] else 0
      col_eff <- if (!is.na(mp$col[c])) mp$col[c] else 0
      auc_df$auc[idx] <- auc_df$auc[idx] - row_eff - col_eff
    }
  }
  auc_df
}

my_design <- read.csv(".", sep = "\t", encoding = "latin1")
files_all  <- list.files("./data/raw/")

files <- files_all[!sapply(files_all, function(f) plate_num_from_file(f) %in% PLATES_EXCLUDE)]
cat(sprintf("excluded: %s\n", paste(PLATES_EXCLUDE, collapse = ", ")))
cat(sprintf("procthatr: %d de %d\n", length(files), length(files_all)))
cat(sprintf(": %s\n\n",
            paste(sort(unique(sapply(files, plate_num_from_file))), collapse = ", ")))

tab_res_final   <- NULL
tab_res_controls <- NULL
diag_controls    <- NULL

pdf(file.path(out_dir, "o.pdf"), width = 28, height = 12)

for (i in seq_along(files)) {

  fname    <- files[i]
  plate_id <- plate_num_from_file(fname)
  run_id   <- sub("\\.csv$", "", fname)  # "plate_2_rep_1"

  imported_widedata <- read_wides(files = paste0("./data/raw/", fname))

  imported_widedata <- imported_widedata[as.numeric(imported_widedata$horas) <= TIME_MAX_H, ]

  for (or in 3:ncol(imported_widedata))
    imported_widedata[, or] <- as.numeric(imported_widedata[, or])

  baseline <- as.matrix(imported_widedata[1, 3:ncol(imported_widedata)])
  imported_widedata <- cbind(
    imported_widedata[, 1:2],
    sweep(as.matrix(imported_widedata[, 3:ncol(imported_widedata)]), 2, baseline)
  )

  dat_mat <- imported_widedata[, 3:ncol(imported_widedata)]
  dat_mat[dat_mat < 0] <- 0
  imported_widedata[, 3:ncol(imported_widedata)] <- dat_mat

  imported_widedata_resh <- melt(imported_widedata, id.vars = c("file", ""))
  colnames(imported_widedata_resh)[3] <- "Well"

  sub_my_design <- my_design[my_design$Well %in% imported_widedata_resh$Well, ]
  ex_dat_mrg    <- merge_dfs(imported_widedata_resh, sub_my_design)
  ex_dat_mrg$horas <- as.numeric(ex_dat_mrg$horas)
  ex_dat_mrg$value <- as.numeric(ex_dat_mrg$value)

  ex_dat_mrg$Well <- factor(ex_dat_mrg$Well, levels = sub_my_design$Well)

  p <- ggplot(ex_dat_mrg, aes(x = horas, and = value)) +
    geom_line() +
    facet_wrap(~Well + Chemical, nrow = 8, ncol = 12) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylab("OD600") + xlab("Time (h)") +
    ggtitle(run_id)
  print(p)

  ex_dat_mrg_sum <- summarize(
    group_by(ex_dat_mrg, Chemical, Well),
    auc = auc(x = horas, and = value)
  )

  ctrl_aucs   <- ex_dat_mrg_sum[ex_dat_mrg_sum$Chemical == "CONTROL", ]
  median_ctrl <- median(ctrl_aucs$auc, na.rm = TRUE)
  mean_ctrl   <- mean(ctrl_aucs$auc, na.rm = TRUE)
  sd_ctrl     <- sd(ctrl_aucs$auc, na.rm = TRUE)
  cv_ctrl     <- 100 * sd_ctrl / mean_ctrl

  diag_controls <- rbind(diag_controls, data.frame(
    run          = run_id,
    plate        = plate_id,
    n_ctrl       = nrow(ctrl_aucs),
    mean_ctrl    = round(mean_ctrl, 3),
    median_ctrl  = round(median_ctrl, 3),
    sd_ctrl      = round(sd_ctrl, 3),
    cv_ctrl_pct  = round(cv_ctrl, 2),
    flag         = ifelse(cv_ctrl > CTRL_CV_WARN, "HIGH_CV", "OK")
  ))

  tab_res_controls <- rbind(tab_res_controls,
                             cbind(ctrl_aucs, run = run_id, plate = plate_id,
                                   median_ctrl = median_ctrl))

  cmpd_aucs <- ex_dat_mrg_sum[ex_dat_mrg_sum$Chemical != "CONTROL", ]
  cmpd_aucs$well_pos <- as.integer(
    sapply(strsplit(as.character(cmpd_aucs$Well), "\\."), `[`, 3)
  )

  cmpd_aucs <- bscore_plate(cmpd_aucs)

  cmpd_aucs$AUC_relative   <- cmpd_aucs$auc / median_ctrl
  cmpd_aucs$AUC_difference <- ifelse(cmpd_aucs$AUC_relative >= 1, "Higher", "Lower")

  cmpd_aucs$run    <- run_id
  cmpd_aucs$plate  <- plate_id

  p2 <- ggplot(cmpd_aucs,
               aes(reorder(Chemical, -AUC_relative), AUC_relative,
                   fill = AUC_difference)) +
    geom_col() +
    theme_bw() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = c("Higher" = "orchid3", "Lower" = "skyblue3")) +
    scale_y_continuous(breaks = seq(0, 2, 0.1)) +
    ggtitle(run_id) + xlab(NULL)
  print(p2)

  tab_res_final <- rbind(tab_res_final, cmpd_aucs)
}

dev.off()

write.table(tab_res_final,
            file.path(out_dir, "f.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(tab_res_controls,
            file.path(out_dir, "c.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(diag_controls,
            file.path(out_dir, "crobust.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== controls (all the plate-runs) ===\n")
print(diag_controls)
cat(sprintf("\nCowith CV > %d%%: %d de %d\n",
            CTRL_CV_WARN,
            sum(diag_controls$flag == "HIGH_CV"),
            nrow(diag_controls)))

pdf(file.path(out_dir, "Boxplotso.pdf"), width = 20, height = 8)

for (pl in sort(unique(tab_res_final$plate))) {
  sub_tab <- tab_res_final[tab_res_final$plate == pl, ]
  p <- ggplot(sub_tab, aes(x = reorder(Chemical, -AUC_relative), and = AUC_relative)) +
    geom_boxplot() +
    theme_bw() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_continuous(breaks = seq(0, 2, 0.25)) +
    ggtitle(paste0("Pl_", pl, "  |  normalisation median + B-score")) +
    xlab(NULL)
  print(p)
}

dev.off()

cat("\n→ ANOVA per compound...\n")

res_significancia <- NULL

pdf(file.path(out_dir, "Boxplots_.pdf"))

for (compound in unique(tab_res_final$Chemical)) {
  sub_res  <- tab_res_final[tab_res_final$Chemical == compound, ]
  pl_id    <- sub_res$plate[1]

  sub_ctrl <- tab_res_controls[tab_res_controls$plate == pl_id, ]
  if (nrow(sub_ctrl) == 0 || nrow(sub_res) == 0) next

  dataset_analysis <- rbind(
    data.frame(Chemical = sub_res$Chemical, Well = sub_res$Well,
               auc = sub_res$auc),
    data.frame(Chemical = "CONTROL", Well = sub_ctrl$Well,
               auc = sub_ctrl$auc)
  )

  model          <- lm(auc ~ Chemical, data = dataset_analysis)
  anova_res       <- anova(model)
  pval            <- anova_res$`Pr(>F)`[1]

  res_significancia <- rbind(res_significancia, data.frame(
    Chemical     = compound,
    ANOVA        = pval,
    Plate        = pl_id,
    AUC_relative = mean(sub_res$AUC_relative, na.rm = TRUE),
    n_reps       = nrow(sub_res)
  ))

  p <- ggplot(dataset_analysis, aes(Chemical, auc)) +
    geom_boxplot() +
    ggtitle(compound) +
    theme_bw()
  print(p)
}

dev.off()

res_significancia$pcat <- cut(res_significancia$ANOVA,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 0.5, Inf),
                              labels = c("<0.001", "<0.01", "<0.05", "<0.1", "<0.5", "ns"))

res_significancia$FDR_BH <- p.adjust(res_significancia$ANOVA, method = "BH")

res_significancia <- res_significancia %>%
  group_by(Plate) %>%
  mutate(FDR_BH_plate = p.adjust(ANOVA, method = "BH")) %>%
  ungroup()

write.table(res_significancia,
            file.path(out_dir, "statistics_robust.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n")
cat("============================================================\n")
cat("Robust pipeline\n")
cat("============================================================\n")
cat(sprintf(":            %d\n", length(files)))
cat(sprintf("excluded:              %s\n", paste(PLATES_EXCLUDE, collapse = ", ")))
cat(sprintf("Computhiss:             %d\n", length(unique(tab_res_final$Chemical))))
cat(sprintf("AUC_rel media / median / DE:  %.3f / %.3f / %.3f\n",
            mean(res_significancia$AUC_relative, na.rm = TRUE),
            median(res_significancia$AUC_relative, na.rm = TRUE),
            sd(res_significancia$AUC_relative, na.rm = TRUE)))
cat(sprintf("Raw p < 0.05:                  %d\n", sum(res_significancia$ANOVA < 0.05, na.rm = TRUE)))
cat(sprintf("FDR (BH) < 0.05:               %d\n", sum(res_significancia$FDR_BH < 0.05, na.rm = TRUE)))
cat(sprintf("FDR (BH) < 0.10:               %d\n", sum(res_significancia$FDR_BH < 0.10, na.rm = TRUE)))
cat("============================================================\n")
cat("============================================================\n")
