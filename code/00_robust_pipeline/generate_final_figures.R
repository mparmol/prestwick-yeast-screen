
library(ggplot2)
library(dplyr)

out_rm  <- "."
out_ae  <- "."
out_fin <- "."


res <- read.delim(file.path(out_rm, "statistics_robust.txt"), encoding = "latin1")


res_clean <- res[res$Chemical != "Chicago sky blue 6B", ]

mean_auc <- mean(res_clean$AUC_relative, na.rm = TRUE)
med_auc  <- median(res_clean$AUC_relative, na.rm = TRUE)
n_inh    <- sum(res_clean$AUC_relative < 1, na.rm = TRUE)
n_enh    <- sum(res_clean$AUC_relative > 1, na.rm = TRUE)

p_hist <- ggplot(res_clean, aes(x = AUC_relative)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 1,        linetype = "dashed",  color = "black",     linewidth = 0.8) +
  geom_vline(xintercept = mean_auc, linetype = "dotted",  color = "firebrick", linewidth = 0.9) +
  geom_vline(xintercept = med_auc,  linetype = "dotted",  color = "darkgreen", linewidth = 0.9) +
  annotate("text", x = mean_auc - 0.015, and = Inf,
           label = sprintf("mean = %.3f", mean_auc),
           vjust = 2, hjust = 1, color = "firebrick", size = 3.5) +
  annotate("text", x = med_auc + 0.015, and = Inf,
           label = sprintf("median = %.3f", med_auc),
           vjust = 2, hjust = 0, color = "darkgreen", size = 3.5) +
  annotate("text", x = 0.2,  and = Inf, label = sprintf("Inhibitors\nn = %d", n_inh),
           vjust = 2, color = "steelblue4", size = 3.2) +
  annotate("text", x = 1.25, and = Inf, label = sprintf("Enhancers\nn = %d", n_enh),
           vjust = 2, color = "orchid4", size = 3.2) +
  scale_x_continuous(breaks = seq(0, 1.6, 0.1)) +
  labs(
    title    = "Distribution of relative AUC across 1,351 Prestwick compounds",
    subtitle = sprintf("n = %d | SD = %.3f | range %.3f–%.3f | Plates 1 & 18 excluded",
                       nrow(res_clean),
                       sd(res_clean$AUC_relative, na.rm = TRUE),
                       min(res_clean$AUC_relative, na.rm = TRUE),
                       max(res_clean$AUC_relative, na.rm = TRUE)),
    x = "AUC_rel  (compound / median control)",
    and = "Number of compounds"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(out_fin, "Fig2A_AUCrel_distribution.pdf"), p_hist, width = 9, height = 5)
cat("-> Fig2A_AUCrel_distribution.pdf\n")


file.copy(file.path(out_ae, "Fig3_plate_effect.pdf"),
          file.path(out_fin, "Fig3_plate_effect_diagnostic.pdf"),
          overwrite = TRUE)
cat("-> Fig3_plate_effect_diagnostic.pdf\n")

file.copy(file.path(out_rm, "Fig_volcano_LMM.pdf"),
          file.path(out_fin, "Fig2B_volcano_LMM.pdf"),
          overwrite = TRUE)
cat("-> Fig2B_volcano_LMM.pdf\n")

file.copy(file.path(out_rm, "Fig_top30_LMM.pdf"),
          file.path(out_fin, "Fig4_top30_hits_LMM.pdf"),
          overwrite = TRUE)
cat("-> Fig4_top30_hits_LMM.pdf\n")

file.copy(file.path(out_rm, "Fig_kinetic_heatmap.pdf"),
          file.path(out_fin, "Fig5_kinetic_heatmap_hits.pdf"),
          overwrite = TRUE)
cat("-> Fig5_kinetic_heatmap_hits.pdf\n")

file.copy(file.path(out_ae, "FigS1_control_CV.pdf"),
          file.path(out_fin, "FigS1_control_CV_per_plate.pdf"),
          overwrite = TRUE)
cat("-> FigS1_control_CV_per_plate.pdf\n")

file.copy(file.path(out_rm, "Fig_kinetic_distributions.pdf"),
          file.path(out_fin, "FigS2_kinetic_distributions.pdf"),
          overwrite = TRUE)
cat("-> FigS2_kinetic_distributions.pdf\n")

file.copy(file.path(out_rm, "Fig_kinetic_hits_scatter.pdf"),
          file.path(out_fin, "FigS3_kinetic_scatter_hits.pdf"),
          overwrite = TRUE)
cat("-> FigS3_kinetic_scatter_hits.pdf\n")

stat_all <- read.delim(".", encoding = "latin1")
fin_all  <- read.delim(".", encoding = "latin1")

stat_all$plate_num <- as.integer(sapply(strsplit(as.character(stat_all$Placa), "\\."), `[`, 1))
fin_all$plate_num  <- as.integer(sapply(strsplit(as.character(fin_all$Well), "\\."), `[`, 1))

excl_stat <- stat_all[stat_all$plate_num %in% c(1, 18),
                       c("Chemical","ANOVA","V3","Placa","AUC_relative","plate_num")]
excl_stat <- excl_stat[!duplicated(excl_stat$Chemical), ]
excl_stat <- excl_stat[order(excl_stat$ANOVA), ]

header_warning <- paste(
  "# TABLE S1 — Compounds from plates 1 and 18 (EXCLUDED from main analysis)",
  "# These compounds were screened but their results are flagged as LOW CONFIDENCE",
  "# Reason: within-plate control coefficient of variation > 45% in both plate runs",
  "# (Plate 1: CV = 48.7%; Plate 18: CV = 45.9%)",
  "# Statistics shown use the ORIGINAL DTW-centroid pipeline (not the robust pipeline).",
  "# These 161 compounds require independent re-screening to draw valid withclusions.",
  "# -------------------------------------------------------------------------",
  sep = "\n"
)
writeLines(header_warning,
           with = file.path(out_fin, "TableS1_plates1_18_excluded_WARNING.txt"))
write.table(excl_stat,
            file = file.path(out_fin, "TableS1_plates1_18_excluded_WARNING.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)

cat(sprintf("-> TableS1_plates1_18_excluded_WARNING.txt  (%d compounds)\n", nrow(excl_stat)))

cat("\n=== Figuras_final/ ===\n")
files_fin <- list.files(out_fin)
for (f in files_fin) cat(sprintf("  %s\n", f))
cat(sprintf("\nTotal: %d files\n", length(files_fin)))
