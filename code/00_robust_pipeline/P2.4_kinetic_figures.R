library(ggplot2)
library(dplyr)
library(tidyr)

out_dir <- "."
lmm <- read.delim(file.path(out_dir,"statistics_LMM.txt"), encoding="latin1")
kin <- read.delim(file.path(out_dir,"Kinetic_forms_summary.txt"), encoding="latin1")

hits <- lmm$Chemical[!is.na(lmm$FDR_BH) & lmm$FDR_BH < 0.05]
ph   <- kin[kin$Chemical %in% hits, ]

## Aggregate per compound (puede aparecer multiples plates)
ph <- ph %>%
  group_by(Chemical) %>%
  summarise(
    mu_max_mean = mean(mu_max_mean, na.rm = TRUE),
    lambda_mean = mean(lambda_mean, na.rm = TRUE),
    A_mean      = mean(A_mean, na.rm = TRUE),
    .groups = "drop"
  )

cat("Hits with kinetic data:", nrow(ph), "\n")
ph$Chemical <- factor(ph$Chemical, levels = ph$Chemical[order(ph$A_mean)])

# Scatter mu_max vs A, coloreado by delta-lambda
p1 <- ggplot(ph, aes(mu_max_mean, A_mean, color = lambda_mean, label = Chemical)) +
  geom_point(size = 3) +
  geom_text(size = 2.8, vjust = -0.7, check_overlap = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_gradient2(low = "steelblue", mid = "grey80", high = "firebrick",
                        midpoint = 0, name = "Delta-lambda (h)") +
  labs(title    = "Kinetic profile of FDR<0.05 hits (LMM)",
       subtitle = "Axes: normalized to plate control median",
       x = "Relative mu_max",
       and = "Relative A (carrying capacity)") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "Fig_kinetic_hits_scatter.pdf"), p1, width = 8, height = 7)
cat("-> Guardado: Fig_kinetic_hits_scatter.pdf\n")

# Heatmap 3 formetros
ph2 <- ph[, c("Chemical","mu_max_mean","lambda_mean","A_mean")]
hlong <- pivot_longer(ph2, -Chemical, names_to = "Param", values_to = "Value")
hlong$Param <- factor(hlong$Param,
                       levels = c("mu_max_mean","lambda_mean","A_mean"),
                       labels = c("mu_max_rel","Delta-lambda (h)","A_rel"))
p2 <- ggplot(hlong, aes(Param, Chemical, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "steelblue3", mid = "grey95", high = "firebrick3",
                       midpoint = 1, name = "Value") +
  labs(title = "Kinetic formeters heatmap - FDR<0.05 hits (LMM)", x = NULL, and = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "Fig_kinetic_heatmap.pdf"), p2,
       width = 6, height = max(5, 0.4 * nrow(ph)))
cat("-> Guardado: Fig_kinetic_heatmap.pdf\n")

# Comforcion ANOVA robusto vs LMM
comp <- read.delim(file.path(out_dir, "ANOVA_vs_LMM_comparison.txt"), encoding = "latin1")
cat("\n=== Concordancia ANOVA robusto vs. LMM (FDR<0.05) ===\n")
cat("Solo ANOVA:", sum(comp$FDR_BH_ANOVA < 0.05 & (is.na(comp$FDR_BH_LMM) | comp$FDR_BH_LMM >= 0.05)), "\n")
cat("Solo LMM:  ", sum(comp$FDR_BH_LMM   < 0.05 & (is.na(comp$FDR_BH_ANOVA) | comp$FDR_BH_ANOVA >= 0.05)), "\n")
cat("En ambos:     ", sum(comp$FDR_BH_ANOVA < 0.05 & comp$FDR_BH_LMM < 0.05, na.rm = TRUE), "\n")

hits_lmm   <- comp$Chemical[!is.na(comp$FDR_BH_LMM)   & comp$FDR_BH_LMM   < 0.05]
hits_anova <- comp$Chemical[!is.na(comp$FDR_BH_ANOVA) & comp$FDR_BH_ANOVA < 0.05]
core_hits  <- intersect(hits_lmm, hits_anova)
cat("\n===  ===\n")
print(data.frame(Compound = core_hits,
                 row.names = NULL))
