
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)

out_dir <- "."

## ── Load data ────────────────────────────────────────────────────────────
fin  <- read.delim(file.path(out_dir, "final_results_robust.txt"),
                   encoding = "latin1", stringsAsFactors = FALSE)
ctrl <- read.delim(file.path(out_dir, "control_results_robust.txt"),
                   encoding = "latin1", stringsAsFactors = FALSE)

## Extract replicate from the Well ID (2ª component: A, B, C, D)
fin$Rep  <- sapply(strsplit(as.character(fin$Well), "\\."), `[`, 2)
ctrl$Rep <- sapply(strsplit(as.character(ctrl$Well), "\\."), `[`, 2)

ctrl_long <- data.frame(
  Chemical = "CONTROL",
  Well     = ctrl$Well,
  auc      = ctrl$auc,
  plate    = ctrl$plate,
  run      = ctrl$run,
  Rep      = ctrl$Rep
)

fin_long <- data.frame(
  Chemical = fin$Chemical,
  Well     = fin$Well,
  auc      = fin$auc,
  plate    = fin$plate,
  run      = fin$run,
  Rep      = fin$Rep
)

long <- rbind(fin_long, ctrl_long)
long$Plate_f <- factor(long$plate)
long$Rep_f   <- factor(long$Rep)

cat(sprintf("Obs",
            nrow(long), length(unique(fin$Chemical)), length(unique(fin$plate))))

cat("\n→ \n")

m_null <- lmer(auc ~ 1 + (1|Plate_f) + (1|Rep_f), data = long, REML = TRUE)
vc <- as.data.frame(VarCorr(m_null))
vc$pct_var <- round(100 * vc$vcov / sum(vc$vcov), 2)

cat("\n=== Components variance (model nulo) ===\n")
print(vc[, c("grp", "vcov", "pct_var")])

write.table(vc[, c("grp", "vcov", "sdcor", "pct_var")],
            file.path(out_dir, "mixed_model_variance_components_robust.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n→  ...\n")
cat("   ...\n")

compounds <- unique(fin$Chemical)
mix_results <- NULL

for (cmpd in compounds) {
  sub <- long[long$Chemical %in% c(cmpd, "CONTROL"), ]
  sub$Chemical_f <- factor(sub$Chemical, levels = c("CONTROL", cmpd))

  n_plates <- length(unique(sub$Plate_f))

  if (n_plates < 2) {
    ## Fallback: ANOVA simple (identical to pipeline original)
    model <- lm(auc ~ Chemical_f, data = sub)
    an     <- anova(model)
    pval   <- an$`Pr(>F)`[1]
    method <- "ANOVA"
    est    <- coef(model)[2]
    se     <- summary(model)$coefficients[2, 2]
  } else {
    fit <- tryCatch(
      lmer(auc ~ Chemical_f + (1|Plate_f), data = sub, REML = FALSE),
      error = function(and) NULL,
      warning = function(w) suppressWarnings(
        lmer(auc ~ Chemical_f + (1|Plate_f), data = sub, REML = FALSE)
      )
    )
    if (is.null(fit)) {
      ## Fallback
      model <- lm(auc ~ Chemical_f, data = sub)
      an     <- anova(model)
      pval   <- an$`Pr(>F)`[1]
      method <- "ANOVA_fallback"
      est    <- coef(model)[2]
      se     <- summary(model)$coefficients[2, 2]
    } else {
      cf     <- coef(summary(fit))
      ## El effect of the compound is the 2ª row (Chemical_f<cmpd>)
      if (nrow(cf) >= 2) {
        est    <- cf[2, "Estimate"]
        se     <- cf[2, "Std. Error"]
        pval   <- cf[2, "Pr(>|t|)"]
      } else {
        est  <- NA; se <- NA; pval <- NA
      }
      method <- "LMM"
    }
  }

  auc_rel_mean <- mean(fin[fin$Chemical == cmpd, "AUC_relative"], na.rm = TRUE)

  mix_results <- rbind(mix_results, data.frame(
    Chemical     = cmpd,
    Plate        = unique(sub$plate[sub$Chemical == cmpd])[1],
    AUC_relative = auc_rel_mean,
    Estimate     = est,
    SE           = se,
    p_value      = pval,
    method       = method,
    stringsAsFactors = FALSE
  ))
}

## BH-FDR
mix_results$FDR_BH <- p.adjust(mix_results$p_value, method = "BH")

mix_results <- mix_results %>%
  group_by(Plate) %>%
  mutate(FDR_BH_plate = p.adjust(p_value, method = "BH")) %>%
  ungroup()

mix_results$Direction <- ifelse(mix_results$AUC_relative >= 1, "Higher", "Lower")

mix_results <- mix_results[order(mix_results$FDR_BH), ]

write.table(mix_results,
            file.path(out_dir, "statistics_LMM.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Distribution FDR (LMM, library-wide) ===\n")
cat(sprintf("FDR < 0.01:  %d\n", sum(mix_results$FDR_BH < 0.01, na.rm = TRUE)))
cat(sprintf("FDR < 0.05:  %d\n", sum(mix_results$FDR_BH < 0.05, na.rm = TRUE)))
cat(sprintf("FDR < 0.10:  %d\n", sum(mix_results$FDR_BH < 0.10, na.rm = TRUE)))
cat(sprintf("FDR < 0.20:  %d\n", sum(mix_results$FDR_BH < 0.20, na.rm = TRUE)))

cat(sprintf("\nLMM used: %d/%d compounds\n",
            sum(mix_results$method == "LMM"),
            nrow(mix_results)))

robust_stat <- tryCatch(
  read.delim(file.path(out_dir, "statistics_robust.txt"), encoding = "latin1"),
  error = function(and) NULL
)

if (!is.null(robust_stat)) {
  comp <- merge(
    mix_results[, c("Chemical", "p_value", "FDR_BH")],
    robust_stat[, c("Chemical", "ANOVA", "FDR_BH")],
    by = "Chemical", suffixes = c("_LMM", "_ANOVA")
  )

  ## Hits by method
  hits_lmm   <- comp$Chemical[comp$FDR_BH_LMM < 0.05]
  hits_anova <- comp$Chemical[comp$FDR_BH_ANOVA < 0.05]

  cat("\n=== Hits FDR<0.05: ANOVA robusto vs. LMM ===\n")
  cat(sprintf("Solo ANOVA:  %d\n", sum(!hits_anova %in% hits_lmm)))
  cat(sprintf("Solo LMM:    %d\n", sum(!hits_lmm   %in% hits_anova)))
  cat(sprintf("En ambos:       %d\n", sum(hits_lmm %in% hits_anova)))

  write.table(comp,
              file.path(out_dir, "ANOVA_vs_LMM_comparison.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

mix_results$logP  <- -log10(mix_results$p_value)
mix_results$logFC <- log2(mix_results$AUC_relative)
mix_results$sig   <- ifelse(mix_results$FDR_BH < 0.05, "FDR<0.05",
                     ifelse(mix_results$FDR_BH < 0.10, "FDR<0.10", "ns"))
mix_results$sig   <- factor(mix_results$sig, levels = c("FDR<0.05", "FDR<0.10", "ns"))

lab_data <- mix_results[!is.na(mix_results$FDR_BH) & mix_results$FDR_BH < 0.05, ]

p_volc <- ggplot(mix_results[!is.na(mix_results$logP), ],
                 aes(logFC, logP, color = sig)) +
  geom_point(aes(size = sig), alpha = 0.6) +
  scale_size_manual(values = c("FDR<0.05" = 2.5, "FDR<0.10" = 1.8, "ns" = 1.2)) +
  scale_color_manual(values = c("FDR<0.05" = "firebrick",
                                "FDR<0.10" = "darkorange",
                                "ns" = "grey70")) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_text(data = lab_data,
            aes(label = Chemical), size = 2.5, vjust = -0.6,
            check_overlap = TRUE) +
  labs(title    = "Volcano plot — Mixed model (LMM, plate as random effect)",
       subtitle = sprintf("n = %d compounds",
                          nrow(mix_results)),
       x = "log2(AUC_rel)",
       and = "-log10(LMM p-value)",
       color = "", size = "") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "Fig_volcano_LMM.pdf"), p_volc, width = 9, height = 7)

## ── 5. Top hits LMM ──────────────────────────────────────────────────────────
top30 <- mix_results[!is.na(mix_results$FDR_BH), ]
top30 <- top30[order(top30$FDR_BH), ][1:min(30, nrow(top30)), ]
top30$Chemical <- factor(top30$Chemical, levels = top30$Chemical[order(top30$AUC_relative)])

p_top <- ggplot(top30, aes(Chemical, AUC_relative, fill = FDR_BH)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_gradientn(colors = c("firebrick", "darkorange", "gold", "grey80"),
                       values = scales::rescale(c(0, 0.05, 0.10, 0.20)),
                       name = "FDR (BH)") +
  coord_flip() +
  labs(title    = "Top 30 hits — Mixed model (LMM)",
       subtitle = "FDR_BH",
       x = NULL, and = "AUC_rel (media)") +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, "Fig_top30_LMM.pdf"), p_top, width = 8, height = 9)

cat("\n→ \n")
cat("  statistics_LMM.txt\n")
cat("  mixed_model_variance_components_robust.txt\n")
cat("  ANOVA_vs_LMM_comparison.txt\n")
cat("  Fig_volcano_LMM.pdf\n")
cat("  Fig_top30_LMM.pdf\n")
