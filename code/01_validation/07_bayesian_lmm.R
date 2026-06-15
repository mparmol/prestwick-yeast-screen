## ============================================================
## 07_bayesian_lmm.R  --  C.3: Bayesian hierarchical model
## Gaussian model with plate random effects (brms/Stan)
## for the 18 core hits; comparison with frequentist LMM.
##
## Run from: d:/Projects/Proyecto_curvas/
## Outputs -> ./data/intermediate/validation/
## ============================================================

library(brms)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(tidyr)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, forllel::detectCores()))

out_dir  <- "./data/intermediate/validation"
defs_dir <- file.path(out_dir, "definitive")
dir.create(defs_dir, showWarnings = FALSE, recursive = TRUE)

## -- Core hits ----------------------------------------------------------------
CORE_HITS <- c(
  "Alexidine dihydrochloride", "Cycloheximide", "Thimerosal",
  "Thonzonium bromide", "Merbromin", "Lanowithazole",
  "Amiodarone hydrochloride", "Amorolfine hydrochloride",
  "Butowithazole nitrate", "Hexachlorophene", "Doxorubicin hydrochloride",
  "Nitroxoline", "Dequalinium dichloride", "Clioquinol",
  "Clomiphene citrate (Z,E)", "Ciclopirox ethanolamine",
  "Oxiwithazole Nitrate", "Tolcapone"
)

## -- Load data ----------------------------------------------------------------
cat("Loading data...\n")

res_pipeline <- read.delim(
  "./data/intermediate/primary_pipeline/final_results_robust.txt",
  encoding = "latin1", stringsAsFactors = FALSE
)
res_ctrl <- read.delim(
  "./data/intermediate/primary_pipeline/control_results_robust.txt",
  encoding = "latin1", stringsAsFactors = FALSE
)
lmm_res <- read.delim(
  "./data/intermediate/primary_pipeline/statistics_LMM.txt",
  encoding = "latin1", stringsAsFactors = FALSE
)

cat(sprintf("Compounds: %d  |  Runs: %d\n",
            length(unique(res_pipeline$Chemical)),
            length(unique(res_pipeline$run))))

## -- Subset: core hits + controls from the same plates -----------------------
hits_data  <- res_pipeline[res_pipeline$Chemical %in% CORE_HITS, ]
hit_plates <- unique(hits_data$plate)
ctrl_data  <- res_ctrl[res_ctrl$plate %in% hit_plates, ]

long_hits <- data.frame(
  Chemical = hits_data$Chemical,
  auc      = hits_data$auc,
  plate    = as.character(hits_data$plate),
  stringsAsFactors = FALSE
)

long_ctrl <- data.frame(
  Chemical = "CONTROL",
  auc      = ctrl_data$auc,
  plate    = as.character(ctrl_data$plate),
  stringsAsFactors = FALSE
)

model_data <- rbind(long_hits, long_ctrl)
model_data$Chemical <- factor(model_data$Chemical,
                               levels = c("CONTROL", CORE_HITS))
model_data$plate    <- factor(model_data$plate)

cat(sprintf("\nModel data:\n"))
cat(sprintf("  Total rows:     %d\n", nrow(model_data)))
cat(sprintf("  Compounds:      %d (+ CONTROL)\n", length(CORE_HITS)))
cat(sprintf("  Plates:         %d\n", nlevels(model_data$plate)))
cat(sprintf("  Control AUC:    mean = %.2f  SD = %.2f\n",
            mean(long_ctrl$auc), sd(long_ctrl$auc)))

## -- Priors -------------------------------------------------------------------
priors <- c(
  prior(normal(0, 20), class = b),           # compound effects
  prior(exponential(0.1), class = sd),        # between-plate SD
  prior(exponential(0.1), class = sigma)      # residual SD
)

## -- Model fitting ------------------------------------------------------------
cat("\nFitting Bayesian model (brms)...\n")
cat("  Family: Gaussian  |  Formula: auc ~ Chemical + (1|plate)\n")
cat("  Chains: 4  |  Iterations: 2000  |  Warmup: 1000\n\n")

set.seed(2026)
fit <- brm(
  formula = auc ~ Chemical + (1 | plate),
  data    = model_data,
  family  = gaussian(),
  prior   = priors,
  chains  = 4L,
  iter    = 2000L,
  warmup  = 1000L,
  seed    = 2026L,
  refresh = 200L,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

cat("\nModel fitted.\n")
print(summary(fit)$fixed[1:5, ])

## -- Posterior extraction -----------------------------------------------------
cat("\nExtracting posteriors...\n")

post_draws <- as.data.frame(fit)

# All fixed-effect compound formeters
all_chem_forms <- names(post_draws)[grepl("^b_Chemical", names(post_draws))]

# brms strips ALL non-alphanumeric characters from factor level names in
# as.data.frame() output. sub() gives the encoded name; decoded maps it back.
bayes_results <- lapply(all_chem_forms, function(par) {
  draws     <- post_draws[[par]]
  cmpd_name <- sub("^b_Chemical", "", par)   # original compound name
  p_inhib   <- mean(draws < 0)               # P(AUC reduction vs CONTROL)

  data.frame(
    form             = par,
    Chemical          = cmpd_name,
    mean_effect       = mean(draws),
    sd_effect         = sd(draws),
    q2.5              = quantile(draws, 0.025),
    q97.5             = quantile(draws, 0.975),
    p_inhibition      = p_inhib,
    significant_bayes = p_inhib >= 0.975,    # equivalent to two-tailed alpha = 0.05
    stringsAsFactors  = FALSE
  )
})

bayes_df <- do.call(rbind, bayes_results)
bayes_df <- bayes_df[order(bayes_df$mean_effect), ]

cat(sprintf("\nHits with P(inhibition) >= 0.975: %d / %d\n",
            sum(bayes_df$significant_bayes), nrow(bayes_df)))

## -- Comparison with frequentist LMM ------------------------------------------
lmm_hits <- lmm_res[, c("Chemical", "p_value", "FDR_BH",
                          "AUC_relative", "method")]
names(lmm_hits)[names(lmm_hits) == "FDR_BH"] <- "FDR_lmm"

# brms strips ALL non-alphanumeric characters from factor level names.
# Rebuild original names via reverse lookup before merging.
decoded <- setNames(CORE_HITS, gsub("[^[:alnum:]]", "", CORE_HITS))
bayes_df$Chemical_orig <- decoded[bayes_df$Chemical]

comparison <- merge(
  bayes_df[, c("Chemical_orig", "mean_effect", "sd_effect",
               "q2.5", "q97.5", "p_inhibition", "significant_bayes")],
  lmm_hits,
  by.x  = "Chemical_orig",
  by.and  = "Chemical",
  all.x = TRUE
)
comparison <- comparison[order(comparison$mean_effect), ]
names(comparison)[1] <- "Chemical"

cat("\n=== Bayesian vs LMM comparison -- core hits ===\n")
print(comparison[, c("Chemical", "AUC_relative", "p_inhibition",
                      "FDR_lmm", "significant_bayes")])

## -- MCMC diagnostics ---------------------------------------------------------
cat("\n=== MCMC diagnostics ===\n")

fit_sum   <- summary(fit)
rhat_vals <- c(fit_sum$fixed[, "Rhat"], fit_sum$random$plate[, "Rhat"])
bulk_ess  <- c(fit_sum$fixed[, "Bulk_ESS"], fit_sum$random$plate[, "Bulk_ESS"])
total_samples <- 4L * (2000L - 1000L)   # chains * post-warmup iterations
ess_ratio <- bulk_ess / total_samples

n_divergences <- sum(nuts_forms(fit, pars = "divergent__")$Value)

cat(sprintf("Max Rhat:           %.4f  (threshold: < 1.01)\n", max(rhat_vals, na.rm = TRUE)))
cat(sprintf("Min ESS ratio:      %.4f  (threshold: > 0.1)\n",  min(ess_ratio, na.rm = TRUE)))
cat(sprintf("Divergences:        %d\n", n_divergences))

## -- Table 1: per-hit Bayesian results + LMM comparison ----------------------
write.table(comparison,
            file.path(out_dir, "Table_bayesian_hits.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## -- Table 2: global model metrics -------------------------------------------
summary_stats <- data.frame(
  statistic = c(
    "N hits fitted",
    "Hits P(inhibition) >= 0.975",
    "Hits 95% CI excludes zero (beta < 0)",
    "Max Rhat",
    "Min ESS ratio",
    "MCMC divergences",
    "Concordance with LMM (FDR < 0.05 in both)"
  ),
  value = c(
    nrow(bayes_df),
    sum(bayes_df$significant_bayes),
    sum(bayes_df$q97.5 < 0),
    round(max(rhat_vals, na.rm = TRUE), 4),
    round(min(ess_ratio, na.rm = TRUE), 4),
    n_divergences,
    sum(comparison$significant_bayes &
        !is.na(comparison$FDR_lmm) & comparison$FDR_lmm < 0.05,
        na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.table(summary_stats,
            file.path(out_dir, "Table_bayesian_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## -- Figure 1: Forest plot of 95% credible intervals -------------------------
plot_df <- comparison[!is.na(comparison$Chemical), ]
plot_df$Chemical <- factor(plot_df$Chemical,
                            levels = plot_df$Chemical[order(-plot_df$mean_effect)])

p_forest <- ggplot(plot_df,
                   aes(x = mean_effect, and = Chemical,
                       xmin = q2.5, xmax = q97.5,
                       color = significant_bayes)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(orientation = "and", width = 0.3, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "grey60"),
    labels = c("TRUE" = "P(inhibition) >= 0.975", "FALSE" = "Not significant"),
    name   = NULL
  ) +
  labs(
    title    = "95% credible intervals -- Bayesian hierarchical model",
    subtitle = sprintf("brms  |  Gaussian  |  AUC ~ compound + (1|plate)  |  %d chains x 1000 iter.",
                       4L),
    x        = "Compound effect (AUC difference vs CONTROL)",
    and        = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Fig_bayesian_forest_draft.pdf"),
       p_forest, width = 9, height = 7)
ggsave(file.path(out_dir, "Fig_bayesian_forest.pdf"),
       p_forest, width = 6.83, height = 7)

## -- Figure 2: Bayesian posterior mean vs empirical AUC_relative --------------
# p_inhibition = 1.000 for all 18 hits (no variance), so scatter of
# Bayesian mean effect vs AUC_relative is more informative for the paper.
comp2 <- comparison[!is.na(comparison$AUC_relative), ]
comp2$neg_log10_fdr <- -log10(pmax(comp2$FDR_lmm, 1e-10))

# Pearson correlation between Bayesian posterior mean and AUC_relative
r_val <- cor(comp2$mean_effect, comp2$AUC_relative, method = "pearson")

p_compare <- ggplot(comp2,
                    aes(x = AUC_relative, and = mean_effect,
                        label = Chemical)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue",
              fill = "steelblue", alpha = 0.15, linewidth = 0.8) +
  geom_point(color = "firebrick", size = 2.5) +
  geom_text(size = 2.8, hjust = -0.05, vjust = 0.5,
            check_overlap = TRUE) +
  annotate("text", x = max(comp2$AUC_relative) * 0.85,
           and = min(comp2$mean_effect) * 0.9,
           label = sprintf("r = %.3f", r_val), size = 4) +
  labs(
    title    = "Bayesian posterior mean vs empirical AUC",
    subtitle = "All 18 core hits: P(inhibition) = 1.000 | Concordance with LMM = 18/18",
    x        = "AUC relative (empirical; 1 = no effect)",
    and        = "Posterior mean compound effect (vs CONTROL)"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "Fig_bayesian_vs_lmm_draft.pdf"),
       p_compare, width = 7, height = 6)
ggsave(file.path(out_dir, "Fig_bayesian_vs_lmm.pdf"),
       p_compare, width = 6.83, height = 5.5)

## -- Figure 3: MCMC trace plots (convergence diagnostics) --------------------
# Rename brms internal formeter names to readable compound names.
# as.array(fit) → 3D array [iterations, chains, formeters] required by mcmc_trace.
draws_arr  <- as.array(fit)
par_names  <- dimnames(draws_arr)[[3]]
par_labels <- par_names
for (i in seq_along(par_names)) {
  if (grepl("^b_Chemical", par_names[i])) {
    enc <- sub("^b_Chemical", "", par_names[i])
    lbl <- decoded[enc]
    if (!is.na(lbl)) par_labels[i] <- lbl
  }
}
dimnames(draws_arr)[[3]] <- par_labels

# Select 4 representative hits (by encoded order) + residual SD
top4        <- all_chem_forms[1:min(4L, length(all_chem_forms))]
top4_labels <- decoded[sub("^b_Chemical", "", top4)]
top4_labels <- top4_labels[!is.na(top4_labels)]

trace_plot <- mcmc_trace(draws_arr, pars = c(top4_labels, "sigma"),
                         facet_args = list(nrow = 5, ncol = 1))
pdf(file.path(out_dir, "Fig_bayesian_trace_draft.pdf"), width = 10, height = 8)
print(trace_plot)
dev.off()
pdf(file.path(out_dir, "Fig_bayesian_trace.pdf"), width = 6.83, height = 8)
print(trace_plot)
dev.off()
## -- Figure 4: Posterior densities for all 18 hits ---------------------------
# Build plain matrix with readable compound names, ordered by mean effect
# (most inhibited compound at top of the plot).
all_labels <- decoded[sub("^b_Chemical", "", all_chem_forms)]
draw_mat   <- as.matrix(post_draws[, all_chem_forms])
colnames(draw_mat) <- all_labels
draw_mat   <- draw_mat[, order(colMeans(draw_mat))]  # ascending: most negative first

p_post <- mcmc_areas(draw_mat,
                     prob      = 0.95,
                     point_est = "mean") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title    = "Posterior distributions -- compound effect vs CONTROL",
       subtitle = "Shaded area: 95% CI  |  red line: beta = 0")
pdf(file.path(out_dir, "Fig_bayesian_posteriors_draft.pdf"), width = 10, height = 8)
print(p_post)
dev.off()
pdf(file.path(out_dir, "Fig_bayesian_posteriors.pdf"), width = 6.83, height = 8)
print(p_post)
dev.off()
## -- Summary ------------------------------------------------------------------
cat("\n")
cat("============================================================\n")
cat("SUMMARY C.3 -- Bayesian hierarchical model\n")
cat("============================================================\n")
cat(sprintf("Hits fitted:                       %d\n", nrow(bayes_df)))
cat(sprintf("P(inhibition) >= 0.975:            %d / %d\n",
            sum(bayes_df$significant_bayes), nrow(bayes_df)))
cat(sprintf("95%% CI excludes zero (beta < 0):   %d / %d\n",
            sum(bayes_df$q97.5 < 0), nrow(bayes_df)))
cat(sprintf("Max Rhat:                          %.4f\n",
            max(rhat_vals, na.rm = TRUE)))
cat(sprintf("Min ESS ratio:                     %.4f\n",
            min(ess_ratio, na.rm = TRUE)))
cat(sprintf("Concordance with LMM (18/18):      %d / %d hits\n",
            sum(comparison$significant_bayes &
                !is.na(comparison$FDR_lmm) & comparison$FDR_lmm < 0.05,
                na.rm = TRUE),
            nrow(comparison)))
file.copy(file.path(out_dir, "Table_bayesian_hits.txt"),
          file.path(defs_dir, "S8_Table_BayesianHits.txt"), overwrite = TRUE)
file.copy(file.path(out_dir, "Table_bayesian_summary.txt"),
          file.path(defs_dir, "S9_Table_BayesianSummary.txt"), overwrite = TRUE)
cat("============================================================\n")
cat("Outputs in ./data/intermediate/validation/:\n")
cat("  Table_bayesian_hits.txt\n")
cat("  Table_bayesian_summary.txt\n")
cat("  Fig_bayesian_forest.pdf\n")
cat("  Fig_bayesian_vs_lmm.pdf\n")
cat("  Fig_bayesian_trace.pdf\n")
cat("  Fig_bayesian_posteriors.pdf\n")
cat("  definitive/: S8_Table_BayesianHits.txt  S9_Table_BayesianSummary.txt\n")
cat("============================================================\n")
