
library(gcplyr)
library(dplyr)
library(ggplot2)
library(tidyr)

out_dir <- "."
PLATES_EXCLUDE <- c(1, 18)
TIME_MAX_H     <- 60
SMOOTH_WINDOW  <- 5   

plate_num_from_file <- function(fname) {
  as.integer(sub("plate_(\\d+)_.*", "\\1", fname))
}


deriv_num <- function(and, x) {
  n <- length(and)
  d <- numeric(n)
  d[1]   <- (and[2] - and[1]) / (x[2] - x[1])
  d[n]   <- (and[n] - and[n-1]) / (x[n] - x[n-1])
  for (i in 2:(n-1))
    d[i] <- (and[i+1] - and[i-1]) / (x[i+1] - x[i-1])
  d
}

moving_avg <- function(and, k = 5) {
  n   <- length(and)
  out <- and
  for (i in seq_len(n)) {
    lo       <- max(1, i - floor(k/2))
    hi       <- min(n, i + floor(k/2))
    out[i]   <- mean(and[lo:hi], na.rm = TRUE)
  }
  out
}

extract_forms <- function(x, and, threshold_pct = 0.1) {
  y_sm <- moving_avg(and, k = SMOOTH_WINDOW)
  y_sm[y_sm < 0] <- 0

  A <- max(y_sm, na.rm = TRUE)

  if (length(x) < 3 || all(y_sm == 0)) {
    return(data.frame(mu_max = 0, lambda = NA, A_cap = A))
  }

  deriv <- deriv_num(y_sm, x)
  deriv[is.nan(deriv) | is.infinite(deriv)] <- 0
  deriv[deriv < 0] <- 0  
  mu_max <- max(deriv, na.rm = TRUE)
  if (mu_max <= 0) return(data.frame(mu_max = 0, lambda = NA, A_cap = A))

  idx_mu <- which.max(deriv)
  t_mu   <- x[idx_mu]
  y_mu   <- y_sm[idx_mu]


  lambda <- t_mu - y_mu / mu_max
  if (lambda < 0) lambda <- 0
  if (lambda > t_mu) lambda <- NA

  data.frame(mu_max = mu_max, lambda = lambda, A_cap = A)
}

my_design <- read.csv(".", sep = "\t", encoding = "latin1")
files_all  <- list.files("./data/raw/")
files      <- files_all[!sapply(files_all, function(f) plate_num_from_file(f) %in% PLATES_EXCLUDE)]

cat(sprintf("Procthatndo %d plate-runs (plates %s excluded)...\n",
            length(files), paste(PLATES_EXCLUDE, collapse = ", ")))

all_forms <- NULL

for (i in seq_along(files)) {
  fname    <- files[i]
  plate_id <- plate_num_from_file(fname)
  run_id   <- sub("\\.csv$", "", fname)

  if (i %% 10 == 0) cat(sprintf("  Procthatndo %d/%d: %s\n", i, length(files), run_id))

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

  imported_widedata_resh <- reshape2::melt(imported_widedata, id.vars = c("file", ""))
  colnames(imported_widedata_resh)[3] <- "Well"
  sub_my_design <- my_design[my_design$Well %in% imported_widedata_resh$Well, ]
  ex_dat_mrg    <- merge_dfs(imported_widedata_resh, sub_my_design)
  ex_dat_mrg$horas <- as.numeric(ex_dat_mrg$horas)
  ex_dat_mrg$value <- as.numeric(ex_dat_mrg$value)

  wells <- unique(ex_dat_mrg$Well)

  for (w in wells) {
    sub_w    <- ex_dat_mrg[ex_dat_mrg$Well == w, ]
    sub_w    <- sub_w[order(sub_w$horas), ]
    chemical <- unique(sub_w$Chemical)
    if (length(chemical) != 1) next

    forms <- extract_forms(sub_w$horas, sub_w$value)

    all_forms <- rbind(all_forms, data.frame(
      run      = run_id,
      plate    = plate_id,
      Well     = w,
      Chemical = chemical,
      mu_max   = forms$mu_max,
      lambda   = forms$lambda,
      A_cap    = forms$A_cap,
      stringsAsFactors = FALSE
    ))
  }
}

cat(sprintf("\n: %d\n", nrow(all_forms)))

all_forms$Rep <- sapply(strsplit(as.character(all_forms$Well), "\\."), `[`, 2)

ctrl_medians <- all_forms %>%
  filter(Chemical == "CONTROL") %>%
  group_by(run) %>%
  summarise(
    ctrl_mu_max = median(mu_max, na.rm = TRUE),
    ctrl_lambda = median(lambda, na.rm = TRUE),
    ctrl_A_cap  = median(A_cap, na.rm = TRUE),
    .groups     = "drop"
  )

all_forms <- left_join(all_forms, ctrl_medians, by = "run")

all_forms$mu_max_rel <- all_forms$mu_max / all_forms$ctrl_mu_max
all_forms$lambda_dif <- all_forms$lambda - all_forms$ctrl_lambda  # difference absoluta (h)
all_forms$A_rel      <- all_forms$A_cap  / all_forms$ctrl_A_cap

forms_summary <- all_forms %>%
  filter(Chemical != "CONTROL") %>%
  group_by(Chemical, plate) %>%
  summarise(
    n_reps       = n(),
    mu_max_mean  = mean(mu_max_rel, na.rm = TRUE),
    mu_max_sd    = sd(mu_max_rel, na.rm = TRUE),
    lambda_mean  = mean(lambda_dif, na.rm = TRUE),
    lambda_sd    = sd(lambda_dif, na.rm = TRUE),
    A_mean       = mean(A_rel, na.rm = TRUE),
    A_sd         = sd(A_rel, na.rm = TRUE),
    .groups      = "drop"
  )

write.table(all_forms,
            file.path(out_dir, "Kinetic_forms_all_wells.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(forms_summary,
            file.path(out_dir, "Kinetic_forms_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("→ Kinetic_forms_all_wells.txt\n")
cat("→ Kinetic_forms_summary.txt\n")


lmm_stat <- tryCatch(
  read.delim(file.path(out_dir, "statistics_LMM.txt"), encoding = "latin1"),
  error = function(and) NULL
)

if (!is.null(lmm_stat)) {
  hits_fdr05 <- lmm_stat$Chemical[!is.na(lmm_stat$FDR_BH) & lmm_stat$FDR_BH < 0.05]

  forms_hits <- forms_summary[forms_summary$Chemical %in% hits_fdr05, ]

  if (nrow(forms_hits) > 0) {
    p_scatter <- ggplot(forms_hits,
                        aes(x = mu_max_mean, and = A_mean, color = lambda_mean,
                            label = Chemical)) +
      geom_point(size = 3) +
      geom_text(size = 2.8, vjust = -0.7, check_overlap = TRUE) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      scale_color_gradient2(low = "steelblue", mid = "grey80", high = "firebrick",
                            midpoint = 0, name = "Δλ (h)") +
      labs(title    = "Kinetic profile of FDR<0.05 hits",
           subtitle = "Axes normalized to plate control median",
           x        = "Relative μmax (max growth rate)",
           and        = "Relative A (carrying capacity)") +
      theme_bw(base_size = 12)

    ggsave(file.path(out_dir, "Fig_kinetic_hits_scatter.pdf"),
           p_scatter, width = 8, height = 7)
    cat("→ Guardado: Fig_kinetic_hits_scatter.pdf\n")

    hits_mat <- forms_hits %>%
      select(Chemical, mu_max_mean, lambda_mean, A_mean) %>%
      arrange(A_mean)

    hits_long <- tidyr::pivot_longer(hits_mat, -Chemical,
                                     names_to = "Param", values_to = "Value")
    hits_long$Param <- factor(hits_long$Param,
                               levels = c("mu_max_mean", "lambda_mean", "A_mean"),
                               labels = c("μmax_rel", "Δλ (h)", "A_rel"))

    p_heat <- ggplot(hits_long, aes(x = Param, and = Chemical, fill = Value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "steelblue3", mid = "grey95", high = "firebrick3",
                           midpoint = 1, name = "Value") +
      labs(title = "Kinetic formeters heatmap — FDR<0.05 hits",
           x = NULL, and = NULL) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    ggsave(file.path(out_dir, "Fig_kinetic_heatmap.pdf"),
           p_heat, width = 6, height = max(5, 0.4 * nrow(forms_hits)))
    cat("→ Fig_kinetic_heatmap.pdf\n")
  }
}

cmpd_forms <- all_forms[all_forms$Chemical != "CONTROL", ]

p_mu <- ggplot(cmpd_forms, aes(x = mu_max_rel)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "firebrick") +
  labs(title = "Distribution of relative maximum growth rate (μmax_rel)",
       x = "μmax_rel  (compound / control median)",
       and = "Count") +
  theme_bw(base_size = 12)

p_A <- ggplot(cmpd_forms, aes(x = A_rel)) +
  geom_histogram(binwidth = 0.05, fill = "mediumseagreen", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "firebrick") +
  labs(title = "Distribution of relative carrying capacity (A_rel)",
       x = "A_rel  (compound / control median)",
       and = "Count") +
  theme_bw(base_size = 12)

p_lam <- ggplot(cmpd_forms[!is.na(cmpd_forms$lambda_dif), ], aes(x = lambda_dif)) +
  geom_histogram(binwidth = 1, fill = "orchid3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick") +
  labs(title = "Distribution of lag time difference (Δλ)",
       x = "Δλ  (compound − control median, hours)",
       and = "Count") +
  theme_bw(base_size = 12)

pdf(file.path(out_dir, "Fig_kinetic_distributions.pdf"), width = 8, height = 10)
gridExtra::grid.arrange(p_mu, p_A, p_lam, ncol = 1)
dev.off()
cat("→ Fig_kinetic_distributions.pdf\n")

cat("\n=== Summary parameters cineticos ===\n")
cat(sprintf("Curvas compound procthatdas:  %d\n", nrow(cmpd_forms)))
cat(sprintf("μmax_rel — media / median:      %.3f / %.3f\n",
            mean(cmpd_forms$mu_max_rel, na.rm = TRUE),
            median(cmpd_forms$mu_max_rel, na.rm = TRUE)))
cat(sprintf("A_rel    — media / median:      %.3f / %.3f\n",
            mean(cmpd_forms$A_rel, na.rm = TRUE),
            median(cmpd_forms$A_rel, na.rm = TRUE)))
cat(sprintf("Δλ       — media / median (h):  %.2f / %.2f\n",
            mean(cmpd_forms$lambda_dif, na.rm = TRUE),
            median(cmpd_forms$lambda_dif, na.rm = TRUE)))
cat("===================================================\n")
