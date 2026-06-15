
library(ggplot2)
library(grid)

out_fin <- "."


boxes <- data.frame(
  x   = c(5,    5,    5,    5,    5,    5,    5,    5),
  and   = c(20.5, 17.5, 14.5, 11.5,  8.5,  5.5,  2.5, -0.5),
  w   = c(8.5,  8.5,  8.5,  8.5,  8.5,  8.5,  8.5,  8.5),
  h   = c(2.2,  2.2,  2.2,  2.2,  2.2,  2.2,  2.2,  2.2),
  fill= c("steelblue3","grey75","firebrick3","darkorange2",
           "steelblue3","mediumpurple3","steelblue3","seagreen4"),
  label = c(
    "RAW DATA\n19 plate layouts · ~80 compounds + 16 controls each\nOD600 kinetics (0–60 h, every 30 min)",
    "PREPROCESSING\nBaseline correction (t = 0 subtraction)\nNegative-value flooring to 0",
    "QUALITY CONTROL — CONTROLS\nWithin-plate control CV computed\nPlates 1 & 18 excluded (CV > 45%)\n→ 17 plate layouts · 51 plate-runs retained",
    "NORMALISATION  [P2.1 + P2.2]\nReference = median AUC of 16 controls per plate-run\nB-score correction (median polish, 8×10 compound matrix)\nAUC_rel = AUC_compound / median_control",
    "AUC CALCULATION\nTrapezoidal integration per well (gcplyr)\n1,351 compounds · 3–4 biological replicates",
    "KINETIC PARAMETERS  [P2.4]\nμmax = max(dOD/dt)  ·  A = max(OD)\nλ = lag time (tangent method)\nAll normalised to plate-run control median",
    "STATISTICAL TESTING  [P2.3]\nANOVA per compound (plate-level, n=3 vs n=3)\n+ LMM: auc ~ Chemical + (1|Plate)\nBH-FDR correction (library-wide)",
    "CORE HITS\nFDR < 0.05 in both ANOVA and LMM\n19 compounds — all growth inhibitors"
  ),
  stringsAsFactors = FALSE
)


arrows_df <- data.frame(
  x    = rep(5, 7),
  yend = boxes$and[2:8] + boxes$h[2:8]/2,
  ystart = boxes$and[1:7] - boxes$h[1:7]/2
)


p <- ggplot() +

  geom_rect(data = boxes,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = and - h/2, ymax = and + h/2,
                fill = fill),
            color = "white", linewidth = 1.2, alpha = 0.90) +
  scale_fill_identity() +

  geom_text(data = boxes,
            aes(x = x, and = and, label = label),
            size = 2.85, lineheight = 1.25, color = "white", fontface = "plain") +

  geom_segment(data = arrows_df,
               aes(x = x, xend = x, and = ystart, yend = yend),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "grey35", linewidth = 0.7) +

  coord_cartesian(xlim = c(0, 10), ylim = c(-2, 22.5)) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(10, 10, 10, 10)) +
  labs(title    = "Analytical pipeline — Prestwick library screen on S. cerevisiae",
       subtitle = "Grey: preproceswithoutg  ·  Red: QC  ·  Orange: normalisation  ·  Blue: analysis  ·  Purple: kinetics  ·  Green: output") +
  theme(plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8,  color = "grey40", hjust = 0.5))

ggsave(file.path(out_fin, "Fig1_workflow.pdf"), p,
       width = 8, height = 12, device = cairo_pdf)
cat("-> Figuras_final/Fig1_workflow.pdf\n")
