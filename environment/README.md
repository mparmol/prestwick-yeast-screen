# Environment files

These files describe the R session and package versions under which
the analyses in this repository were executed. They are required by
PLOS Computational Biology's data-availability policy.

## How to (re)generate them

Open R at the repository root and run:

```r
# 1. Session info
writeLines(
  capture.output(sessionInfo()),
  con = "environment/sessionInfo.txt"
)

# 2. Loaded-package table (after running the full pipeline once)
pkgs <- subset(
  as.data.frame(installed.packages()[, c("Package", "Version")],
                stringsAsFactors = FALSE),
  Package %in% c(
    "gcscreen", "gcplyr", "lme4", "lmerTest", "brms", "rstan",
    "refund", "growthcurver", "mgcv", "irr", "pROC",
    "dplyr", "tidyr", "ggplot2", "cowplot", "gridExtra",
    "reshape2", "zoo"
  )
)
write.csv(pkgs, "environment/packages.csv", row.names = FALSE)
```

The current placeholder files are committed so that the directory
is not empty; please overwrite them with the output of the commands
above before the final submission.
