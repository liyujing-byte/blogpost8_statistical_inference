# Simulation, Bootstrap Confidence Intervals, and Hypothesis Testing

## Project overview

A base R simulation explains bootstrap percentile intervals, repeated-sample coverage, and two-sided hypothesis testing. All spending observations are synthetic, independent normal draws (mean 100, SD 20), not real household data.

## Repository structure

- code/blogpost8_all_process.R: complete seeded analysis, figures, article and validation generation
- data/raw/: synthetic sample and source/model documentation
- data/processed/: observation-indexed analysis table
- data/results/: bootstrap draws, coverage intervals, null draws, numerical summaries and R session information
- figures/: three PNG figures, interactive browser lab and standalone article preview
- figures/post/: generated English article in Markdown and Quarto formats
- figures/evidence/: process notes, rubric map, run log, numerical checks, output checksums and browser validation record

This follows the directory organization of Blog Post 7.

## Replication

Open blogpost8_statistical_inference.Rproj in RStudio. From the project root run:

```r
source("code/blogpost8_all_process.R")
```

Or run `Rscript code/blogpost8_all_process.R`. Base R only; no package installation or dataset download is needed. R seed: 640008. There are 80 observations, 2,000 bootstrap resamples per interval, 1,000 independent coverage experiments and 50,000 null simulations. The coverage experiment may take a minute.

Open figures/index.html for the local article preview, or figures/interactive.html for the standalone interactive lab. GitHub file preview does not execute JavaScript. Browser simulations use a separate seeded generator and differ numerically from R.

## Findings and assumptions

The 95% bootstrap interval is [94.28, 101.81]; repeated coverage is 94.9% (Monte Carlo SE 0.70 percentage points). The two-sided test of population mean 95 gives Z = 1.370 and simulated p = 0.1718, so it does not reject at 5%. The test assumes known SD 20. Percentile bootstrap coverage is approximate; synthetic normal independent observations do not establish validity for biased, dependent or strongly skewed samples.

## Article and publication

Website source: https://github.com/liyujing-byte/website/tree/main/blog/posts/post8

Blog URL: https://liyujing-byte.github.io/website/blog/posts/post8/

Analysis repository: https://github.com/liyujing-byte/blogpost8_statistical_inference

The website source contains the article, display figures and interactive lab. This repository contains all analysis data, code and process/evidence. The blog URL becomes live after the website deployment completes.
