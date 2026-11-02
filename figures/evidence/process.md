# Analysis process

1. Define an independent normal population: mean 100, standard deviation 20.
2. Set R seed 640008 and draw 80 synthetic spending observations.
3. Generate 2,000 resampled means and a 95% percentile bootstrap interval.
4. Repeat fresh sampling and bootstrap interval construction 1,000 times to evaluate coverage.
5. Generate 50,000 samples under H0: population mean 95, treating SD 20 as known.
6. Compute observed Z, two-sided rejection boundaries and a plus-one Monte Carlo p-value.
7. Save numerical outputs, figures and R session information programmatically.
8. Generate English narrative from actual results; provide a standalone seeded browser lab.
9. Place website article at blog/posts/post8/index.qmd, with its display assets alongside it.

## Replication

From the repository root in RStudio run:

```r
source("code/blogpost8_all_process.R")
```

No external dataset or additional R package is required. The browser companion uses a different random generator, so its illustrative numbers differ from R.

## Observed validation

The R analysis completed successfully. The browser lab rendered and updating sample size from 80 to 40 changed the computed results. Website article metadata and all three image paths were checked locally. The website has not yet been rendered with Quarto or published from this session.
