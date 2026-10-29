# Final-production figures

This folder contains the exact October 29, 2026 production scripts and aggregate input CSVs for Figures 3–5 of **When negative control exposures fail to detect bias**.

These scripts use the same simulation results and plotting definitions as the original scripts in `figures/`, with the submitted dimensions, Arial typography, and 12-point panel titles. No simulations are run by these scripts.

## Requirements

- R with macOS Quartz support and Arial available. These exact production scripts are macOS-specific.
- ggplot2. Reproduction was checked with R 4.5.2 and ggplot2 4.0.1.
- For other operating systems, the original `figures/MakeFigure_D0.R`, `MakeFigure_D1.R`, and `MakeFigure_D2.R` scripts use the standard R PDF device and reproduce the plotted results with different dimensions and typography; see the main README for their input-file workflow.

## Reproduce Figures 3–5

From a terminal in the repository root (or `CodeForPublication/` after extracting Appendix S4):

```sh
mkdir -p figures/scratch/production
Rscript --vanilla figures/production/MakeFigure_D0_final.R figures/production figures/scratch/production
Rscript --vanilla figures/production/MakeFigure_D1_final.R figures/production figures/scratch/production
Rscript --vanilla figures/production/MakeFigure_D2_final.R figures/production figures/scratch/production
```

The first argument is the directory containing the three input CSVs. The second is the output directory; create it before invoking the scripts. The output names are `Figure_3_final.pdf`, `Figure_4_final.pdf`, and `Figure_5_final.pdf`. PDF metadata, such as creation dates, can differ between runs even when the figures render identically.

| Paper figure | Script | Included input | PDF size (inches) |
| --- | --- | --- | --- |
| Figure 3: Design 0 | `MakeFigure_D0_final.R` | `D0_agg.csv` | 7 × 6 |
| Figure 4: Design 1 | `MakeFigure_D1_final.R` | `D1_agg.csv` | 7 × 4.2 |
| Figure 5: Design 2A | `MakeFigure_D2_final.R` | `D2_agg.csv` | 7 × 4.2 |

The included CSVs are byte-for-byte copies of these saved runs in the GitHub repository:

| Included input | Original result file |
| --- | --- |
| `D0_agg.csv` | `results/D0/D0_norm_n1000_it1000_seed314_20260605_174455_agg.csv` |
| `D1_agg.csv` | `results/D1/D1_norm_n1000_it1000_seed314_20260605_181841_agg.csv` |
| `D2_agg.csv` | `results/D2/D2_norm_n1000_it1000_seed314_20260605_182834_agg.csv` |

These fixed inputs reproduce the submitted results. Running a new simulation does not update them. The original plotting scripts' `run_prefix` setting selects newly generated results, as described in the [main README](../../README.md).
