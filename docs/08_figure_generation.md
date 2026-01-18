# 08_figure_generation.R

## Purpose

Assembles analysis results into main and supplementary figures.

## Input

| File | Source |
|------|--------|
| `results/sample_metadata.csv` | Script 01 |
| `results/deg/deg_summary.csv` | Script 03 |
| `results/deg/axis_enrichment.csv` | Script 03 |
| `results/axis_scoring/axis_activation_stats.csv` | Script 04 |
| `results/temporal/gene_temporal_classification.csv` | Script 05 |
| `results/preservation/module_preservation_summary.csv` | Script 06 |
| `results/ha_analysis/gsea_all_results.csv` | Script 07 |
| `results/expr_gene_level.rds` | Script 02 |
| `results/qc/pca_coordinates.csv` | Script 01 |

## Output

### Main Figure
| File | Description |
|------|-------------|
| `figures/main/Figure_Main.pdf` | 6-panel main figure |
| `figures/main/Figure_Main.png` | PNG version (300 dpi) |

### Supplementary Figures
| File | Description |
|------|-------------|
| `figures/supplementary/SuppFig_S1_AllAxes.pdf` | All ECM axes across timepoints |
| `figures/supplementary/SuppFig_S2_GSEA_All.pdf` | Complete GSEA heatmap |
| `figures/supplementary/SuppFig_S3_PCA.pdf` | Sample QC PCA |
| `figures/supplementary/SuppFig_S4_DEG_Trajectory.pdf` | DEG counts over time |

## Main Figure Panels

| Panel | Content | Data Source |
|-------|---------|-------------|
| A | Pathway enrichment ranking | gsea_all_results.csv |
| B | HA-DAMP signaling NES over time | gsea_all_results.csv |
| C | HA gene expression trajectories | expr_gene_level.rds |
| D | Literature signature validation | gsea_all_results.csv |
| E | Temporal gene classification | gene_temporal_classification.csv |
| F | Module preservation | module_preservation_summary.csv |

## Figure Specifications

```r
# Main figure dimensions
width = 15
height = 11

# Panel layout
top_row <- pA + pB + pC + plot_layout(widths = c(1.2, 0.9, 1))
bottom_row <- pD + pE + pF + plot_layout(widths = c(1.1, 0.8, 1.1))
main_figure <- top_row / bottom_row
```

## Theme

Uses custom `theme_publication()` from `R/theme_publication.R`:
- Base font size: 11pt
- White background
- Minimal gridlines
- Color palettes defined for axes, conditions, timepoints

## Dependencies

```r
library(tidyverse)
library(patchwork)
library(ggplot2)
source("analysis/R/theme_publication.R")
```

