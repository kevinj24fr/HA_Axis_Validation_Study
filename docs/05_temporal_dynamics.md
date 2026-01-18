# 05_temporal_dynamics.R

## Purpose

Classifies differentially expressed genes into temporal patterns based on when they are significant across the time course.

## Input

| File | Source |
|------|--------|
| `results/deg/deg_results_list.rds` | Script 03 |
| `results/sample_metadata.csv` | Script 01 |
| `results/reference/ecm_axes_rat.rds` | Script 02 |

## Output

### Data Files
| File | Description |
|------|-------------|
| `results/temporal/gene_temporal_classification.csv` | Gene-level pattern assignments |
| `results/temporal/axis_temporal_patterns.csv` | Axis-level pattern counts |
| `results/temporal/axis_pattern_matrix.csv` | Matrix of axis × pattern counts |

### Figures
| File | Description |
|------|-------------|
| `figures/05_gene_temporal_patterns.pdf` | Bar plot of pattern distribution |
| `figures/05_axis_temporal_heatmap.pdf` | Heatmap of axis × pattern |

## Methods

### Classification Criteria

Genes must be significant (FDR < 0.05, |logFC| ≥ 1) at ≥2 timepoints to be classified.

| Pattern | Definition |
|---------|------------|
| Resolving | Significant at Day 0-14, NOT at Day 126 |
| Persistent | Significant at early (Day 0-14) AND late (Day 126) |
| Late-emerging | NOT significant Day 0-14, significant Day 28 or 126 |

### Code

```r
gene_classification <- deg_results %>%
  group_by(gene) %>%
  summarize(
    n_sig = sum(is_significant),
    early_sig = any(is_significant & day <= 14),
    late_sig = any(is_significant & day >= 126)
  ) %>%
  filter(n_sig >= 2) %>%
  mutate(
    pattern = case_when(
      early_sig & !late_sig ~ "Resolving",
      early_sig & late_sig ~ "Persistent",
      !early_sig & late_sig ~ "Late-emerging",
      TRUE ~ "Other"
    )
  )
```

## Output Columns

### gene_temporal_classification.csv
| Column | Description |
|--------|-------------|
| gene | Gene symbol |
| n_sig_timepoints | Number of significant timepoints |
| early_significant | TRUE if DEG at Day 0-14 |
| late_significant | TRUE if DEG at Day 126 |
| pattern | Classification result |
| max_logFC | Maximum absolute logFC |
| peak_timepoint | Timepoint with max |logFC| |

## Dependencies

```r
library(tidyverse)
library(ggplot2)
library(pheatmap)
```
