# Code for Figure Generation

> **Paper:** A phase II peri-operative study of pembrolizumab plus lenvatinib for mucosal melanoma

This repository contains the R code used to generate figures in the above manuscript. All visualization scripts are organized as modular, reusable code templates in the `code/` directory.

---

## Repository Structure

```text
├── code/                               # Modular R visualization scripts
│   ├── Heatmap.R                       # Heatmap (pheatmap)
│   ├── Box_Plot.R                      # Box plot with statistical comparison
│   ├── Paired_Line_Chart.R             # Paired line chart (Pre vs Post)
│   ├── Volcano_Plot.R                  # Volcano plot (DEG results)
│   ├── GSEA_Enrichment_Plot.R          # GSEA enrichment plot (GseaVis)
│   ├── Enrichment_Bar_Chart.R          # Enrichment bar chart (NES ranking)
│   ├── Kaplan-Meier_Survival_Curve.R   # K-M survival curve with risk table
│   ├── Tile_Plot.R                     # Tile plot (sample annotation)
│   ├── Stacked_Column_Chart.R          # Stacked column chart (TCR clonal ratio)
│   ├── Oncoprint.R                     # Oncoprint (mutation landscape)
│   ├── Error_Bound_Based_Sample_Size_Calculation.R
│   └── Hypothesis_Test-based_Sample_Size_Calculation.R
└── README.md                           # This file
```

---

## Code Modules

### 1. Heatmap (`code/Heatmap.R`)

- **Package:** `pheatmap`
- **Input:** TSV file with row names, a `Group` column for annotation, and additional annotation columns
- **Figure panels:** Fig 3a, Fig 4c
- **Statistical method:** Not applicable
- **Parameters to adjust:** `gaps_col`, `gaps_row`, annotation colors in `ann_colors`, color palette

### 2. Box Plot (`code/Box_Plot.R`)

- **Packages:** `tidyverse`, `reshape2`, `ggpubr`
- **Input:** TSV file with `Group` column and numeric columns
- **Figure panels:** Fig 3b, Fig 3j, Fig 5a/c/e, Fig 5b/c/d (box plots), Fig S3
- **Statistical method:** Configurable via the `method` argument in `stat_compare_means()`
  - Mann–Whitney U test (default): `method = "wilcox.test"` (equivalent to two-sided Mann–Whitney U)
  - Unpaired t-test: `method = "t.test"` — **required for Fig 3b and Fig 5e**
- **Parameters to adjust:** `comparisons` list (group pairs), `Group` factor levels, `method` argument

### 3. Paired Line Chart (`code/Paired_Line_Chart.R`)

- **Packages:** `ggplot2`, `ggpubr`
- **Input:** TSV file with `Timepoint`, `Value`, `Group`, `ID2` columns
- **Figure panels:** Fig 3c, Fig 4e, Fig 4g, Fig 5f, Fig 5g
- **Statistical method:** Two-sided Wilcoxon signed-rank test (`method = "wilcox.test", paired = T`)
- **Parameters to adjust:** Y-axis label (`expression`), comparison groups

### 4. Volcano Plot (`code/Volcano_Plot.R`)

- **Packages:** `tidyverse`, `openxlsx`, `ggrepel`, `ggtext`
- **Input:** XLSX file with `logFC` and `P.Value` columns
- **Figure panels:** Fig 4a
- **Statistical method:** Not applicable
- **Parameters to adjust:** Significance thresholds (`|logFC| >= 1`, `P.Value < 0.05`)

### 5. GSEA Enrichment Plot (`code/GSEA_Enrichment_Plot.R`)

- **Package:** `GseaVis`
- **Input:** Pre-loaded R object `GSEA_result` (generated from GSEA analysis)
- **Figure panels:** Fig 4b, Fig 4d
- **Statistical method:** Not applicable
- **Parameters to adjust:** `geneSetID` (pathway name), `subPlot` (plot layout)

### 6. Enrichment Bar Chart (`code/Enrichment_Bar_Chart.R`)

- **Package:** `ggplot2`
- **Input:** R data frame `gsea_go` with `pathway` and `NES` columns; `title` variable
- **Figure panels:** Fig 3d
- **Statistical method:** Not applicable
- **Parameters to adjust:** `title`, `gsea_go` data, NES color scale

### 7. K-M Survival Curve (`code/Kaplan-Meier_Survival_Curve.R`)

- **Packages:** `tidyverse`, `openxlsx`, `survival`, `survminer`
- **Input:** XLSX file with `RFS_Month`, `RFS_Event`, `Group` (high/low) columns
- **Figure panels:** Fig 3k, Fig 3m, Fig 3o, Fig 3q
- **Statistical method:** Log-rank test (computed via `survdiff`)
- **Parameters to adjust:** `xlim`, `break.x.by`, `break.y.by`, palette colors

### 8. Tile Plot (`code/Tile_Plot.R`)

- **Packages:** `tidyverse`, `reshape2`
- **Input:** TSV file with `ID` column and categorical annotation variables
- **Figure panels:** Fig 5b (combined with stacked column chart)
- **Statistical method:** Not applicable
- **Parameters to adjust:** Color mapping in `cols`, annotation variables

### 9. Stacked Column Chart (`code/Stacked_Column_Chart.R`)

- **Packages:** `tidyverse`
- **Input:** TSV file with `name` and `aaSeqCDR3` columns; separate `group.txt` for group annotation
- **Figure panels:** Fig 5b, Fig 5d
- **Statistical method:** Not applicable
- **Parameters to adjust:** `levels` in `factor()` (sample names), expansion threshold (`0.001`)

### 10. Oncoprint (`code/Oncoprint.R`)

- **Packages:** `ComplexHeatmap`, `GetoptLong`
- **Input:** 
  - Mutation data file (TSV): sample ID, gene name, variant classification
  - Phenotype file (TSV): sample ID with TMB and TNB annotations
  - Group file (TSV): sample ID with group labels (Pre/Post × Responder/Non-responder)
  - Sample file (optional): ordered sample list for column arrangement
  - Classification file (optional): gene classification for row splitting
- **Figure panels:** Fig 2
- **Statistical method:** Not applicable
- **Parameters to adjust:** `cut` (minimum sample count per gene), `width`, `height`, `colname` (show/hide column names), `atype`/`btype` (Cosmic/IntOGene tumor type annotations), mutation color palette (`mutcol`), gene order (`SortSmutCountPer`), group colors

---

## Environment

**R version:** 4.5.1

**Required R packages:**

| Package            | Version  | Usage                             |
| ------------------ | -------- | --------------------------------- |
| `tidyverse`        | 2.0.0    | Data manipulation and ggplot2     |
| ↳ `ggplot2`        | 4.0.2    | Data visualization                |
| ↳ `dplyr`          | 1.2.1    | Data manipulation                 |
| ↳ `tidyr`          | 1.3.2    | Data tidying                      |
| ↳ `readr`          | 2.2.0    | Reading rectangular data          |
| ↳ `purrr`          | 1.2.2    | Functional programming            |
| ↳ `tibble`         | 3.3.1    | Modern data frames                |
| ↳ `stringr`        | 1.6.0    | String manipulation               |
| ↳ `forcats`        | 1.0.1    | Factor handling                   |
| ↳ `lubridate`      | 1.9.5    | Date/time handling                |
| `ggpubr`           | 0.6.3    | Statistical annotations on plots  |
| `ggrepel`          | 0.9.8    | Label placement (volcano plot)    |
| `ggtext`           | 0.1.2    | Markdown rendering in axis labels |
| `pheatmap`         | 1.0.13   | Heatmap visualization             |
| `ComplexHeatmap`   | 2.24.1   | Oncoprint and complex heatmaps    |
| `GetoptLong`       | 1.2.5    | Command-line argument parsing     |
| `clusterProfiler`  | 4.18.4   | GSEA and functional enrichment    |
| `GseaVis`          | 0.0.5    | GSEA enrichment plots             |
| `survival`         | 3.8-6    | Survival analysis                 |
| `survminer`        | 0.5.2    | K-M survival curve plotting       |
| `openxlsx`         | 4.2.8.1  | Reading XLSX files                |
| `reshape2`         | 1.4.5    | Data reshaping (melt)             |

---

## Reproducing Figures

### Step-by-step workflow

1. **Install required R packages** (see list above).
2. **Prepare input data** — processed data are provided as separate sheets in the Supplementary Data file, each sheet corresponding to a figure panel.
3. **Load the appropriate code module** from `code/` and adjust parameters as needed:
   - For **box plots** and **line charts**, modify the `method` argument to match the statistical test required for each panel (see module descriptions above).
   - For the **oncoprint**, provide mutation data, phenotype annotations, and group assignments via command-line arguments.
   - For all modules, update input file paths and group/sample names as needed.
4. **Run the script** in R to generate the figure.
