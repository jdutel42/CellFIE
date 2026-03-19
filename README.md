<div align="center">

# 🔬 CellFIE
### Cells & Features Interactive Explorer

**An interactive R Shiny application for single-cell RNA-seq data exploration and visualization**

[![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?style=flat-square&logo=r)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-1.8+-blue?style=flat-square)](https://shiny.posit.co/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18+-87B13F?style=flat-square)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](LICENSE)

</div>

---

## Overview

**CellFIE** (Cells & Features Interactive Explorer) is a modular R Shiny application designed to facilitate interactive exploration of single-cell RNA sequencing (scRNA-seq) datasets. It provides biologists with an intuitive, no-code interface to visualize gene expression, inspect cell populations, assess quality control metrics, and export publication-ready figures — all from a `SingleCellExperiment` or `Seurat` object.

CellFIE was built with non-bioinformaticians in mind: every plot is fully customizable (color scales, themes, labels, transparency) and exportable in **PNG**, **PDF**, and **SVG** formats at user-defined resolutions.

---

## Features

| Module | Description |
|---|---|
| 🏠 **Home** | Data loading, object summary (cells, genes, assays, reductions) |
| 🔬 **FeaturePlot** | Gene expression overlay on any low-dimensional embedding |
| 🎻 **ViolinPlot** | Expression distribution across cell groups, with optional boxplot and jitter |
| 🔵 **DotPlot** | Mean expression × percent expressed per group, with optional gene clustering |
| 🧬 **Cell Explorer** | Any metadata visualized on embedding, with centroid labels and cell composition |
| 🧪 **QC** | Quality control metrics: violin distributions and pairwise scatter with regression |
| 📋 **Metadata Table** | Filterable, exportable cell-level metadata (CSV / Excel) |

**All plots include:**
- Configurable ggplot2 theme (Classic, Minimal, BW, Void, Light)
- Custom title, font size, legend toggle
- Color picker for continuous scales, palette selector for categorical data
- Export to **PNG / PDF / SVG** with user-defined width × height

---

## Input Data

CellFIE accepts the following file formats:

| Format | Object type |
|---|---|
| `.rds` | `SingleCellExperiment` (Bioconductor) |
| `.rds` | `Seurat` object (auto-converted to SCE on load) |
| `.qs` | `SingleCellExperiment` or `Seurat` (via `qs2`) |

> **Maximum upload size:** 5 GB (configurable in `ui.R`)

The object must contain at least:
- One or more **assays** (e.g., `counts`, `logcounts`, `normcounts`)
- At least one **reduced dimension** (e.g., `UMAP`, `TSNE`, `PCA`)
- **Cell metadata** (`colData`) for grouping and annotation

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/jdutel42/CellFIE.git
cd CellFIE
```

### 2. Install R dependencies

Open R (≥ 4.2) and run:

```r
# CRAN packages
install.packages(c(
  "shiny", "bslib", "ggplot2", "dplyr", "tidyr",
  "shinycssloaders", "shinyjs", "shinyFeedback", "shinyalert",
  "colourpicker", "DT", "openxlsx", "reshape2", "scales"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "SingleCellExperiment",
  "SummarizedExperiment"
))

# qs2 (fast serialization)
install.packages("qs2")

# Seurat (optional — only required if loading Seurat objects)
install.packages("Seurat")
```

---

## Usage

### Launch locally

```r
shiny::runApp(
  appDir = "path/to/CellFIE",
  launch.browser = TRUE
)
```

### Launch on a server (headless)

```r
shiny::runApp(
  appDir = "~/CellFIE",
  host   = "0.0.0.0",
  port   = 3838,
  launch.browser = FALSE
)
```

The app will then be accessible at `http://<server-ip>:3838` from any machine on the same network.

---

## Application Structure

```
CellFIE/
├── ui.R          # User interface — layout, panels, input widgets
├── server.R      # Server logic — reactives, plot generation, download handlers
└── README.md
```

> CellFIE follows a standard two-file Shiny architecture. The `ui.R` and `server.R` files are fully modular: each visualization module (FeaturePlot, ViolinPlot, etc.) is self-contained with its own reactive data pipeline and download handlers.

---

## Module Details

### 🔬 FeaturePlot

Visualizes gene expression (or any numeric/categorical metadata) projected onto a 2D embedding (UMAP, t-SNE, PCA, etc.).

- Up to **6 genes** displayed simultaneously via `facet_wrap`
- Continuous scale: custom low/high color via color picker
- Categorical scale: RColorBrewer palette selector
- Point size and alpha transparency controls

### 🎻 ViolinPlot

Displays the distribution of expression values across user-defined cell groups.

- Up to **8 genes** in faceted layout
- Optional **boxplot** overlay (IQR + median)
- Optional **jitter** (individual cell points)
- Axis flip for readability with many groups

### 🔵 DotPlot

Classic scRNA-seq summary plot: dot **size** encodes the percentage of expressing cells; dot **color** encodes mean expression.

- Supports up to **30 genes**
- Optional **hierarchical clustering** of genes by expression profile
- Multiple continuous color palettes (Blues, RdYlBu, Viridis, Plasma, YlOrRd)

### 🧬 Cell Explorer

Projects any cell metadata (cluster labels, sample origin, cell type annotation, etc.) onto the embedding.

- **Centroid labels** auto-computed from median coordinates per group
- Sidebar **composition barplot** and **count table** per group
- Supports both categorical and continuous metadata

### 🧪 QC

Assesses dataset quality using numeric metadata columns (typically `nCounts`, `nFeatures`, `pct_mito`, etc.).

- **Violin + boxplot** for each selected metric
- **Scatter plot** between any two metrics with linear regression line
- Color-by option using any categorical column (e.g., sample, batch)

### 📋 Metadata Table

Fully interactive table of `colData` metadata.

- Column selection, column filtering, free-text search
- Export to **CSV** or **Excel** (`.xlsx`)
- Built with `DT` (DataTables) with copy-to-clipboard support

---

## Export Options

Every visualization module provides three download formats with configurable dimensions:

| Format | Use case |
|---|---|
| **PNG** | Presentations, lab meetings, quick sharing |
| **PDF** | Vector format for manuscripts and figures |
| **SVG** | Fully editable vector format (Inkscape, Illustrator) |

Width and height are set in pixels (PNG) and auto-converted to inches for PDF/SVG.

---

## Dependencies

| Package | Source | Role |
|---|---|---|
| `shiny` | CRAN | Web application framework |
| `bslib` | CRAN | Bootstrap 5 theming |
| `ggplot2` | CRAN | Plot generation |
| `dplyr` / `tidyr` | CRAN | Data manipulation |
| `SingleCellExperiment` | Bioconductor | SCE object handling |
| `SummarizedExperiment` | Bioconductor | `colData` / `assay` access |
| `qs2` | CRAN | Fast `.qs` file I/O |
| `Seurat` | CRAN | Optional Seurat → SCE conversion |
| `colourpicker` | CRAN | Interactive color selection |
| `DT` | CRAN | Interactive data tables |
| `openxlsx` | CRAN | Excel export |
| `shinycssloaders` | CRAN | Loading spinners |
| `shinyjs` | CRAN | JavaScript helpers |
| `shinyFeedback` | CRAN | Input validation feedback |
| `shinyalert` | CRAN | Alert dialogs |
| `reshape2` | CRAN | Data reshaping (DotPlot) |
| `scales` | CRAN | Color scale utilities |

---

## Recommended Workflow

```
1. Prepare your SingleCellExperiment object in R
   └── Ensure it contains: assays (logcounts), reducedDims (UMAP/TSNE), colData (metadata)

2. Save as .rds or .qs
   └── qs2::qsave(sce, "my_dataset.qs")   # faster I/O for large objects

3. Launch CellFIE and upload your file

4. Explore:
   ├── FeaturePlot  → marker genes on embedding
   ├── ViolinPlot   → expression by cluster
   ├── DotPlot      → gene panel × cluster heatmap
   ├── Cell Explorer → annotation, batch, sample composition
   └── QC           → nCounts, nFeatures, pct_mito distributions

5. Customize and download figures (PNG / PDF / SVG)
```

---

## Notes for Bioinformaticians

- **Seurat objects** are automatically converted to `SingleCellExperiment` on load using `Seurat::as.SingleCellExperiment()`. All assays, reductions, and metadata are preserved.
- **Large datasets** (> 100k cells): enable the *Rasterize* option in FeaturePlot to avoid rendering overhead. Consider pre-downsampling for interactive use.
- **Assay selection**: CellFIE dynamically lists all available assays. For expression visualization, `logcounts` (log-normalized) is recommended over raw `counts`.
- **colData requirements**: QC and grouping features rely entirely on `colData(sce)`. Ensure meaningful columns are present (cell type, cluster, sample, batch, QC metrics).

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Author

Developed as part of a single-cell analysis project.  
Contributions, issues and pull requests are welcome.

---

<div align="center">
<sub>Built with R · Shiny · Bioconductor · ggplot2</sub>
</div>
