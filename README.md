# Single-Cell Analysis Workflow for KL-TRP63 Lung Adenocarcinoma-to-Squamous Carcinoma Transdifferentiation (AST)

An R/Seurat-based single-cell RNA sequencing (scRNA-seq) analysis workflow for lung organoids. This project investigates the dynamics of TRP63-positive cells during lung adenocarcinoma (LUAD)-to-lung squamous cell carcinoma (LUSC) transdifferentiation (Adeno-to-Squamous Transition, AST), identifies the phenotypic transition critical point (Tipping Point) using Dynamic Network Biomarker (DNB) theory, and evaluates the clinical value of NRG1 and the `trajectory_up` signature.

## Core Concepts

| Concept | Description |
|------|------|
| `tdTomato_ratio` | Proportion of TRP63 expression (0 = LUAD phenotype, 1 = LUSC phenotype), calculated from the Amplicon-Seq fluorescent barcode |
| `GFP` / `tdTomato` | Cell populations with low / high TRP63 expression |
| `AST` | Trajectory subpopulation undergoing adeno-to-squamous transition |
| `DNB` | Dynamic network biomarker, CI = (PCC_in / PCC_out) × SD_in |
| `Tipping Point` | Critical state of phenotypic transition; the states before and after it are Pre-Tipping and Post-Tipping |

## Directory Structure

```
repository/
├── README.md                        # This file
├── R/                               # Utility functions and color schemes
│   ├── custom_functions.R           # General plotting and helper functions
│   └── custom_colors.R              # Unified color configuration
├── scripts/                         # Analysis scripts (run in numbered order)
│   ├── 01_KLTN_process.R            # Process the KLTN dataset from this experiment
│   ├── 02_integrated_analysis.R     # Integrate single-cell data from other studies
│   ├── 03_trajectory_subset_analysis.R  # AST trajectory subpopulation analysis (including top-100 associated-gene clustering/GO/network analysis for FigS4)
│   ├── 04_DNB_analysis.R            # DNB critical-point analysis
│   ├── 05_survival_analysis.R       # Clinical value analysis (including EGFR/KRAS signature boxplots for FigS5)
│   └── cellrank.ipynb               # CellRank trajectory-potential analysis (Python, used with script 03)
├── top100_associate_cytoscape/       # FigS4 B top-100 associated-gene GeneMANIA/Cytoscape network workflow
├── data/                            # Input data (see "Input Data" below)
└── output/                          # Figures and RDS outputs
    ├── figures/                     # Figures (organized by analysis module)
    ├── object/                      # RDS objects
    └── tables/                      # Analysis tables
```

## Analysis Workflow and Execution Order

The scripts assume that they are run from the repository root (`DIR_ROOT <- getwd()`). Do not source them from a different working directory.

| Order | Script | Content | Key Outputs |
|------|------|------|----------|
| 1 | `01_KLTN_process.R` | Read 10X h5 → QC → remove doublets with DoubletFinder → clustering → detailed immune/epithelial annotation → integrate Amplicon-Seq fluorescence (`tdTomato_ratio`) → CytoTRACE2 → inferCNV | `se.rds`, subset RDS objects, CNV scores |
| 2 | `02_integrated_analysis.R` | Unified QC of the KLTN + SCC + KP-LUAD datasets → SCT + Harmony integration → reference annotation mapping (symphony K-NN) → AST transitional-state validation → GMT pathway scoring | Integrated RDS objects and `cell_annotation` figures |
| 3 | `03_trajectory_subset_analysis.R` | AST trajectory subpopulation: tdRatio temporal genes → slingshot trajectory → tradeSeq pseudotime → SCENIC regulatory network → signature scoring → top-100 associated-gene clustering/GO/network analysis (FigS4) | `se_trajectory_subset.rds`, SCENIC results, top-100 heatmaps/GO figures |
| 3.5 | `cellrank.ipynb` (Python) | Combined PseudotimeKernel + ConnectivityKernel → GPCCA macrostates/fate probabilities → C1/C2 lineage-driver genes | `C1/C2_drivers.csv` (used by Chapter 5 of script 03) |
| 4 | `04_DNB_analysis.R` | DNB sliding windows → Tipping Point identification → Pre/Tipping/Post-state DEGs → dynamic TF-target-gene networks | `dnb_state` column and `DNB_final_df.csv` |
| 5 | `05_survival_analysis.R` | KM survival analysis of NRG1 and the `trajectory_up` signature (TCGA-LUAD/LUSC, KRYSTAL-1) → EGFR/KRAS pathway correlations → EGFR/KRAS signature boxplots (FigS5) | Survival curves, Cox results, and signature boxplots |

### Data Flow Between Scripts (`generated_data` Notes)

The project has two core intermediate data objects (paths originally recorded in `src/generated_data.R`):

- `output/object/seurat/se.rds` — produced by script 01 (the complete KLTN annotation object, including fluorescence and CNV information) → used by scripts 02 and 03
- `output/object/seurat/se_trajectory_subset.rds` — produced by script 03 (the AST trajectory subpopulation) → used by script 04, with additional score columns added by subsequent analyses (script 02 Chapter 10 and script 04)

Note the execution order: some chapters of script 02 depend on outputs from script 03 (`association_results.csv` and `se_trajectory_subset.rds`). The recommended complete order is 01 → 03 → 02 → 04 → 05, or follow the “Input” notes at the beginning of each script.

Also note that `cellrank.ipynb` (Python) must be run after the slingshot chapter of script 03. It reads `se_trajectory_subset.h5ad` and the slingshot pseudotime CSV. Its output, `C1/C2_drivers.csv`, is an input to Chapter 5 of script 03 (TCGA survival analysis).

## Input Data

**See [data/README.md](data/README.md) for the complete data inventory, md5 checksums, and large-file download links.**

```
data/
├── README.md                        # Data inventory and download instructions (uses/size/checksums for 38 files)
├── MD5SUMS.txt                     # md5 checksums for all data files
├── 10xgenomics/                    # Raw 10X matrix from the KLTN experiment
├── GFP_algin_score.txt             # GFP amplicon alignment results (>100MB, externally hosted)
├── tdTomato_algin_score.txt        # tdTomato amplicon alignment results (>100MB, externally hosted)
├── geneset.hs.gmt                  # Human gene sets (GMT format)
├── scc_tong_2026/                  # SCC dataset
├── gene_orthology/                 # Human/mouse signature lists (script 05 Chapters 2/5)
├── clinical_validation/            # KRAS/EGFR clinical validation data + mouse organoid TPM (script 05 Chapters 3/5)
└── precomputed/                    # Intermediate objects from upstream processing (mirrors output/object/;
                                    # copy them to output/object/ to continue from intermediate results)
```

Data acquisition notes:

- **Small files (<100MB, 29 files):** distributed directly with the repository
- **Large files (>100MB, 9 files, approximately 14.3GB):** exceed GitHub’s per-file limit and are listed in `.gitignore`; download them from external hosting (Zenodo links are provided in `data/README.md`) and place them in the corresponding paths under `data/`
- **Using intermediate objects:** after running `cp -r data/precomputed/* output/object/`, each script (with switches set to `FALSE`) can continue from the intermediate results without rerunning upstream steps

External public data (not distributed with the repository; official download links are in `data/README.md`):

- KPTracer: `PATH_KP_LUAD_RDS` (declared at the top of script 02) — `adata_processed.nt.rds` from KPTracer-Data
- SCENIC mm10 ranking database (feather): the `SCENIC_ref` path configured in script 03
- SCENIC helper `aux_createViewerApp.R` and KRYSTAL-1 `mmc2/mmc3.xlsx`: included in `data/precomputed/`

## Logic Switches

Set the corresponding switch to `TRUE` when rerunning time-consuming sections. With the default value `FALSE`, previously saved results are read:

| Script | Switch | Function |
|------|------|------|
| 01 | `RUN_DOUBLETFINDER` | Fully rerun QC + DoubletFinder from h5 |
| 02 | `RUN_SCC_PREPROCESS` | Fully preprocess the raw SCC h5 data |
| 03 | `RUN_TCGA_DOWNLOAD` | Download the TCGA data again |
| 05 | `RUN_KRYSTAL_PREPROCESS` / `RUN_EXPORT_SIGNATURE_LONG_CSV` | KRYSTAL-1 preprocessing / export the signature |

## Dependencies

Core: Seurat, harmony, tidyverse (dplyr/ggplot2/tidyr/stringr/purrr), patchwork, clustree, DoubletFinder, future  
Trajectory: slingshot, tradeSeq, CytoTRACE2  
Regulatory networks: SCENIC (GENIE3/AUCell/RcisTarget), igraph  
Pathways: GSVA/ssGSEA, clusterProfiler  
Clinical: survival, survminer, timeROC, rstatix  
Other: pheatmap, ggalluvial, reshape2, SCP, symphony, presto, scico, viridis, ggsci, RColorBrewer, org.Mm.eg.db

> Version note: The original environment used R 4.4. If igraph and Seurat have version conflicts, load igraph first (this is handled in script 02). Some packages (DoubletFinder, SCP, CytoTRACE2) must be installed from GitHub.

Python (`scripts/cellrank.ipynb`): cellrank, scanpy, scvelo, anndata, pandas, numpy; `cr.models.GAMR` requires rpy2 and calls R’s mgcv package.
