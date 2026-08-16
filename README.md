# R utilities

An incubator for reusable R functions developed across analysis projects. Code
is grouped by purpose while it is refined, tested, and prepared for possible
promotion into a stable package.

## Collections

| Folder | Intended contents and examples |
| --- | --- |
| `linear-modelling/` | Linear and mixed models, contrasts, model checks, and diagnostics; for example, fitting an LMM and extracting estimated marginal means. |
| `imputation/` | Imputation methods and their evaluation; for example, benchmarking imputation error at different missing-data rates. |
| `distribution-testing/` | Tests and diagnostics for feature distributions; for example, normality and pre/post-processing distribution shifts. |
| `correlation/` | Correlation and dependence analysis; for example, estimating correlations robustly across repeated subsamples. |
| `enrichment/` | Gene-set, pathway, and ranked-set enrichment; for example, GO enrichment or an mHG test. |
| `tf-activity/` | Transcription-factor activity analysis and wrappers for methods such as BITFAM. |
| `cell-annotation/` | Cell-label assignment, comparison, confidence measures, and annotation diagnostics. |
| `data-integration/` | Batch-effect assessment and integration of datasets or representations. |
| `quality-control/` | Sample, cell, and feature quality control, including outlier detection and filtering. |
| `clustering/` | Clustering, cluster assessment, and related visualisation helpers. |
| `data-conversion/` | Conversion between data structures; for example, Seurat and `SingleCellExperiment` objects. |
| `trajectory-analysis/` | Trajectory inference and downstream trajectory-associated testing, including Slingshot and tradeSeq helpers. |
| `reference-data/` | Functions for retrieving, loading, and preparing reference datasets. |
| `marker-sets/` | Curated marker genes and reusable biological signature definitions. |
| `plotting/` | General plotting functions, themes, palettes, and plot composition helpers. |
| `environment/` | R library, Conda, reticulate, and session-environment setup. |
| `parallel-computing/` | Parallel backend and worker configuration. |
| `file-handling/` | Reusable directory, filename, and file-management helpers. |
| `general/` | Small, genuinely general helpers that do not yet belong to a more specific collection. |

## Current contents

The initial collections cover cell annotation, integration, quality control,
reference datasets, functional enrichment, trajectory analysis, plotting,
environment setup, marker sets, and transcription-factor activity. Mixed
scripts remain intact until their functions can be reviewed and separated
without changing behaviour.

Interfaces in this repository may change while utilities are being developed.
Projects consuming them should pin a known Git commit or submodule revision.
