# TODO

## Split mixed utility collections

This restructuring is deferred because it is not required by Feature Flow.
Before moving functions, audit consumers in SLE and other repositories so
existing `source()` calls can be updated without silently changing behaviour.

### `general/utils.R`

- Move label-specificity checks to `cell-annotation/` and consolidate the
  duplicate implementation in `reference-data/reference.R`.
- Move factor-ordering helpers to a focused data-manipulation collection.
- Move Cell Ranger FASTQ naming helpers to `file-handling/`.
- Consolidate palette definitions under `plotting/`.
- Remove commented scratch examples or replace them with maintained examples.

### `general/wrapper.R`

- Move Seurat and `SingleCellExperiment` conversions to `data-conversion/`.
- Move droplet and single-cell diagnostics to `quality-control/`.
- Move cell-cycle scoring to a feature-scoring collection.
- Move denoised PCA helpers to a dimensionality-reduction collection.
- Move clustering and clustree helpers to `clustering/`.
- Move UpSet input conversion to `data-conversion/`.
- Replace source-time package attachment and Python imports with explicit,
  function-level dependencies.

### Other mixed scripts

- Split file, plotting, and parallel helpers currently held in
  `environment/env.R`.
- Review generic testing and plotting helpers in the trajectory scripts.
- Separate generic subsampling from reference-dataset loaders.
- Refactor `environment/reticulate.R` into portable functions or a maintained
  example without cluster-specific paths.

## Validation for future splits

- Preserve public function names initially.
- Update every known consumer in the same change as a function move.
- Add focused loading and behaviour tests before removing the original files.
- Keep `general/` small and remove it once all remaining functions have a
  clear, stable owner collection.
