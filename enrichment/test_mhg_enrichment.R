# Test whether a term set is enriched at the top of a ranked query list using
# the minimum hypergeometric (mHG) test (Eden et al., 2007, PLoS Comp Biol).
# At each prefix cutoff k = 1..max_k, a hypergeometric p-value measures whether
# the top-k items contain more term members than expected by chance. The minimum
# p-value across cutoffs is the mHG statistic, Bonferroni-corrected for the
# number of cutoffs tested.
test_mhg_enrichment <- function(query, term, min_hits = 2L, max_frac = 0.5) {
  # Validate inputs
  stopifnot(is.character(query), length(query) > 0, !anyNA(query), !anyDuplicated(query))
  stopifnot(is.character(term), length(term) > 0, !anyNA(term), !anyDuplicated(term))
  stopifnot(is.numeric(min_hits), length(min_hits) == 1, min_hits >= 1, min_hits == as.integer(min_hits))
  stopifnot(is.numeric(max_frac), length(max_frac) == 1, max_frac > 0, max_frac <= 1)

  # Mark which positions in the ranked query belong to the term set
  membership <- query %in% term
  B <- length(query)
  N <- sum(membership)
  n_term_total <- length(term)
  n_term_missing <- n_term_total - N

  # Template for early-exit results
  empty <- data.frame(
    status = NA_character_,
    n_query = B, n_term_total = n_term_total,
    n_term_in_query = N, n_term_missing = n_term_missing,
    min_hits = min_hits, max_k = NA_integer_,
    n_cutoffs_tested = 0L,
    best_k = NA_integer_, hits_at_k = NA_integer_,
    expected_at_k = NA_real_, fold_enrichment = NA_real_,
    scan_pvalue = NA_real_, bonferroni_pvalue = NA_real_,
    stringsAsFactors = FALSE
  )

  if (N < min_hits) {
    empty$status <- "too_few_term_in_query"
    return(empty)
  }

  # Determine how far down the ranking to scan (clamped to list length)
  max_k <- min(B, as.integer(floor(max_frac * B)))
  if (max_k < 1L) {
    empty$status <- "max_k_too_small"
    return(empty)
  }

  # Count cumulative hits at each cutoff and keep only those with enough hits
  cumhits <- cumsum(membership)
  ks <- seq_len(max_k)
  valid <- cumhits[ks] >= min_hits

  if (!any(valid)) {
    empty$status <- "no_valid_cutoffs"
    empty$max_k <- max_k
    return(empty)
  }

  valid_ks <- ks[valid]
  valid_hits <- cumhits[valid_ks]

  # Compute hypergeometric p-value at each valid cutoff
  pvals <- phyper(
    valid_hits - 1L,
    N,
    B - N,
    valid_ks,
    lower.tail = FALSE
  )

  # Select the cutoff with the minimum p-value (the mHG statistic)
  best_idx <- which.min(pvals)
  scan_pvalue <- pvals[best_idx]
  best_k <- valid_ks[best_idx]
  hits_at_k <- valid_hits[best_idx]
  expected_at_k <- best_k * N / B
  fold_enrichment <- hits_at_k / expected_at_k

  # Bonferroni-correct for scanning over multiple cutoffs
  n_cutoffs_tested <- length(valid_ks)
  bonferroni_pvalue <- min(1, scan_pvalue * n_cutoffs_tested)

  data.frame(
    status = "ok",
    n_query = B, n_term_total = n_term_total,
    n_term_in_query = N, n_term_missing = n_term_missing,
    min_hits = min_hits, max_k = max_k,
    n_cutoffs_tested = n_cutoffs_tested,
    best_k = best_k, hits_at_k = hits_at_k,
    expected_at_k = expected_at_k, fold_enrichment = fold_enrichment,
    scan_pvalue = scan_pvalue, bonferroni_pvalue = bonferroni_pvalue,
    stringsAsFactors = FALSE
  )
}
