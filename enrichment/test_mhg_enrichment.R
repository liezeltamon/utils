# Test whether a term set is enriched at the top of a ranked query list using
# either the current manual mHG-style scan with Bonferroni correction or the
# exact XL-mHG implementation from mhg::mhg_test.
test_mhg_enrichment <- function(
    query,
    term,
    min_hits = 2L,
    max_frac = 0.5,
    method = c("xlmhg_package", "manual"),
    boundary_policy = c("set_p_to_1", "asis"),
    # Refer to https://github.com/slowkow/xlmhg/blob/master/src/mhg.cpp for xlmhg_package parameters:
    upper_bound = FALSE, # Instead of running a dynamic programming algorithm, return the upper bound for the p-value.
    tol = 1e-16
  ) {
  method <- match.arg(method)
  boundary_policy <- match.arg(boundary_policy)

  # Validate inputs
  stopifnot(is.character(query), length(query) > 0, !anyNA(query), !anyDuplicated(query))
  stopifnot(is.character(term), length(term) > 0, !anyNA(term), !anyDuplicated(term))
  stopifnot(is.numeric(min_hits), length(min_hits) == 1, min_hits >= 1, min_hits == as.integer(min_hits))
  stopifnot(is.numeric(max_frac), length(max_frac) == 1, max_frac > 0, max_frac <= 1)
  stopifnot(is.logical(upper_bound), length(upper_bound) == 1, !is.na(upper_bound))
  stopifnot(is.numeric(tol), length(tol) == 1, tol > 0)

  method_upper_bound <- if (method == "xlmhg_package") upper_bound else NA

  # Mark which positions in the ranked query belong to the term set
  membership <- query %in% term
  B <- length(query)
  N <- sum(membership)
  n_term_total <- length(term)
  n_term_missing <- n_term_total - N

  empty <- data.frame(
    status = NA_character_,
    method = method,
    boundary_policy = boundary_policy,
    upper_bound = method_upper_bound,
    n_query = B,
    n_term_total = n_term_total,
    n_term_in_query = N,
    n_term_missing = n_term_missing,
    min_hits = min_hits,
    max_k = NA_integer_,
    n_cutoffs_tested = 0L,
    best_k = NA_integer_,
    hits_at_k = NA_integer_,
    expected_at_k = NA_real_,
    fold_enrichment = NA_real_,
    mhg_statistic = NA_real_,
    pvalue = NA_real_,
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
    empty$max_k <- max_k
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
  n_cutoffs_tested <- length(valid_ks)

  if (method == "manual") {
    # Compute hypergeometric p-value at each valid cutoff
    pvals <- phyper(
      valid_hits - 1L, # Number of successes among those drawn genes, -1L to include the current number of hits in the p-value calculation (i.e. "at least this many")
      N, # Number of successes in the population (term genes in query)
      B - N, # Number of failures in the population (non-term genes in query)
      valid_ks, # Number drawn
      lower.tail = FALSE # Probability of seeing at least this many hits among those drawn genes
    )

    # Select the cutoff with the minimum p-value (the mHG statistic)
    best_idx <- which.min(pvals)
    mhg_statistic <- pvals[best_idx]
    best_k <- valid_ks[best_idx]
    hits_at_k <- valid_hits[best_idx]

    # Bonferroni-correct for scanning over multiple cutoffs
    pvalue <- min(1, mhg_statistic * n_cutoffs_tested)
  } else if (method == "xlmhg_package") {
    if (!requireNamespace("mhg", quietly = TRUE)) {
      stop("method = 'xlmhg_package' requires the mhg package")
    }
    x <- as.integer(membership)
    res <- mhg::mhg_test(
      x,
      B,
      N,
      max_k,
      min_hits,
      upper_bound = upper_bound,
      tol = tol
    )
    best_k <- as.integer(res$threshold)
    hits_at_k <- if (best_k > 0L) sum(x[seq_len(best_k)]) else 0L
    mhg_scores <- as.numeric(res$mhg)
    mhg_statistic <- if (best_k > 0L && length(mhg_scores) >= best_k) {
      mhg_scores[[best_k]]
    } else if (length(mhg_scores) > 0L) {
      mhg_scores[[1]]
    } else {
      NA_real_
    }
    pvalue <- as.numeric(res$pvalue)
  } else {
    stop("test_mhg_enrichment(): Invalid method.")
  }

  # If term genes make up (N / B)% of the whole query list, then by chance you’d expect about 10% of any top slice (e.g. top best_k) to be term genes too, if members were randomly distributed.
  expected_at_k <- if (best_k > 0L) best_k * N / B else NA_real_
  fold_enrichment <- if (!is.na(expected_at_k) && expected_at_k > 0) {
    hits_at_k / expected_at_k
  } else {
    NA_real_
  }

  if (boundary_policy == "set_p_to_1" && best_k == max_k) {
    pvalue <- 1
  }

  data.frame(
    status = "ok",
    method = method,
    boundary_policy = boundary_policy,
    upper_bound = method_upper_bound,
    n_query = B,
    n_term_total = n_term_total,
    n_term_in_query = N,
    n_term_missing = n_term_missing,
    min_hits = min_hits,
    max_k = max_k,
    n_cutoffs_tested = n_cutoffs_tested,
    best_k = best_k,
    hits_at_k = hits_at_k,
    expected_at_k = expected_at_k,
    fold_enrichment = fold_enrichment,
    mhg_statistic = mhg_statistic,
    pvalue = pvalue,
    stringsAsFactors = FALSE
  )
}

# %%
