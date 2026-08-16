compute_group_prop_complete <- function(
    df,
    grouping_keys,
    counted_key,
    as_percentage = FALSE,
    counted_universe = NULL,
    remove_na = TRUE # Will include NA in category if remove_na = FALSE
) {
    required_cols <- c(grouping_keys, counted_key)
    missing_cols <- setdiff(required_cols, names(df))

    if (length(missing_cols) > 0) {
        stop(
            "Missing required columns: ",
            paste(missing_cols, collapse = ", ")
        )
    }

    if (remove_na) {
        df_count <- df %>%
            dplyr::filter(!is.na(.data[[counted_key]]))
    } else {
        df_count <- df
    }

    universe_supplied <- !is.null(counted_universe)

    if (universe_supplied) {
        if (anyDuplicated(counted_universe)) {
            stop("counted_universe contains duplicate categories.")
        }
    } else {
        counted_universe <- unique(df_count[[counted_key]])
    }

    if (length(counted_universe) == 0) {
        stop("No retained categories found in ", counted_key, ".")
    }

    grouping_df <- df_count %>%
        dplyr::distinct(dplyr::pick(dplyr::all_of(grouping_keys)))

    if (nrow(grouping_df) == 0) {
        stop("No grouping combinations remain after filtering.")
    }

    # Match the Python categorical behavior: when a universe is supplied,
    # values outside it are ignored.
    if (universe_supplied) {
        count_input <- df_count %>%
            dplyr::filter(.data[[counted_key]] %in% counted_universe)
    } else {
        count_input <- df_count
    }

    counts_df <- count_input %>%
        dplyr::count(
            dplyr::pick(
                dplyr::all_of(c(grouping_keys, counted_key))
            ),
            name = "count"
        )

    prop_df <- tidyr::crossing(
        grouping_df,
        !!counted_key := counted_universe
    ) %>%
        dplyr::left_join(
            counts_df,
            by = c(grouping_keys, counted_key)
        ) %>%
        dplyr::mutate(
            count = tidyr::replace_na(count, 0L),
            prop = count / sum(count),
            .by = dplyr::all_of(grouping_keys)
        )

    category_group_counts <- prop_df %>%
        dplyr::count(
            .data[[counted_key]],
            name = "n_groups"
        )

    if (dplyr::n_distinct(category_group_counts$n_groups) != 1L) {
        stop(
            "Error in completing the data frame: not all ",
            counted_key,
            " levels appear in all groups."
        )
    }

    if (anyNA(prop_df$prop) || any(is.nan(prop_df$prop))) {
        stop(
            "Some groups have zero total count across all ",
            counted_key,
            " categories."
        )
    }

    if (as_percentage) {
        prop_df <- prop_df %>%
            dplyr::mutate(prop = prop * 100)
    }

    n_unique_groupings <- nrow(grouping_df)
    n_unique_counted <- length(counted_universe)

    expected_rows <- n_unique_groupings * n_unique_counted
    actual_rows <- nrow(prop_df)

    if (expected_rows != actual_rows) {
        stop(
            "Error in computing proportions: expected ",
            expected_rows,
            " rows (",
            n_unique_groupings,
            " unique groupings * ",
            n_unique_counted,
            " unique ",
            counted_key,
            "), but got ",
            actual_rows,
            " rows."
        )
    }

    prop_df
}
