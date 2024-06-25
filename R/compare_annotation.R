# Function/s to compare reference-based annotation using reference dataset from SingleR and Seurat label transfer methods
# Goal is to output a robust reference-based annotation combining results from those methods

library(cowplot)
library(DelayedArray)
library(DelayedMatrixStats)
library(flipPlots)
library(ggplot2)
library(reshape2)

# Function 1: Get difference between highest value (or value for the given label if labels provided) and median value

# Based from SingleR::getDeltaFromMedian() but parameters changed to generalise for any matrix e.g. prediction scores from Seurat transfer data, and to break down the code to make it more understandable
# Confirmed working by comparing output to output of SingleR::getDeltaFromMedian()
# Allowed elements for labels vector are the column names of scores_mx

#library(DelayedArray)
get_delta_frommedian <- function(scores_mx, labels = NULL) {
  
  if (!is.matrix(scores_mx)) {
    stop("scores_mx is not a matrix.")
  } else if (is.null(labels)) {
    message("Provide labels argument. Still need to add code that works without labels provided i.e. it determines max value per row. But the current code does allow you to provide any set of labels (doesn't have to be the ones with max value, although this is the expected usage, and compare their value with median.")
  } else {
    # 2-column matrix - 1st column is row index, 2nd column is column index of label
    index_mx <- cbind(seq_along(labels), match(labels, colnames(scores_mx)))
    # The matrix is used to vectorise the selection of score per label
    assigned <- scores_mx[index_mx]
    return(assigned - DelayedMatrixStats::rowMedians(DelayedArray(scores_mx)))
  }
  
}

# Function 2: Get difference between highest value and the second highest value

# Confirmed working when using iris[,1:2] i.e. delta equal to absolute difference between 2 columns
get_delta_fromnext <- function(scores_mx) {
  
  first_second_scores_mx <- apply(scores_mx, 1, FUN = function(row_scores){
    sorted <- sort(unname(row_scores), decreasing = TRUE)[1:2]
    names(sorted) <- c("first", "second")
    return(sorted)
  })
  #first_second_scores_mx <- t(first_second_scores_mx)
  # Output of apply is a transposed version of scores_mx, costly to do t() so just index accordingly:
  return(abs(first_second_scores_mx[1,] - first_second_scores_mx[2,]))
  
}

# Function 3: Generate a dataframe with features from both methods that can be used to get a robust annotation

# Confirmed working by manually calculating each column from source then comparing with output columns of function
# Prediction quality features:
#   a. **delta.med** and **delta.next** - difference of highest prediction score to median and second highest, then determine those with really low deltas based on outlier approach
#   b. **prediction.score.max** - Also use a threshold for prediction score i.e. a high confidence prediction should have high prediction score and high delta.next

combine_annotation_features <- function(singler_obj, transfer_obj) {
  
  # Break down objects into parts for downstream processing
  
  singler_scores_mx <- as.matrix(singler_obj$scores)
  transfer_scores_mx <- as.matrix(transfer_obj[, setdiff(colnames(transfer_obj), c("predicted.id", "prediction.score.max"))])
  rownames(transfer_scores_mx) <- NULL
  
  singler_labels <- singler_obj$labels
  singler_labels_pruned_nareplaced <- singler_obj$pruned.labels
  singler_labels_pruned_nareplaced[is.na(singler_obj$pruned.labels)] <- "pruned"
  
  transfer_labels <- paste0("prediction.score.", transfer_obj$predicted.id)
  transfer_labels_prefix_removed <- vapply(strsplit(transfer_labels, split = "prediction.score.", fixed = TRUE),
                                           function(strsplit_out) strsplit_out[2], FUN.VALUE = character(1))
  
  # Check input objects and objects derived from them are valid
  
  is_valid_objs <- c(identical(rownames(singler_obj), rownames(transfer_obj)), # Identical barcodes
                     # Make sure reference labels from both methods are from the same source
                     all(singler_labels %in% colnames(singler_scores_mx)),
                     all(transfer_labels_prefix_removed %in% colnames(singler_scores_mx))
                     )
  
  if (!all(is_valid_objs)) {
    stop("combine_annotation_features(): Invalid input objects.")
  }
  
  # Build dataframe with features to diagnose annotation results
  
  features_df <- data.frame(
    
    labels_singler_pruned = singler_labels_pruned_nareplaced,
    labels_singler = singler_labels,
    labels_transfer = transfer_labels_prefix_removed,
    labels_singler_transfer = paste0(singler_labels, "-", transfer_labels_prefix_removed),
    labels_singler_transfer_sorted = apply(cbind(singler_labels, transfer_labels_prefix_removed), MARGIN = 1, function(labels_vector) {
      paste(sort(labels_vector, decreasing = FALSE), collapse = "-")
    }),
    labels_combined = ifelse(singler_labels == transfer_labels_prefix_removed, yes = singler_labels, no = "Ambiguous"),
    
    delta.med_singler = get_delta_frommedian(singler_scores_mx, singler_labels),
    delta.med_transfer = get_delta_frommedian(transfer_scores_mx, transfer_labels),
    
    delta.next_singler_existing = singler_obj$delta.next,
    delta.next_singler = get_delta_fromnext(singler_scores_mx),
    delta.next_transfer = get_delta_fromnext(transfer_scores_mx),
    
    max_score_singler = apply(singler_scores_mx, 1, max), # These are correlation-based values that are very sensitive to technical factors as detailed in SingleR book that is why they used delta features instead to diagnose annotation
    max_score_transfer = transfer_obj$prediction.score.max
    
  )
  #summary(features_df)
  
  return(features_df)
  
}

# Function 4: Function to plot feature distributions

# library(cowplot)
# library(ggplot2)
# library(reshape2)
plot_features_distribution <- function(df, grouping_var){
  
  numeric_colnames <- colnames(df)[unlist(vapply(df, is.numeric, logical(1)))]
  tidy_df <- reshape2::melt(df[,c(grouping_var, numeric_colnames)], id = grouping_var)
  p_lst <- NULL
  
  p <- ggplot(tidy_df, aes(.data[[grouping_var]], value)) +
    geom_violin() +
    geom_boxplot() +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~variable, ncol = 2)
  
  #print(plot_grid(plotlist = p_lst, ncol = 1))
  return(p)
  
}

# Function 5: Function to identify robust annotations based on thresholds for annotation features to be used for selection

# min_thresholds can be a function (with numeric(1) as output) or vector of threshold values each corresponding to features in features_touse
# library(flipPlots)
#apply_out <- apply_feature_thresholds(combine_df[combine_df$labels_combined == "Microglia",], features_touse = c("delta.med_singler", "max_score_transfer"), min_thresholds = quantile,probs = 0.15)
apply_feature_thresholds <- function(df, features_touse = NULL, min_thresholds = quantile, ...){
  
  # Prepare matrix of features for selection
  numeric_columns_mx <- as.matrix(df[,unlist(vapply(df, is.numeric, logical(1)), use.names = FALSE)])
  #numeric_columns_mx <- numeric_columns_mx[df[[grouping_var]] == group,]
  
  if (!is.null(features_touse)) {
    numeric_columns_mx <- numeric_columns_mx[,features_touse, drop = FALSE]
  } else {
    message("Using all numeric columns for selection...")
  }
  
  # Calculate threshold for each feature based on min_thresholds
  if (is.function(min_thresholds)) {
    min_thresholds_vals <- apply(numeric_columns_mx, 2, FUN = function(x) {
      return(vapply(list(x), min_thresholds, numeric(1), ...))
    })
  } else {
    min_thresholds_vals <- setNames(min_thresholds, nm = features_touse)
    message("Using input min_thresholds as minimum thresholds...")
  }
  
  message("The minimum thresholds used are:")
  print(min_thresholds_vals)
  
  if (length(min_thresholds_vals) != length(features_touse)) {
    stop("apply_feature_thresholds: Number of minimum thresholds not equal to number of features to use.")
  }
 
  # Identify which annotations are robust and not based on minimum threshold value per feature
  
  # comparison gives me wrong behavior, 1s not equal to 1s - some numbers might appear 1 but they could be 0.999.. with more than seven 9s so the value appear as 1 i.e. just trust the resulting comparison but do check it works for other decimal values
  # **TRY with other threshold**
  min_thresholds_vals_mx <- matrix(unname(min_thresholds_vals), nrow = nrow(numeric_columns_mx), ncol = ncol(numeric_columns_mx), byrow = TRUE,
                                   dimnames = dimnames(numeric_columns_mx))
  is_geqthreshold_mx <- numeric_columns_mx >= min_thresholds_vals_mx
 
  # Use Sankey diagram to show how many have robust annotation based on input thresholds for each feature
  edges <- as.data.frame(table(as.data.frame(is_geqthreshold_mx)))
  if (length(features_touse) > 1) {
    
    criteria_len <- ncol(edges)
    for (i in 1:(criteria_len - 1)) {
      edges[[i]] <- as.character(edges[[i]])
    }
    plot_sankey <- SankeyDiagram(edges[, -(criteria_len)], link.color = "Source", weights = edges$Freq, max.categories = 23, hovertext.show.percentages = TRUE, label.show.counts = TRUE, font.size = 10, variables.share.values = TRUE, label.show.percentages = TRUE)
    plot_sankey
    
  } else {
    plot_sankey <- NULL
  }
 
  # Output
  is_allgeqthreshold <- rowSums(is_geqthreshold_mx) == ncol(is_geqthreshold_mx)
  output <- list(is_allgeqthreshold = is_allgeqthreshold, is_geqthreshold_mx = is_geqthreshold_mx, min_thresholds_vals = min_thresholds_vals, tally = edges, plot_sankey = plot_sankey)
  return(output)
  
}
