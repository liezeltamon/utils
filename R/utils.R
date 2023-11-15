# No fix for levels of same size
as_factor_bySize <- function(x, decreasing = TRUE){
  x <- as.character(x)
  ordered_levels <- sort(table(x), decreasing = decreasing)
  ordered_levels <- names(ordered_levels)
  x <- factor(x, levels = ordered_levels)
  return(x)
}