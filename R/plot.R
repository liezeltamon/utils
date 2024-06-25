# Reference for making Sankey diagram using flipPlots package
edges <- as.data.frame(table(as.data.frame(is_geqthreshold_mx)))
  criteria_len <- ncol(edges)
  for (i in 1:(criteria_len - 1)) {
    edges[[i]] <- as.character(edges[[i]])
  }
  SankeyDiagram(edges[, -(criteria_len)], link.color = "Source", weights = edges$Freq, max.categories = 23, hovertext.show.percentages = TRUE, label.show.counts = TRUE, font.size = 10, variables.share.values = TRUE, label.show.percentages = TRUE)
  
