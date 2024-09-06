# Extending palettes
extendPaletteFUN <- colorRampPalette(RColorBrewer::brewer.pal(7, "Spectral"))

# Reference for making Sankey diagram using flipPlots package
edges <- as.data.frame(table(as.data.frame(is_geqthreshold_mx)))
  criteria_len <- ncol(edges)
  for (i in 1:(criteria_len - 1)) {
    edges[[i]] <- as.character(edges[[i]])
  }
  SankeyDiagram(edges[, -(criteria_len)], link.color = "Source", weights = edges$Freq, max.categories = 23, hovertext.show.percentages = TRUE, label.show.counts = TRUE, font.size = 10, variables.share.values = TRUE, label.show.percentages = TRUE)

# Collection of colours mostly distinguishable  
# From https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pie(rep(1, 25), col = c25)
