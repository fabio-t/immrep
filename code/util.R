
set.seed(42)

library(methods)
library(divo)
library(vegan)
library(gplots)
library(venneuler)
library(plyr)
library(randomcoloR)
library(fmsb)
library(circlize)
library(RColorBrewer)
library(Rtsne)
library(matrixStats)

library(devtools)
install_github("fabio-t/mixed.utils")
library(mixed.utils)

heatmaps <- function(regexp, prefix = NULL, invert = F, cex = 0.9, tsne = F) {
  if (is.null(regexp)) {
    ## All
    make_full_heatmap(mids_counts, "horn", mid_names, prefix = "full", binary = F, types = "square", cexRow = cex, cexCol = cex, tsne = tsne)
    make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = cex, cexCol = cex)
    make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = cex, cexCol = cex)
  } else {
    # prefix is expected to be set, now
    columns_i <- grep(regexp, mid_names, perl = T, invert = invert)
    columns_n <- mid_names[columns_i]
    data <- mids_counts[columns_i]
    make_full_heatmap(data, "horn", columns_n, prefix = prefix, binary = F, types = "square", cexRow = cex, cexCol = cex, tsne = tsne)
    make_full_heatmap(data, "intersect", columns_n, prefix = prefix, types = "square", vlim = c(0, NA), cexRow = cex, cexCol = cex)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0(prefix, "_top0.85"), topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = cex, cexCol = cex)
  }
}

radarcharts <- function(regexp, dirname, save = F, invert = F, colours = NULL) {
  print(dirname)

  column_indices <- grep(regexp, mid_names, perl = T, invert = invert)

  column_names <- colnames(mids_counts)[column_indices]
  column_names_2 <- mid_names[column_indices]

  print(column_names_2)

  data <- mids_counts[column_names]

  dirname <- paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname, "/")

  l <- save_multiple_radarcharts(data, column_names, labels = column_names_2,
                                 topPerc = tp, percentage = T, use_log = F,
                                 dirname = dirname, use_colours = colours,
                                 show_axis = F, show_legend = T)

  if (save && (length(l$common_clones) > 0)) {
    # common clones
    d1 <- sweep(data, 2, colSums(data), "/")[l$common_clones, , drop=F]
    d2 <- data[l$common_clones, , drop=F]
    d3 <- as.data.frame(colRanks(as.matrix(d2), preserveShape = T, ties.method = "average"))
    colnames(d1) <- column_names_2
    colnames(d2) <- column_names_2
    colnames(d3) <- column_names_2
    rownames(d3) <- l$common_clones
    write.csv(d1, file = paste0(dirname, "common_clones_perc.csv"))
    write.csv(d2, file = paste0(dirname, "common_clones_counts.csv"))
    write.csv(d3, file = paste0(dirname, "common_clones_ranks.csv"))
  }

  if (save && (length(l$total_clones) > 0)) {
    # common clones
    d1 <- sweep(data, 2, colSums(data), "/")[l$total_clones, , drop=F]
    d2 <- data[l$total_clones, , drop=F]
    d3 <- as.data.frame(colRanks(as.matrix(d2), preserveShape = T, ties.method = "average"))
    colnames(d1) <- column_names_2
    colnames(d2) <- column_names_2
    colnames(d3) <- column_names_2
    rownames(d3) <- l$total_clones
    write.csv(d1, file = paste0(dirname, "total_clones_perc.csv"))
    write.csv(d2, file = paste0(dirname, "total_clones_counts.csv"))
    write.csv(d3, file = paste0(dirname, "total_clones_ranks.csv"))
  }
}

make_diversity <- function(mids_counts, mid_names) {
  divers_df <- as.data.frame(mid_names)
  rownames(divers_df) <- colnames(mids_counts)
  colnames(divers_df) <- "Sample"
  divers_df$Shannon <- diversity(t(mids_counts), index = "shannon")
  divers_df$ExpShannon <- exp(divers_df$Shannon)
  divers_df$InvSimpson <- diversity(t(mids_counts), index = "invsimpson")
  write.csv(divers_df, file = "diversity.csv")
}

elbow_plot <- function(pca, main = NULL) {
  plot(pca@R2, ylim = c(0, 1), xaxt = "n", xlab = NA, ylab = "R2", main = main)
  lines(pca@R2)
  axis(1, at = 1:length(pca@R2), labels = paste0("PC", 1:length(pca@R2)))
}

darken <- function(color, factor = 1.4) {
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  col
}

lighten <- function(color, factor = 1.4) {
  col <- col2rgb(color)
  col <- col * factor
  col <- rgb(t(col), maxColorValue = 255)
  col
}

load_denise_matrix <- function(filename = "matrix.csv", darken_on = NULL) {
  # rownames must be in MIDX format
  matrix <- read.csv(filename, row.names = 1)
  matrix$joined <- paste0(
    ifelse(is.na(matrix$Mouse), "", paste0(matrix$Mouse, "_")),
    matrix$Organ,
    ifelse(is.null(matrix$Day) | is.na(matrix$Day), "", paste0("_d", matrix$Day))
  )
  matrix$colours <- rgb(169, 169, 169, maxColorValue = 255)
  matrix[which(matrix$Organ == "Liver"), "colours"] <- rgb(255, 0, 0, maxColorValue = 255)
  matrix[which(matrix$Organ == "SPL"), "colours"] <- rgb(34, 139, 34, maxColorValue = 255)
  matrix[which(matrix$Organ == "MLN"), "colours"] <- rgb(30, 144, 255, maxColorValue = 255)
  matrix[which(matrix$Organ == "PP"), "colours"] <- rgb(218, 165, 32, maxColorValue = 255)
  matrix[which(matrix$Organ == "CC"), "colours"] <- rgb(148, 0, 211, maxColorValue = 255)

  # FIXME: needs to support more "organs".. or do it in a more general way
  # (we must still guarantee unique colour per organ across experiments, though)

  if (!is.null(darken_on)) {
    matrix[grep(darken_on, matrix$Mouse), "colours"] <- darken(matrix[grep(darken_on, matrix$Mouse), "colours"], 2)
  }

  print(matrix)
}

make_denise_plots <- function(exp, matrix, day = NULL, mouse = NULL, organ = NULL, dirname = NULL, heatmap = T, radarchart = T, barchart = T) {
  print(paste(day, mouse, organ))

  if (is.null(day)) {
    # if day is not specified, we get everything but not the Inoculum(s)
    logical <- !grepl("Inoculum", matrix$Organ, ignore.case = T)

    if (is.null(dirname)) {
      dirname <- paste0("plots/all")
    }
  }
  else {
    logical <- matrix$Day %in% day

    if (is.null(dirname)) {
      dirname <- paste0("plots/d", paste0(day, collapse = "-"))
    }
  }

  if (is.null(mouse) && is.null(organ)) {
    prefix <- "full"

    # we don't touch the logical vector
  }
  else if (is.null(mouse)) {
    prefix <- paste0(organ, collapse = "-")

    logical <- logical & (matrix$Organ %in% organ)
  }
  else if (is.null(organ)) {
    prefix <- paste0(mouse, collapse = "-")

    logical <- logical & (matrix$Mouse %in% mouse)
  }
  else {
    prefix <- paste0(paste0(mouse, collapse = "-"), "_", paste0(organ, collapse = "-"))
    if (length(prefix) > 1) {
      prefix <- paste0(prefix, collapse = "+")
    }

    logical <- logical & (matrix$Mouse %in% mouse) & (matrix$Organ %in% organ)
  }

  print(logical)

  mids <- rownames(matrix[which(logical), ])
  data <- exp[, mids]
  labels <- matrix[mids, "joined"]
  colours <- matrix[mids, "colours"]

  print(mids)
  print(data)
  print(matrix[mids, ])
  print(labels)

  if (heatmap) {
    make_full_heatmap(data, method = "horn", labels = labels, dirname = dirname, prefix = prefix, types = "nonsquare", scale = "columnperc")
    make_full_heatmap(data, method = "horn", labels = labels, dirname = dirname, prefix = prefix, types = "square")
  }

  # make_chord_diag(exp[, mids], dirname=paste0(dirname, "/chord"), prefix=prefix)

  if (radarchart) {
    save_radarchart(data, mids,
      labels = labels, dirname = paste0(dirname, "/radar"),
      filename = prefix, percentage = T, use_log = T,
      vlabels = rownames(data), sort_by_sum = F,
      use_clones = rownames(data), colours = colours
    )
  }
}

make_chord_diag <- function(data, dirname, prefix, topPerc = NULL) {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  data_copy <- data

  if (!is.null(topPerc)) {
    rows <- c()

    # for each column, take the topPerc clones (non-null rows)
    # and merge them together.
    for (column in colnames(data_copy))
    {
      df <- data_copy[column]
      df <- subset(df, df > 0)
      df <- df[order(df[, c(1)], decreasing = T), c(1), drop = F]
      abundance <- sum(df) * topPerc
      df <- df[cumsum(df) <= abundance, , drop = F]

      # fix for issue identified on 29/05/2017:
      # if we take the union, some of the clones left out will
      # essentially "come back" for each organ.
      data_copy[setdiff(rownames(data_copy), rownames(df)), column] <- 0

      rows <- union(rows, rownames(df))
    }
  }
  else {
    rows <- rownames(data_copy)
  }

  data_copy <- data_copy[rows, ]

  print(dim(data_copy))

  d <- sweep(data_copy, 2, colSums(data_copy), "/")
  rownames(d) <- 1:nrow(d)
  cairo_pdf(paste0(dirname, "/", prefix, "_chord.pdf"))
  chordDiagram(as.matrix(d))
  dev.off()
}

abundance_matrix <- function(df, topPerc = NULL) {
  # get the total number of reads by MID
  df_sums <- colSums(df)

  # extract the top % of clones and
  # the modified data (eg, for each organ the discarded
  # clones' abudance is set to zero)
  data <- get_top_clones(df, topPerc = topPerc)

  n <- ncol(df)
  mat <- matrix(0, ncol = n, nrow = n)
  colnames(mat) <- rownames(mat) <- colnames(df)

  for (i in 1:nrow(mat))
  {
    for (j in 1:ncol(mat))
    {
      # get the clones in the intersection, using the corrected data
      # (eg, by artificially setting to zero the excluded clones, by MID)
      clones <- which(data[, i] != 0 & data[, j] != 0)

      # now we need the sum of the abundance of these clones
      # divided by the total for this MID
      rel_abund <- sum(data[clones, j]) / df_sums[j]

      mat[i, j] <- rel_abund
    }
  }

  return(mat)
}

make_abundance_matrix <- function(df, topPerc = NULL, dirname = "dissimilarity", prefix = "full") {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  write.csv(abundance_matrix(df, topPerc = topPerc), file = paste0(dirname, "/", prefix, "_abund-inters.csv"))
}

logseq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}

get_freqs <- function(vec) {
  vec / sum(vec)
}

get_present_clones <- function(column) {
  column <- subset(column, column > 0)

  return(rownames(column))
}

get_top_clones <- function(data, topN = NULL, topPerc = NULL, min_clones = 3, exclude = NULL) {
  columns <- setdiff(colnames(data), exclude)

  if (!is.null(topN)) {
    rows <- c()

    topN <- max(topN, min_clones)

    # for each column, take the topN clones (non-null rows)
    # and merge them together.
    for (column in columns)
    {
      df <- data[column]
      df <- subset(df, df > 0)
      df <- df[order(df[, c(1)], decreasing = T), c(1), drop = F]
      df <- df[1:min(nrow(df), topN), , drop = F]

      rows <- union(rows, rownames(df))
    }
  }
  else if (!is.null(topPerc)) {
    rows <- c()

    # for each column, take the topPerc clones (non-null rows)
    # and merge them together.
    for (column in columns)
    {
      df <- data[column]
      df <- subset(df, df > 0)
      df <- df[order(df[, c(1)], decreasing = T), c(1), drop = F]
      abundance <- sum(df) * topPerc
      df_rows <- rownames(df)[cumsum(df) <= abundance]

      # we force a minimum number of clones
      n_clones <- max(length(df_rows), min_clones)
      df_rows <- rownames(df)[1:n_clones]

      # fix for issue identified on 29/05/2017:
      # when we take the union (couple of lines below),
      # some of the clones left out will "come back" for each organ (but they are NOT among the top for that organ!).
      # so we set them to zero column by column, to make sure that the non-selected clones completely disappear.
      data[setdiff(rownames(data), df_rows), column] <- 0

      rows <- union(rows, df_rows)
    }
  }
  else {
    rows <- rownames(data)
  }

  return(data[rows, ])
}

make_full_heatmap <- function(data, method, labels = NULL, binary = F, dirname = "heatmaps",
                              prefix = "full", clustering = F, topN = NULL, topPerc = NULL, cexRow = 0.5,
                              cexCol = 0.5, types = c("square", "nonsquare"), vlim = NULL, scale = "none", tsne = F) {
  # the returned "data" must overwrite the parameter "data",
  # because we modified some of the cells (see the function body)
  data <- get_top_clones(data, topN = topN, topPerc = topPerc)

  if (!is.null(labels)) {
    colnames(data) <- labels
  }

  print(head(data))
  print(dim(data))

  if ("square" %in% types) {
    k <- ifelse(clustering, 4, 0)
    l <- make_heatmap(data,
      dirname = dirname, prefix = paste0(prefix, "_", method, "_heatmap-square"),
      threshold = 0, k = k, dist_method = method, abs = F, method = "ward.D2", scale = scale,
      square = T, output = cairo_pdf, cexRow = cexRow, cexCol = cexRow, vlim = vlim
    )

    if (tsne) {
      set.seed(42)
      cairo_pdf(paste0(dirname, "/", prefix, "_", method, "_tsne-square.pdf"))
      dist <- vegdist(t(data), method = "horn", diag = T, upper = T, binary = F)
      res <- Rtsne(dist, perplexity = ceiling(ncol(data) / 6), theta = 0.0, max_iter = 5000, eta = 1, is_distance = T)
      plot(res$Y, cex = 0.1, xlim = c(min(res$Y[, 1]) - 0.1, max(res$Y[, 1]) + 15), main = "Tsne (with MHI distance)")
      text(res$Y, colnames(data), cex = 0.5, pos = 4)
      dev.off()
    }
  }

  if ("nonsquare" %in% types) {
    k <- ifelse(clustering, 4, 0)
    k2 <- ifelse(clustering, 4, 0)
    l <- make_heatmap(data,
      dirname = dirname, prefix = paste0(prefix, "_", method, "_heatmap"),
      threshold = 0, k = k, k2 = k2, dist_method = method, abs = F, method = "ward.D2", scale = scale,
      square = F, output = cairo_pdf, cexRow = cexRow, cexCol = cexRow, vlim = vlim
    )
  }

  write.csv(l$data, file = paste0(dirname, "/", prefix, "_", method, "_matrix.csv"))
}

# ## alternative clustering for heatmap
# d=hclust(vegdist(data, method=method, diag=T, upper=F, binary=binary), method="average")
# heatmap.2(1-as.matrix(vegdist(data, method = method, diag = T, upper = T, binary = binary)), margins=c(7,7), dendrogram="row", Rowv=as.dendrogram(d), density.info="none", trace="none", cexRow=0.6, cexCol=0.6, col=colours, breaks=breaks, lhei = c(1, 5))

make_input_histograms <- function(data, common, uncommon, columns, dirname = "distributions", filename = "input") {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  svg(paste0(dirname, "/", filename, ".svg"))

  layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))

  plot(sort(rowSums(data[common, columns]), decreasing = T), type = "l", xlab = "clones", ylab = "# reads", main = "Common clones (input)", col = "blue")
  plot(sort(rowSums(data[uncommon, columns]), decreasing = T), type = "l", xlab = "clones", ylab = "# reads", main = "Non-common clones (input)", col = "red")

  plot(sort(rowSums(data[common, columns]), decreasing = T), type = "l", log = "xy", xlab = "clones", ylab = "# reads", main = "Common (blue) and Uncommon (red) clones, log scale", col = "blue")
  lines(sort(rowSums(data[uncommon, columns]), decreasing = T), type = "l", col = "red")

  dev.off()
}

make_rarefaction_plots <- function(data, input_columns = NULL, dirname = "distributions", suffix = "", rarefy = T) {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  if (!is.null(input_columns)) {
    if (rarefy) {
      min_count <- min(colSums(data[, input_columns]))
      print(rarefy(t(data[, input_columns]), sample = min_count, se = T))
      df <- rrarefy(t(data[, input_columns]), sample = min_count)
      svg(paste0(dirname, "/", "rarefied_inputs", suffix, ".svg"))
      rarecurve(df, ylab = "Clones", xlab = "# Reads", col = rainbow(3))
      dev.off()
      svg(paste0(dirname, "/", "rarefied_inputs_log", suffix, ".svg"))
      rarecurve(df, ylab = "Clones", xlab = "# Reads", col = rainbow(3), log = "xy")
      dev.off()
    }

    svg(paste0(dirname, "/", "unrarefied_inputs", suffix, ".svg"))
    rarecurve(t(data[, input_columns]), step = 10, ylab = "Clones", xlab = "# Reads", col = rainbow(3))
    dev.off()
    svg(paste0(dirname, "/", "unrarefied_inputs_log", suffix, ".svg"))
    rarecurve(t(data[, input_columns]), step = 10, ylab = "Clones", xlab = "# Reads", col = rainbow(3), log = "xy")
    dev.off()
  }

  svg(paste0(dirname, "/", "unrarefied_all", suffix, ".svg"))
  rarecurve(t(data), step = 10, ylab = "Clones", xlab = "# Reads", col = rainbow(ncol(data)))
  dev.off()

  svg(paste0(dirname, "/", "unrarefied_all_log", suffix, ".svg"))
  rarecurve(t(data), step = 10, ylab = "Clones", xlab = "# Reads", col = rainbow(ncol(data)), log = "xy")
  dev.off()

  if (rarefy) {
    min_count <- min(colSums(data))
    print(rarefy(t(data), sample = min_count, se = T))
    df <- rrarefy(t(data), sample = min_count)
    svg(paste0(dirname, "/", "rarefied_all", suffix, ".svg"))
    rarecurve(df, ylab = "Clones", xlab = "# Reads", col = rainbow(3))
    dev.off()
    svg(paste0(dirname, "/", "rarefied_all_log", suffix, ".svg"))
    rarecurve(df, ylab = "Clones", xlab = "# Reads", col = rainbow(3), log = "xy")
    dev.off()
  }
}

make_venn_2sets <- function(A, B) {
  a <- length(setdiff(A, B))
  b <- length(setdiff(B, A))

  ab <- length(intersect(A, B))

  return(c(A = a, B = b, "A&B" = ab))
}

make_venn_3sets <- function(A, B, C) {
  a <- length(setdiff(A, union(B, C)))
  b <- length(setdiff(B, union(A, C)))
  c <- length(setdiff(C, union(A, B)))

  ab <- length(setdiff(intersect(A, B), C))
  ac <- length(setdiff(intersect(A, C), B))
  bc <- length(setdiff(intersect(B, C), A))

  abc <- length(intersect(intersect(A, B), C))

  return(c(A = a, B = b, C = c, "A&B" = ab, "A&C" = ac, "B&C" = bc, "A&B&C" = abc))
}

make_venn_4sets <- function(A, B, C, D) {
  # venneuler needs disjoint sets

  a <- length(setdiff(A, union(union(B, C), D)))
  b <- length(setdiff(B, union(union(A, C), D)))
  c <- length(setdiff(C, union(union(A, B), D)))
  d <- length(setdiff(D, union(union(A, B), C)))

  ab <- length(setdiff(intersect(A, B), union(C, D)))
  ac <- length(setdiff(intersect(A, C), union(B, D)))
  ad <- length(setdiff(intersect(A, D), union(B, C)))
  bc <- length(setdiff(intersect(B, C), union(A, D)))
  bd <- length(setdiff(intersect(B, D), union(A, C)))
  cd <- length(setdiff(intersect(C, D), union(A, B)))

  abc <- length(setdiff(intersect(intersect(A, B), C), D))
  abd <- length(setdiff(intersect(intersect(A, B), D), C))
  acd <- length(setdiff(intersect(intersect(A, C), D), B))
  bcd <- length(setdiff(intersect(intersect(B, C), D), A))

  abcd <- length(intersect(intersect(intersect(A, B), C), D))

  return(c(A = a, B = b, C = c, D = d, "A&B" = ab, "A&C" = ac, "A&D" = ad, "B&C" = bc, "B&D" = bd, "C&D" = cd, "A&B&C" = abc, "A&B&D" = abd, "A&C&D" = acd, "B&C&D" = bcd, "A&B&C&D" = abcd))
}

save_venn_diagram <- function(venn_sets, labels, filename, dirname = "venn") {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  svg(paste0(dirname, "/", filename, ".svg"))
  vd <- venneuler(venn_sets)
  vd$labels <- labels
  print(vd)
  plot(vd)
  dev.off()
}

save_heatmap <- function(df, col = df$m1, cols = c(1, 2, 3), filename, dirname = "heatmaps_old") {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  svg(paste0(dirname, "/", filename, ".svg"))
  matrix <- as.matrix(df[order(-col)[1:100], cols])
  # matrix = as.matrix(df[order(-col),cols])
  heatmap.2(matrix, col = colorpanel(75, low = "black", high = "green"), labRow = F, Colv = NA, Rowv = NA, cexCol = 0.7, margins = c(5, 5), dendrogram = "none", density.info = "none", trace = "none", hclustfun = hclustfun)
  dev.off()
}

save_piechart <- function(df1, df2, main1, main2, filename, df3 = NULL, main3 = NULL, df4 = NULL, main4 = NULL, common_clones = NULL, show_black = F, dirname = "./piecharts", topN = NULL, topPerc = NULL, order_by_size = T) {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  if (is.null(topN)) {
    topN <- Inf
  }

  df1 <- subset(df1, df1 > 0)
  df1 <- df1[order(df1[, c(1)], decreasing = T), c(1), drop = F]
  # print(head(df1))

  df2 <- subset(df2, df2 > 0)
  df2 <- df2[order(df2[, c(1)], decreasing = T), c(1), drop = F]
  # print(head(df2))

  count <- 2

  common <- intersect(rownames(df1), rownames(df2))
  tot <- union(rownames(df1), rownames(df2))

  if (!is.null(df3)) {
    df3 <- subset(df3, df3 > 0)
    df3 <- df3[order(df3[, c(1)], decreasing = T), c(1), drop = F]
    # print(head(df3))

    common <- intersect(common, rownames(df3))
    tot <- union(tot, rownames(df3))

    count <- count + 1
  }

  if (!is.null(df4)) {
    df4 <- subset(df4, df4 > 0)
    df4 <- df4[order(df4[, c(1)], decreasing = T), c(1), drop = F]
    # print(head(df4))

    common <- intersect(common, rownames(df4))
    tot <- union(tot, rownames(df4))

    count <- count + 1
  }

  if (!is.null(common_clones)) {
    df1_rows <- intersect(rownames(df1), rownames(common_clones))

    if (is.null(topPerc)) {
      df1_rows <- df1_rows[1:min(length(df1_rows), topN)]
    }
    else {
      abundance <- sum(df1[df1_rows, ]) * topPerc
      print(abundance)
      print(sum(df1[df1_rows, ]))
      print(head(df1[df1_rows, ]))
      print(cumsum(df1[df1_rows, ]))
      df1_rows <- df1_rows[cumsum(df1[df1_rows, ]) <= abundance]
      print(df1_rows)
    }

    if (!order_by_size) {
      # now that we selected the N most abundant clones,
      # we re-order them as they appear in common_clones (using intersect as a dirty trick here)
      df1_rows <- intersect(rownames(common_clones), df1_rows)
    }

    df2_rows <- intersect(rownames(df2), rownames(common_clones))

    if (is.null(topPerc)) {
      df2_rows <- df2_rows[1:min(length(df2_rows), topN)]
    }
    else {
      abundance <- sum(df2[df2_rows, ]) * topPerc
      df2_rows <- df2_rows[cumsum(df2[df2_rows, ]) <= abundance]
      print(df2_rows)
    }

    if (!order_by_size) {
      # now that we selected the N most abundant clones,
      # we re-order them as they appear in common_clones (using intersect as a dirty trick here)
      df2_rows <- intersect(rownames(common_clones), df2_rows)
    }

    df1_cols <- as.character(common_clones[df1_rows, c("colours")])
    df2_cols <- as.character(common_clones[df2_rows, c("colours")])

    if (show_black) {
      df1_cols <- c(df1_cols, "black")
      df2_cols <- c(df2_cols, "black")
    }

    png(paste0(dirname, "/", filename, ".png"), height = 400, width = 800, pointsize = 14, bg = "white")
    # svg(paste0("piecharts/", filename, ".svg"))
    par(mfrow = c(1, count))

    if (show_black) {
      pie(c(df1[df1_rows, c(1)], sum(df1[setdiff(rownames(df1), df1_rows), c(1)])), col = df1_cols, labels = NA, density = NA, clockwise = T, main = main1)
      pie(c(df2[df2_rows, c(1)], sum(df2[setdiff(rownames(df2), df2_rows), c(1)])), col = df2_cols, labels = NA, density = NA, clockwise = T, main = main2)
    }
    else {
      pie(df1[df1_rows, c(1)], col = df1_cols, labels = NA, density = NA, clockwise = T, main = main1)
      pie(df2[df2_rows, c(1)], col = df2_cols, labels = NA, density = NA, clockwise = T, main = main2)
    }

    if (!is.null(df3)) {
      df3_rows <- intersect(rownames(df3), rownames(common_clones))

      if (is.null(topPerc)) {
        df3_rows <- df3_rows[1:min(length(df3_rows), topN)]
      }
      else {
        abundance <- sum(df3[df3_rows, ]) * topPerc
        df3_rows <- df3_rows[cumsum(df3[df3_rows, ]) <= abundance]
      }

      if (!order_by_size) {
        # now that we selected the N most abundant clones,
        # we re-order them as they appear in common_clones (using intersect as a dirty trick here)
        df3_rows <- intersect(rownames(common_clones), df3_rows)
      }

      if (show_black) {
        df3_cols <- c(as.character(common_clones[df3_rows, c("colours")]), "black")

        pie(c(df3[df3_rows, c(1)], sum(df3[setdiff(rownames(df3), df3_rows), c(1)])), col = df3_cols, labels = NA, density = NA, clockwise = T, main = main3)
      }
      else {
        df3_cols <- as.character(common_clones[df3_rows, c("colours")])

        pie(df3[df3_rows, c(1)], col = df3_cols, labels = NA, density = NA, clockwise = T, main = main3)
      }
    }

    if (!is.null(df4)) {
      df4_rows <- intersect(rownames(df4), rownames(common_clones))

      if (is.null(topPerc)) {
        df4_rows <- df4_rows[1:min(length(df4_rows), topN)]
      }
      else {
        abundance <- sum(df4[df4_rows, ]) * topPerc
        df4_rows <- df4_rows[cumsum(df4[df4_rows, ]) <= abundance]
      }

      if (!order_by_size) {
        # now that we selected the N most abundant clones,
        # we re-order them as they appear in common_clones (using intersect as a dirty trick here)
        df4_rows <- intersect(rownames(common_clones), df4_rows)
      }

      if (show_black) {
        df4_cols <- c(as.character(common_clones[df4_rows, c("colours")]), "black")

        pie(c(df4[df4_rows, c(1)], sum(df4[setdiff(rownames(df4), df4_rows), c(1)])), col = df4_cols, labels = NA, density = NA, clockwise = T, main = main4)
      }
      else {
        df4_cols <- as.character(common_clones[df4_rows, c("colours")])

        pie(df4[df4_rows, c(1)], col = df4_cols, labels = NA, density = NA, clockwise = T, main = main4)
      }
    }

    dev.off()
  }
  else {
    # if common_clones is not provided, than we colour "black" those clones that are not in common
    # among df1, df2 and df3.

    colours <- c(rainbow(length(common)), "black")

    png(paste0(dirname, "/", filename, ".png"), height = 400, width = 800, pointsize = 14, bg = "white")
    # svg(paste0("piecharts/", filename, ".svg"))
    par(mfrow = c(1, count))

    pie(c(df1[common, c(1)], sum(df1[setdiff(rownames(df1), common), c(1)])), col = colours, labels = NA, density = NA, clockwise = T, main = main1)
    pie(c(df2[common, c(1)], sum(df2[setdiff(rownames(df2), common), c(1)])), col = colours, labels = NA, density = NA, clockwise = T, main = main2)
    if (!is.null(df3)) {
      pie(c(df3[common, c(1)], sum(df3[setdiff(rownames(df3), common), c(1)])), col = colours, labels = NA, density = NA, clockwise = T, main = main3)
    }
    if (!is.null(df4)) {
      pie(c(df4[common, c(1)], sum(df4[setdiff(rownames(df4), common), c(1)])), col = colours, labels = NA, density = NA, clockwise = T, main = main4)
    }

    dev.off()
  }
}

make_radarchart <- function(data, column_names, labels = NULL, topN = NULL, topPerc = NULL, minmax = "global", percentage = T, use_clones = NULL, show_axis = T,
                            use_union = T, sort_by_sum = T, legend_x = "topright", legend_inset = c(0, 0), legend_counts = F, add_clones = NULL, noplot = F, use_log = F,
                            colours_border = c(rgb(0.2, 0.5, 0.6, 0.8), rgb(1.0, 0.2, 0, 0.8), rgb(0.9, 0.7, 0.1, 0.8), rgb(0.2, 0.6, 0.3, 0.8)),
                            colours_in = c(rgb(0.2, 0.5, 0.6, 0.4), rgb(1.0, 0.2, 0, 0.4), rgb(0.9, 0.7, 0.1, 0.4), rgb(0.2, 0.6, 0.3, 0.4)), vlabels = NA) {
  if (is.null(labels)) {
    labels <- column_names
  }

  if (is.null(topN)) {
    topN <- Inf
  }

  data <- data[column_names]

  # backup for later
  orig_data <- data

  if (percentage) {
    # takes the sum of all reads for each organ,
    # then divides the single clonal abudance by the total (again by organ).
    # the output presents relative abundance values (hence, each organ sums up to 1)
    data_sums <- colSums(data)
    data <- sweep(data, 2, data_sums, "/")
  }

  common_rows <- c()

  if (is.null(use_clones)) {
    # since clones aren't provided, we use either the union or the intersection
    # of all our organs' clones
    for (column in column_names)
    {
      df <- subset(data[column], data[column] > 0)
      df <- df[order(df[, , drop = F], decreasing = T), , drop = F]

      if (is.null(topPerc)) {
        df_rows <- rownames(df[1:min(nrow(df), topN), , drop = F])
      }
      else {
        print(sum(df))
        abundance <- sum(df) * topPerc
        print(abundance)
        df_rows <- rownames(df)[cumsum(df[df_rows, ]) <= abundance]
        print(df[df_rows, ])
      }

      if (use_union) {
        common_rows <- union(common_rows, df_rows)
      }
      else {
        common_rows <- intersect(common_rows, df_rows)
      }
    }
  }
  else {
    # clones are provided, so we use them
    # note: topN is ignored

    common_rows <- use_clones
  }

  print(paste0("#clones: ", length(common_rows)))

  if (!is.null(add_clones)) {
    # trick: putting add_clones first compels the ordering to be based on add_clones
    common_rows <- union(add_clones, common_rows)
  }

  # order selected clones by total abundance
  if (sort_by_sum) {
    # for each of the selected clones, take the sum across organs
    clone_abund <- rowSums(data[common_rows, ])
    # sort the clones by their row-sum, bigger-to-smaller
    clone_abund <- sort(clone_abund, decreasing = T)

    # trick: intersect fixes the order using the first argument
    common_rows <- intersect(names(clone_abund), common_rows)
  }

  data <- as.data.frame(t(data[common_rows, ]))
  colnames(data) <- 1:ncol(data)

  # if counts are needed in the legend, we must go back
  # to the original data (ie, before getting percentanges)
  orig_data <- as.data.frame(t(orig_data[common_rows, ]))
  colnames(orig_data) <- 1:ncol(orig_data)
  organ_abund_counts <- rowSums(orig_data)
  organ_abund_perc <- rowSums(data)

  print("organ abundance:")
  print(organ_abund_counts)
  print(organ_abund_perc)

  if (minmax == "global") {
    # global minimum and maximum: this makes it easier to see
    # print(data)
    print(max(data))
    print(min(data))
    data <- rbind(max(data), min(data), data)
  }
  else if (minmax == "by-column") {
    # by-organ max/min: might be hard to visualise
    data <- rbind(colMax(data), colMin(data), data)
  }
  else if (minmax == "by-row") {
    # by-organ max/min: might be hard to visualise
    data <- rbind(rowMax(data), rowMin(data), data)
  }
  else if (length(minmax) == 2) {
    # [min, max]

    data <- rbind(rep(minmax[2], ncol(data)), rep(minmax[1], ncol(data)), data)
  }
  else {
    stop("ERROR: minmax argument wrong")
  }

  print(data)

  if (!noplot) {
    if (use_log) {
      if (percentage) {
        data <- data * 100
      }

      # zero becomes one, so when log'ed it is still zero
      data <- log(data + 1)
    }

    print(data)

    if (minmax == "global" || length(minmax == 2)) {
      center_labels <- seq(data[2, 1], data[1, 1], length.out = 5)

      # un-log it
      if (use_log) {
        center_labels <- sprintf("%.5f", exp(center_labels) - 1)
      }
      else {
        center_labels <- sprintf("%.5f", center_labels)
      }

      # add the percent symbol
      if (percentage) {
        center_labels <- paste0(center_labels, "%")
      }
    }
    else {
      center_labels <- NULL
    }

    print(center_labels)

    title <- paste("#clones:", length(common_rows))

    if (show_axis) {
      axistype <- 1
    }
    else {
      axistype <- 0

      title <- paste(title, "| axis:", paste(center_labels, collapse = " "))
    }

    # radarchart(data, axistype=axistype, pcol=colours_border , pfcol=colours_in , plwd=4 , plty=1, cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, vlcex=0.8)
    radarchart(data,
      axistype = axistype, pcol = colours_border, pfcol = colours_in, plwd = 2, plty = 1, caxislabels = center_labels,
      cglcol = "grey", cglty = 1, axislabcol = "darkslategrey", cglwd = 0.8, vlcex = 0.8, maxmin = T, vlabels = vlabels, title = title
    )

    if (!is.null(legend_x)) {
      if (legend_counts || !percentage) {
        legend_values <- paste0(labels, " (", round(organ_abund_counts, digits = 3), ")")
      }
      else {
        legend_values <- paste0(labels, " (", round(organ_abund_perc * 100, digits = 3), "%)")
      }

      legend(legend_x, inset = legend_inset, legend = legend_values, bty = "n", pch = 20, col = colours_in, text.col = "grey", cex = 0.85, pt.cex = 2)
    }
  }

  return(list(
    clones = common_rows,
    data = data,
    organ_abund_counts = organ_abund_counts,
    organ_abund_perc = organ_abund_perc
  ))
}

save_radarchart <- function(data, column_names, labels = NULL, topN = NULL, dirname = "./radarcharts", filename = "radar", minmax = "global", percentage = T, sort_by_sum = T, use_log = F, show_axis = T, legend_x = "topright", vlabels = NA, use_clones = NULL, colours = NULL) {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  if (!is.null(colours)) {
    colours_border <- add.alpha(colours, 0.8)
    colours_in <- add.alpha(colours, 0.08)
  }
  else {
    # FIXME: this default is pretty crappy
    colours_border <- c(rgb(0.2, 0.5, 0.6, 0.8), rgb(1.0, 0.2, 0, 0.8), rgb(0.9, 0.7, 0.1, 0.8), rgb(0.2, 0.6, 0.3, 0.8))
    colours_in <- c(rgb(0.2, 0.5, 0.6, 0.4), rgb(1.0, 0.2, 0, 0.4), rgb(0.9, 0.7, 0.1, 0.4), rgb(0.2, 0.6, 0.3, 0.4))
  }

  svg(paste0(dirname, "/", filename, ".svg"))

  make_radarchart(data, column_names, labels = labels, topN = topN, minmax = minmax, percentage = percentage, sort_by_sum = sort_by_sum, use_log = use_log, show_axis = show_axis, legend_x = legend_x, vlabels = vlabels, use_clones = use_clones, colours_border = colours_border, colours_in = colours_in)

  dev.off()
}

save_double_radarchart <- function(data, column_names, data2, column_names2, labels = NULL, labels2 = NULL, topN = NULL, dirname = "./radarcharts/multiple", filename = "radar", minmax = "global", percentage = T, sort_by_sum = T, use_log = F) {
  if (!dir.exists(dirname)) {
    dir.create(dirname, mode = "0755", recursive = T)
  }

  svg(paste0(dirname, "/", filename, ".svg"))

  par(xpd = TRUE, mar = c(0, 1, 0, 2), mfrow = c(1, 2))

  ## ugly trick coming.. brace yourself

  # first we calculate the clones for data2, without actually plotting anything
  clones <- make_radarchart(data2, column_names2, labels = labels2, topN = topN, minmax = minmax, percentage = percentage, sort_by_sum = sort_by_sum, noplot = T)

  # then we pass the above-calculated clones to create a "joined" radarchart for data1, that we plot...
  clones <- make_radarchart(data, column_names, labels = labels, topN = topN, minmax = minmax, percentage = percentage, sort_by_sum = sort_by_sum, legend_x = "top", legend_inset = c(0, 0.15), add_clones = clones, use_log = use_log)

  # ...then we pass the joined clones (relevant for both data1 and data2) to the radarchart for data2, that we plot (right-hand side).
  # of course, we cannot sort the clones by sum here, or we screw up the matching between data1 and data2. So we are stuck with whatever
  # ordering was found for data1.
  clones2 <- make_radarchart(data2, column_names2, labels = labels2, topN = topN, minmax = minmax, percentage = percentage, sort_by_sum = F, legend_x = "top", legend_inset = c(0, 0.15), add_clones = clones, use_log = use_log)

  # and we write down the image
  dev.off()

  # now, clones should be equal to clones2 - otherwise, something's gone amiss
  stopifnot(identical(clones, clones2))
}

## Add an alpha value to a colour
add.alpha <- function(col, alpha = 1) {
  if (missing(col)) {
    stop("Please provide a vector of colours.")
  }

  apply(sapply(col, col2rgb) / 255, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

save_multiple_radarcharts <- function(data, column_names, labels = NULL, topN = NULL, topPerc = NULL, dirname = "./radarcharts/one-vs-all", minmax = "global", percentage = T, use_log = F, use_colours = NULL, join_plots = F, show_legend = T, show_axis = T, legend_counts = F) {
  print("save_multiple_radarcharts:")
  print(column_names)

  dirname <- make_path(dirname)

  if (is.null(labels)) {
    labels <- column_names
  }

  if (is.null(topN)) {
    topN <- Inf
  }

  if (is.null(use_colours) || length(use_colours) < length(column_names)) {
    colours <- distinctColorPalette(length(column_names))
  }
  else {
    colours <- use_colours
  }

  colours_border <- add.alpha(colours, 0.8)
  colours_in <- add.alpha(colours, 0.08)

  if (join_plots) {
    svg(paste0(dirname, "all.svg"))

    par(xpd = TRUE, mar = c(0, 1, 0, 2), mfcol = c(ceiling(length(column_names) / 2), 2))
  }

  all_counts <- data.frame()

  ret <- list()

  # start with rownames to keep the ordering, then subset with "intersect"
  common_clones <- rownames(data)

  # we accumulate the selected clones of each radarchart, then at the end
  # we intersect with rownames to keep the ordering
  total_clones <- c()

  for (i in 1:length(column_names)) {
    column <- column_names[i]
    label <- labels[i]

    print(label)

    if (!join_plots) {
      svg(paste0(dirname, label, ".svg"))
    }

    df <- subset(data[column], data[column] > 0)
    df <- df[order(df[, , drop = F], decreasing = T), , drop = F]

    if (is.null(topPerc)) {
      df_rows <- rownames(df[1:min(nrow(df), topN), , drop = F])
    }
    else {
      abundance <- sum(df) * topPerc
      df_rows <- rownames(df)[cumsum(df) <= abundance]
    }

    if (length(df_rows) < 3) {
      print(paste0("ERROR: ", column, " has ", length(df_rows), " selected, but needs at least 3. Skipping."))
      next()
    }

    if (show_legend) {
      # each radarchart consists of the i-th organ topN clones sorted by abudance,
      # with the other organs following through
      l <- make_radarchart(data, column_names,
        labels = labels, use_clones = df_rows, sort_by_sum = F, minmax = minmax, percentage = percentage,
        use_log = use_log, colours_border = colours_border, colours_in = colours_in, show_axis = show_axis, legend_counts = legend_counts
      )
    }
    else {
      # each radarchart consists of the i-th organ topN clones sorted by abudance,
      # with the other organs following through
      l <- make_radarchart(data, column_names,
        labels = labels, use_clones = df_rows, sort_by_sum = F, minmax = minmax, percentage = percentage,
        use_log = use_log, colours_border = colours_border, colours_in = colours_in, legend_x = NULL, show_axis = show_axis
      )
    }

    organ_abund <- data.frame(
      Index = rep(label, length(labels)),
      Sample = labels,
      Counts = l$organ_abund_counts,
      Percentage = l$organ_abund_perc
    )
    write.csv(organ_abund, file = paste0(dirname, "/", label, ".csv"), row.names = F)

    all_counts <- rbind(all_counts, organ_abund)

    if (!join_plots) {
      dev.off()
    }

    common_clones <- intersect(common_clones, l$clones)
    total_clones <- union(total_clones, l$clones)
  }

  write.csv(all_counts, file = paste0(dirname, "all.csv"), row.names = F)

  if (join_plots) {
    dev.off()
  }

  print(paste0("common clones: ", length(common_clones)))
  print(paste0("total clones: ", length(total_clones)))

  ret$common_clones <- common_clones
  ret$total_clones <- intersect(rownames(data), total_clones)

  return(ret)
}
