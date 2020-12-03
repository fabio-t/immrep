
library(xlsx)
library(circlize)
library(gtools)

cols <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# must be run inside an analysis directory

if (!dir.exists("chords")) {
  dir.create("chords", mode = "0755", recursive = T)
}

if (!dir.exists("barplots")) {
  dir.create("barplots", mode = "0755", recursive = T)
}

trav_names <- c()
traj_names <- c()
df_list <- list()
df_list2 <- list()

for (f in Sys.glob("mid*clones.csv"))
{
  name <- gsub(".csv$", "", basename(f))

  d <- read.csv(f, sep = "\t")
  d$bestVFamily <- sapply(strsplit(as.character(d$V.segments), "[,\\*-]"), `[`, 1)
  d$bestJFamily <- sapply(strsplit(as.character(d$J.segments), "[,\\*-]"), `[`, 1)
  d$V.segments <- NULL
  d$J.segments <- NULL

  d$bestVFamily <- gsub("[DN]?(/DV.+)?$", "", d$bestVFamily)
  d$bestVFamily <- gsub("/OR.*$", "", d$bestVFamily)

  families <- d$bestVFamily
  table_f <- as.data.frame(table(families))
  # table_f$num = as.numeric(gsub("(TR[AB]|IGH)V","",table_f$families))
  # table_f = table_f[order(table_f$num), ]
  table_f <- table_f[mixedorder(table_f$families), ]
  print(table_f)

  svg(paste0("barplots/", name, ".svg"))
  par(mar = c(4, 8, 2, 2))
  barplot(table_f$Freq, names.arg = table_f$families, horiz = T, las = 2, cex.names = 0.7, col = cols[table_f$num])
  dev.off()

  # V
  df <- data.frame(aggregate(Percentage ~ bestVFamily, data = d, sum), row.names = 1)
  # df$num = as.numeric(gsub("(TR[AB]|IGH)V", "", rownames(df)))
  # df = df[order(df$num), ]
  df <- df[mixedorder(rownames(df)), , drop = F]
  print(df)

  svg(paste0("barplots/", name, "_abund.svg"))
  par(mar = c(4, 8, 2, 2))
  barplot(df$Percentage, names.arg = rownames(df), horiz = T, las = 2, cex.names = 0.7, col = cols[df$num])
  dev.off()

  # J
  df2 <- data.frame(aggregate(Percentage ~ bestJFamily, data = d, sum), row.names = 1)
  # df2$num = as.numeric(gsub("(TR[AB]|IGH)J", "", rownames(df2)))
  # df2 = df2[order(df2$num), ]
  df2 <- df2[mixedorder(rownames(df2)), , drop = F]
  print(df2)

  trav_names <- c(trav_names, rownames(df))
  traj_names <- c(traj_names, rownames(df2))

  df_list[[toupper(gsub("_clones", "", name))]] <- df
  df_list2[[toupper(gsub("_clones", "", name))]] <- df2

  svg(paste0("chords/", name, "_chord.svg"))
  # d$bestVFamily = gsub("(TR[AB]|IGH)", "", d$bestVFamily)
  # d$bestJFamily = gsub("(TR[AB]|IGH)", "", d$bestJFamily)
  tab <- table(d$bestVFamily, d$bestJFamily)
  chorddf <- as.data.frame(tab)
  grid.cols <- c(rep("grey", length(unique(chorddf$Var2))), rainbow(length(unique(chorddf$Var1))))
  labels <- c(names(sort(colSums(tab), decreasing = T)), names(sort(rowSums(tab), decreasing = T)))
  chordDiagram(chorddf, link.zindex = rank(chorddf[[3]]), order = labels, transparency = 0.05, grid.col = grid.cols, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(chorddf))))), link.sort = T, link.decreasing = F, annotationTrack = "grid", annotationTrackHeight = 0.03)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, cex = 0.5, facing = "clockwise", niceFacing = T, adj = c(0, 0.5))
  }, bg.border = NA)
  dev.off()
}

# V

trav_names <- unique(trav_names)
# nums = as.numeric(gsub("(TR[AB]|IGH)V","", trav_names))
# trav_names = trav_names[order(nums)]
trav_names <- mixedsort(trav_names)
df <- data.frame(row.names = trav_names)

for (name in names(df_list))
{
  temp <- df_list[[name]]

  df[rownames(temp), name] <- temp$Percentage
}

df[is.na(df)] <- 0
print(df)

write.xlsx(t(df), file = "table.xlsx", sheetName = "V Family")

# J
traj_names <- unique(traj_names)
# nums = as.numeric(gsub("(TR[AB]|IGH)J","", traj_names))
# traj_names = traj_names[order(nums)]
traj_names <- mixedsort(traj_names)
df <- data.frame(row.names = traj_names)

for (name in names(df_list2))
{
  temp <- df_list2[[name]]

  df[rownames(temp), name] <- temp$Percentage
}

df[is.na(df)] <- 0
print(df)

write.xlsx(t(df), file = "table.xlsx", sheetName = "J Family", append = T)
