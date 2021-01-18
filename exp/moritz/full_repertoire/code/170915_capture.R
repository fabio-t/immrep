
library(rprojroot)
root <- is_vcs_root$make_fix_file()
source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

organs = c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")

mid_names = c("M1_input", "M2_input", "M3_input",
              paste0("M1_", organs),
              paste0("M2_", organs),
              paste0("M3_", organs)
          )

# rarefaction
make_rarefaction_plots(mids_counts, 1:3)

# diversity
make_diversity(mids_counts, mid_names)

# heatmaps
heatmaps(NULL, NULL, cex = 0.5, tsne = T)
heatmaps("\\+$", "only-plus")

heatmaps("^M1", "M1")
heatmaps("^M1_.*\\+", "M1_only-plus")
heatmaps("^M2", "M2")
heatmaps("^M2_.*\\+", "M2_only-plus")
heatmaps("^M3", "M3")
heatmaps("^M3_.*\\+", "M3_only-plus")

for (organ in organs) {
  heatmaps(paste0("M\\d+_", gsub("+", "\\+", organ, fixed=T), "$"), organ)
}

# radarcharts
for (tp in c(1, 0.85)) {
  radarcharts("^M1", "M1", save = T)
  radarcharts("^M1_.*\\+", "M1_only-plus", save = T)
  radarcharts("^M2", "M2", save = T)
  radarcharts("^M2_.*\\+", "M2_only-plus", save = T)
  radarcharts("^M3", "M3", save = T)
  radarcharts("^M3_.*\\+", "M3_only-plus", save = T)

  for (organ in organs) {
    radarcharts(paste0("M\\d+_", gsub("+", "\\+", organ, fixed=T), "$"), organ, save = T)
  }
}
