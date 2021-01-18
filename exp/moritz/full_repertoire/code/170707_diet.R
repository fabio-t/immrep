
library(rprojroot)
root <- is_vcs_root$make_fix_file()
source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c("input1", "input2", "input3",

              "M1_IFN-_IL17-", "M1_IFN+",  "M1_IL17+", "M1_FoxP3+",
              "M2_IFN-_IL17-", "M2_IFN+",  "M2_IL17+", "M2_FoxP3+",
              "M3_IFN-_IL17-", "M3_IFN+",  "M3_IL17+", "M3_FoxP3+"
              )

# rarefaction
# make_rarefaction_plots(mids_counts, 1:3)

# diversity
make_diversity(mids_counts, mid_names)

# heatmaps
heatmaps(NULL, NULL, cex = 0.5, tsne = T)
heatmaps("input", "full-noinput", invert = T)
heatmaps("input", "inputs")
heatmaps("input", "inputs")
heatmaps("IFN-_IL17-", "ifn-_il17-")
heatmaps("IFN\\+", "ifn+")
heatmaps("IL17\\+", "il17+")
heatmaps("FoxP3\\+", "foxp3+")
heatmaps("^M1", "M1")
heatmaps("^M2", "M2")
heatmaps("^M3", "M3")

# radarcharts
for (tp in c(1, 0.85)) {
  radarcharts("input", "inputs", save = T)
  radarcharts("IFN-_IL17-", "ifn-_il17-", save = T)
  radarcharts("IFN\\+", "ifn+", save = T)
  radarcharts("IL17\\+", "il17+", save = T)
  radarcharts("FoxP3\\+", "foxp3+", save = T)
  radarcharts("^M1", "M1", save = T)
  radarcharts("^M2", "M2", save = T)
  radarcharts("^M3", "M3", save = T)
}
