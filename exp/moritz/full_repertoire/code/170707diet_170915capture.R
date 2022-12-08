
library(rprojroot)
root <- is_vcs_root$make_fix_file()
source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

organs = c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")

samples = c(
              "M1_input", "M2_input", "M3_input",
              paste0("M1_", organs),
              paste0("M2_", organs),
              paste0("M3_", organs)
          )

mid_names = c(
                paste0("170707_", samples),
                paste0("170915_", samples)
            )

# diversity
make_diversity(mids_counts, mid_names)

# heatmaps
heatmaps(NULL, NULL, cex = 0.5, tsne = T)
heatmaps("input", "full-noinput", invert = T)
heatmaps("input", "inputs")
heatmaps("IFNy-_IL17-", "ifn-_il17-")
heatmaps("IFNy\\+", "ifn+")
heatmaps("IL17\\+", "il17+")
heatmaps("Foxp3\\+", "foxp3+")

# radarcharts
for (tp in c(1, 0.85)) {
  radarcharts("input", "inputs", save = T)
  radarcharts("IFNy-_IL17-", "ifn-_il17-", save = T)
  radarcharts("IFNy\\+", "ifn+", save = T)
  radarcharts("IL17\\+", "il17+", save = T)
  radarcharts("Foxp3\\+", "foxp3+", save = T)
}
