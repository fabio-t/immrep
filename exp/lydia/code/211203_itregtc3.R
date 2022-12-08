
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names = 1)
print(mid_labels)

idata <- immload(which = "not full")

mids_counts <- read.csv("mids_counts.csv", sep = "\t", row.names = 1)
mid_names <- mid_labels[colnames(mids_counts), "label", drop = T]

overview(immdata = idata)
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("m173", "m173")
heatmaps("m174", "m174")
heatmaps("input", "input")
heatmaps("pool", "pool")

track_clones(NULL, NULL, immdata = idata)
track_clones("m173", "m173", immdata = idata)
track_clones("m174", "m174", immdata = idata)
track_clones("input", "input", immdata = idata)
track_clones("pool", "pool", immdata = idata)
