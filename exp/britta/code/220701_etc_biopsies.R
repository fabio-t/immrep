
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

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("m12", "m12_")
heatmaps("m13", "m13_")
heatmaps("m14", "m14_")
heatmaps("m53", "m53_")
heatmaps("m32", "m32_")
heatmaps("m42", "m42_")
heatmaps("m43", "m43_")
heatmaps("m30", "m30_")

track_clones(NULL, NULL, immdata = idata)
track_clones("m12", "m12_", immdata = idata)
track_clones("m13", "m13_", immdata = idata)
track_clones("m14", "m14_", immdata = idata)
track_clones("m53", "m53_", immdata = idata)
track_clones("m32", "m32_", immdata = idata)
track_clones("m42", "m42_", immdata = idata)
track_clones("m43", "m43_", immdata = idata)
track_clones("m30", "m30_", immdata = idata)


### TODO: join the biopsy experiments
### TODO: add d50 to diversity, and also fix it / improve it (no quotes, 1-evenness)
