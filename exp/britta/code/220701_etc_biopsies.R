
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names = 1)
print(mid_labels)

idata <- immload(which = "not full")

mids_counts <- read.csv("mids_counts.csv", sep = "\t", row.names = 1)
mid_names <- mid_labels[colnames(mids_counts), "label", drop = T]

# overview(immdata = idata)
make_diversity(mids_counts, mid_names)

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("m11", "m11")
heatmaps("m12", "m12")
heatmaps("m13", "m13")
heatmaps("m14", "m14")
# heatmaps("m28", "m28")
# heatmaps("m29", "m29")
heatmaps("m30", "m30")
heatmaps("m31", "m31")
heatmaps("m32", "m32")
# heatmaps("m36", "m36")
heatmaps("m37", "m37")
# heatmaps("m38", "m38")
heatmaps("m39", "m39")
heatmaps("m42", "m42")
heatmaps("m43", "m43")
heatmaps("m53", "m53")


track_clones(NULL, NULL, immdata = idata)
track_clones("m11", "m11", immdata = idata)
track_clones("m12", "m12", immdata = idata)
track_clones("m13", "m13", immdata = idata)
track_clones("m14", "m14", immdata = idata)
# track_clones("m28", "m28", immdata = idata)
# track_clones("m29", "m29", immdata = idata)
track_clones("m30", "m30", immdata = idata)
track_clones("m31", "m31", immdata = idata)
track_clones("m32", "m32", immdata = idata)
# track_clones("m36", "m36", immdata = idata)
track_clones("m37", "m37", immdata = idata)
# track_clones("m38", "m38", immdata = idata)
track_clones("m39", "m39", immdata = idata)
track_clones("m42", "m42", immdata = idata)
track_clones("m43", "m43", immdata = idata)
track_clones("m53", "m53", immdata = idata)

### TODO: join the biopsy experiments
### TODO: add d50 to diversity, and also fix it / improve it (no quotes, 1-evenness)
