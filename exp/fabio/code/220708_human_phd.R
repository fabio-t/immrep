
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

track_clones(NULL, NULL, immdata = idata)
# cross-organ
track_clones(NULL, "c1", indices = which(idata$meta$cross == "c1"), immdata = idata)
track_clones(NULL, "c2", indices = which(idata$meta$cross == "c2"), immdata = idata)
track_clones(NULL, "c3", indices = which(idata$meta$cross == "c3"), immdata = idata)
track_clones(NULL, "c4", indices = which(idata$meta$cross == "c4"), immdata = idata)
track_clones(NULL, "c5", indices = which(idata$meta$cross == "c5"), immdata = idata)

# by organ
track_clones(NULL, "liver", indices = which(idata$meta$organ == "liver"), immdata = idata)
track_clones(NULL, "si", indices = which(idata$meta$organ == "si"), immdata = idata)

# by patient
track_clones(NULL, "p10", indices = which(idata$meta$patient == "p10"), immdata = idata)
track_clones(NULL, "p11", indices = which(idata$meta$patient == "p11"), immdata = idata)
track_clones(NULL, "p7", indices = which(idata$meta$patient == "p7"), immdata = idata)
track_clones(NULL, "p8", indices = which(idata$meta$patient == "p8"), immdata = idata)
track_clones(NULL, "p9", indices = which(idata$meta$patient == "p9"), immdata = idata)
