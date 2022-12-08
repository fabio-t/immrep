
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)

idata <- immload(which="not full")

mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)
mid_names <- mid_labels[colnames(mids_counts), "label", drop=T]

overview(immdata=idata)
# make_rarefaction_plots(mids_counts, rarefy=T)
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("SPF", "spf")
heatmaps("MM12", "mm12")

for (tp in c(1, 0.85)) {
  radarcharts("m1_", "m1", save=T)
  radarcharts("m2_", "m2", save=T)
  radarcharts("m3_", "m3", save=T)
  radarcharts("SPF", "spf", save=T)
  radarcharts("MM12", "mm12", save=T)
}

track_clones(NULL, NULL, immdata=idata)
track_clones("m1_", "m1", immdata=idata)
track_clones("m2_", "m2", immdata=idata)
track_clones("m3_", "m3", immdata=idata)
track_clones("SPF", "spf", immdata=idata)
track_clones("MM12", "mm12", immdata=idata)
