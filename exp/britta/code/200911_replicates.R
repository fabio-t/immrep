
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)

idata <- immload(which="not full")

mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)
mid_names <- mid_labels[colnames(mids_counts), "label", drop=T]

overview(immdata=idata)
make_rarefaction_plots(mids_counts, rarefy=T)
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("rep7", "no_m2rep7", cex = 0.5, invert = T)
heatmaps("m1", "m1")
heatmaps("m2", "m2")
heatmaps("m2 rep[^7]", "m2_norep7")

## radarcharts

for (tp in c(1, 0.85)) {
  radarcharts("m1", "m1",  save=T)
  radarcharts("m2", "m2",  save=T)
  radarcharts("m2 rep[^7]", "m2_norep7", save=T)
  radarcharts("rep7", "m1m2_norep7", save=T, invert=T)
  radarcharts("m[12]", "m1-m2",  save=T)

  radarcharts("(m1|m2 rep[37])", "m1_plus_m2-rep3and7",  save=T)
}

track_clones(NULL, NULL, immdata=idata)
track_clones("rep7", "no_m2rep7", invert = T, immdata=idata)
track_clones("m1", "m1", immdata=idata)
track_clones("m2", "m2", immdata=idata)
track_clones("m2 rep[^7]", "m2_norep7", immdata=idata)
track_clones("rep7", "m1m2_norep7", invert=T, immdata=idata)
track_clones("(m1|m2 rep[37])", "m1_plus_m2-rep3and7", immdata=idata)
