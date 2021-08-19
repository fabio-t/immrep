
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

## radarcharts

for (tp in c(1, 0.85)) {
  for (p_i in c(1, 5, 10, 24, 25, 29, 44, 53, 58, 59)) {
    radarcharts(paste0("^p", p_i, "_"), paste0("p", p_i), save=T)
  }
  radarcharts("^p(24|10|5)_", "p24_p10_p5", save=T)
}

track_clones(NULL, NULL, immdata=idata)
