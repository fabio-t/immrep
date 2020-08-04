
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = colnames(mids_counts) # should be substituted with proper labels, but were not provided

# make_rarefaction_plots(mids_counts, rarefy=T)

## diversity

make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)

## radarcharts

for (tp in c(1, 0.85)) {
    radarcharts("MID[234]$", "mid2-3-4", save = T)
}
