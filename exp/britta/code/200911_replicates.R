
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)
mid_names <- mid_labels[colnames(mids_counts), "label", drop=T]

# make_rarefaction_plots(mids_counts, rarefy=T)

## diversity

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

library(immunarch)
immdata <- repLoad(paste0("full/", rownames(mid_labels), ".csv"))
make_path("tracking")
for (mid in rownames(mid_labels)) {
  tc1 <- trackClonotypes(immdata$data, list(mid, 15), .col = "v+nt")
  vis(tc1) + theme_minimal(base_size=12)
  ggsave(filename = paste0("tracking/", mid, ".pdf"), width = 16, height = (9/16) * 16)
}
