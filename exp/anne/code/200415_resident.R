
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)

idata <- immload(which="full")

mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)
mid_names <- mid_labels[colnames(mids_counts), "label", drop=T]

overview(immdata=idata)
# make_rarefaction_plots(mids_counts, rarefy=T)
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
for (m in 1:3) {
  ms = paste0("m", m)
  heatmaps(ms, ms)
  heatmaps(paste0("Foxp3\\+HD2\\+.*_", ms), paste0(ms, "_plus-plus"))
  # heatmaps(paste0("Foxp3\\+HD2-.*_", ms), paste0(ms, "_plus-minus"))
  heatmaps(paste0("Foxp3-HD2-.*_", ms),     paste0(ms, "_minus-minus"))
  # heatmaps(paste0("Foxp3\\+HD2-.*_", ms),     paste0(ms, "_plus-minus"))
}

colours = brewer.pal(6, "Set1")
for (tp in c(1, 0.85)) {
  radarcharts("_ILN-L_m1", "M1_ILN-L", save=T)
  radarcharts("_ILN-R_m1", "M1_ILN-R", save=T)
  radarcharts("_CoMLN_m1", "M1_CoMLN", save=T)

  radarcharts("_ILN-L_m2", "M2_ILN-L", save=T)
  radarcharts("_ILN-R_m2", "M2_ILN-R", save=T)
  radarcharts("_CoMLN_m2", "M2_CoMLN", save=T)

  radarcharts("_ILN-L_m3", "M3_ILN-L", save=T)
  radarcharts("_ILN-R_m3", "M3_ILN-R", save=T)
  radarcharts("_CoMLN_m3", "M3_CoMLN", save=T)

  radarcharts("Foxp3\\+HD2\\+_.*_m1", "M1_plus-plus", save=T)
  # radarcharts("Foxp3\\+HD2-_.*_m1",   "M1_plus-minus", save=T)
  # radarcharts("Foxp3-HD2\\+_.*_m1",   "M1_minus-plus", save=T)
  radarcharts("Foxp3-HD2-_.*_m1",     "M1_minus-minus", save=T)

  radarcharts("Foxp3\\+HD2\\+_.*_m2", "M2_plus-plus", save=T)
  # radarcharts("Foxp3\\+HD2-_.*_m2",   "M2_plus-minus", save=T)
  # radarcharts("Foxp3-HD2\\+_.*_m2",   "M2_minus-plus", save=T)
  radarcharts("Foxp3-HD2-_.*_m2",     "M2_minus-minus", save=T)

  radarcharts("Foxp3\\+HD2\\+_.*_m3", "M3_plus-plus", save=T)
  # radarcharts("Foxp3\\+HD2-_.*_m3",   "M3_plus-minus", save=T)
  # radarcharts("Foxp3-HD2\\+_.*_m3",   "M3_minus-plus", save=T)
  radarcharts("Foxp3-HD2-_.*_m3",     "M3_minus-minus", save=T)
}

track_clones(NULL, NULL, immdata=idata)
