
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c(
  "Teff_HZI_1",
  "Teff_HZI_2",
  "Teff_HZI_3",

  "Teff_Prev_1",
  "Teff_Prev_2",
  "Teff_Prev_3",

  "Treg_HZI_1",
  "Treg_HZI_2",
  "Treg_HZI_3",
  "Treg_HZI_4",

  "Treg_Prev_1",
  "Treg_Prev_1"
)

# make_rarefaction_plots(mids_counts, rarefy=T)
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("Teff", "Teff")
heatmaps("Treg", "Treg")
heatmaps("HZI",  "HZI")
heatmaps("Prev", "Prevotella")

## radarcharts

for (tp in c(1, 0.85)) {
    radarcharts("Teff", "Teff", save = T)
    radarcharts("Treg", "Treg", save = T)
    radarcharts("HZI",  "HZI", save = T)
    radarcharts("Prev", "Prev", save = T)
}
