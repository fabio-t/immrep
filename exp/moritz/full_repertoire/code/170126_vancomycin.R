
library(rprojroot)
root <- is_vcs_root$make_fix_file()
source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c("input1", "input2", "input3",

              "Ctrl1_LI_FoxP3-",  "Ctrl2_LI_FoxP3-",  "Ctrl3_LI_FoxP3-",
              "Ctrl1_LI_FoxP3+",  "Ctrl2_LI_FoxP3+",  "Ctrl3_LI_FoxP3+",

              "Van1_LI_FoxP3-",  "Van2_LI_FoxP3-",  "Van3_LI_FoxP3-",
              "Van1_LI_FoxP3+",  "Van2_LI_FoxP3+",  "Van3_LI_FoxP3+")


# rarefaction
# make_rarefaction_plots(mids_counts, 1:3)

# diversity
make_diversity(mids_counts, mid_names)

# heatmaps
heatmaps(NULL, NULL, cex = 0.5, tsne = T)
heatmaps("Ctrl", "ctrl")
heatmaps("Van", "van")
heatmaps("LI", "li")
heatmaps("LI_FoxP3-", "li_foxp3-")

# radarcharts
for (tp in c(1, 0.85)) {
  radarcharts("input", "inputs", save = T)
  radarcharts("LI_FoxP3-", "li_foxp3-", save = T)
  radarcharts("LI_FoxP3\\+", "li_foxp3+", save = T)
  radarcharts("Ctrl._LI", "ctrl_li", save = T)
}
