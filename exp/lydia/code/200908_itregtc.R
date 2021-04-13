
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

library(matrixStats)

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c("FoxP3+_Pool1_Input_1", "FoxP3+_Pool1_Input_2",
              "FoxP3-_Pool1_Input_1", "FoxP3-_Pool1_Input_2",
              "FoxP3+_Pool2_Input_1", "FoxP3-_Pool2_Input_1",

              "FoxP3+_Pool1_m91_1", "FoxP3-_Pool1_m91_1",
              "FoxP3+_Pool1_m92_1", "FoxP3-_Pool1_m92_1",
              "FoxP3+_Pool2_m93_1", "FoxP3-_Pool2_m93_1",

              "FoxP3+_Pool1_m94_1", "FoxP3+_Pool1_m94_2", "FoxP3+_Pool1_m94_3",
              "FoxP3-_Pool1_m94_1", "FoxP3-_Pool1_m94_2", "FoxP3-_Pool1_m94_3",

              "FoxP3+_Pool2_m95_1", "FoxP3-_Pool2_m95_1"
          )

# rarefaction
make_rarefaction_plots(mids_counts, rarefy=T)

# diversity
make_diversity(mids_counts, mid_names)

# heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("Pool1",   "pool1")
heatmaps("Pool1_m", "pool1-noinput")
heatmaps("Pool2",   "pool2")
heatmaps("Pool2_m", "pool2-noinput")
heatmaps("m94_",    "m94")

# radarcharts

for (tp in c(1, 0.85)) {
    radarcharts("Pool1",   "pool1", save = T)
    radarcharts("Pool1_m", "pool1_noinput", save = T)
    radarcharts("Pool2",   "pool2", save = T)
    radarcharts("Pool2_m", "pool2_noinput", save = T)
    radarcharts("FoxP3\\+_Pool1", "foxp3+_pool1", save = T)
    radarcharts("FoxP3-_Pool1", "foxp3-_pool1", save = T)
    radarcharts("FoxP3\\+_Pool2", "foxp3+_pool2", save = T)
    radarcharts("FoxP3-_Pool2", "foxp3-_pool2", save = T)
}
