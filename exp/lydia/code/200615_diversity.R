
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c(
              "1_FoxP3+_10K_m1",
              "2_FoxP3+_10K_m1",
              "3_FoxP3+_10K_m1",
              "4_FoxP3+_10K_m1",
              "5_FoxP3+_10K_m1",
              "6_FoxP3+_10K_m1",
              "7_FoxP3+_30K_m1",
              "8_FoxP3+_30K_m1",
              "9_FoxP3+_30K_m1",
              "10_FoxP3+_80K_m1",
              "11_FoxP3+_80K_m1",
              "12_FoxP3+_80K_m1",
              "13_FoxP3+_80K_m1",

              "17_FoxP3-_10K_m1",
              "18_FoxP3-_10K_m1",
              "19_FoxP3-_10K_m1",
              "20_FoxP3-_10K_m1",
              "21_FoxP3-_10K_m1",
              "22_FoxP3-_10K_m1",
              "23_FoxP3-_80K_m1",
              "24_FoxP3-_80K_m1",
              "25_FoxP3-_80K_m1",
              "26_FoxP3-_80K_m1",

              "1_FoxP3+_10K_m2",
              "2_FoxP3+_10K_m2",
              "3_FoxP3+_10K_m2",
              "4_FoxP3+_10K_m2",
              "5_FoxP3+_10K_m2",
              "6_FoxP3+_10K_m2",
              "7_FoxP3+_30K_m2",
              "8_FoxP3+_30K_m2",
              "9_FoxP3+_30K_m2",
              "10_FoxP3+_30K_m2",
              "11_FoxP3+_80K_m2",
              "12_FoxP3+_80K_m2",

              "15_FoxP3+_LI_m2",

              "17_FoxP3-_10K_m2",
              "18_FoxP3-_10K_m2",
              "19_FoxP3-_10K_m2",
              "20_FoxP3-_10K_m2",
              "21_FoxP3-_10K_m2",
              "22_FoxP3-_10K_m2",
              "23_FoxP3-_80K_m2",
              "24_FoxP3-_80K_m2",
              "25_FoxP3-_80K_m2"
          )

# make_rarefaction_plots(mids_counts, rarefy=T)

## diversity

make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("_m1$", "M1")
heatmaps("_m2$", "M2")
heatmaps("FoxP3\\+_.*_m1$", "FoxP3+_M1")
heatmaps("FoxP3\\+_.*_m2$", "FoxP3+_M2")
heatmaps("FoxP3-_.*_m1$", "FoxP3-_M1")
heatmaps("FoxP3-_.*_m2$", "FoxP3-_M2")
heatmaps("FoxP3\\+_80K_", "Foxp3+_80K")
heatmaps("FoxP3-_80K_",   "Foxp3-_80K")
heatmaps("_80K_",   "80K_m1-2")

## radarcharts

for (tp in c(1, 0.85)) {
    radarcharts("FoxP3\\+_.*_m1$", "foxp3+_m1", save = T)
    radarcharts("FoxP3-_.*_m1$", "foxp3-_m1", save = T)
    radarcharts("FoxP3\\+_.*_m2$", "foxp3+_m2", save = T)
    radarcharts("FoxP3-_.*_m2$", "foxp3-_m2", save = T)

    radarcharts("FoxP3\\+_10K_m1$", "foxp3+_m1_10K", save = T)
    radarcharts("FoxP3\\+_30K_m1$", "foxp3+_m1_30K", save = T)
    radarcharts("FoxP3\\+_80K_m1$", "foxp3+_m1_80K", save = T)

    radarcharts("FoxP3-_10K_m1$", "foxp3-_m1_10K", save = T)
    radarcharts("FoxP3-_80K_m1$", "foxp3-_m1_80K", save = T)

    radarcharts("FoxP3\\+_10K_m2$", "foxp3+_m2_10K", save = T)
    radarcharts("FoxP3\\+_30K_m2$", "foxp3+_m2_30K", save = T)
    radarcharts("FoxP3\\+_80K_m2$", "foxp3+_m2_80K", save = T)

    radarcharts("FoxP3-_10K_m2$", "foxp3-_m2_10K", save = T)
    radarcharts("FoxP3-_80K_m2$", "foxp3-_m2_80K", save = T)

    radarcharts("_10K_m1", "m1_10K", save = T)
    radarcharts("_80K_m1", "m1_80K", save = T)

    radarcharts("_10K_m2", "m2_10K", save = T)
    radarcharts("_80K_m2", "m2_80K", save = T)

    radarcharts("FoxP3\\+_80K_", "foxp3+_80K", save = T)
    radarcharts("FoxP3-_80K_",   "foxp3-_80K", save = T)
    radarcharts("_80K_",         "80K_m1-2",   save = T)
}
