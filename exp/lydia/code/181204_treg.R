
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

library(matrixStats)

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

celltypes = c("hCD2+", "hCD2-", "GFP+")

mid_names = c("input_r1",
              paste0(paste0(celltypes, "_r1"), "_m1"),
              paste0(paste0(celltypes, "_r1"), "_m2"),
              paste0(paste0(celltypes, "_r1"), "_m3"),

              "input_r2_a", "input_r2_b",
              paste0(paste0(celltypes, "_r2"), "_m1"),
              paste0(paste0(celltypes, "_r2"), "_m2"),
              paste0(paste0(celltypes, "_r2"), "_m3"),

              paste0(paste0(celltypes, "_r3"), "_m1"),
              paste0(paste0(celltypes, "_r3"), "_m2"),
              paste0(paste0(celltypes, "_r3"), "_m3"),

              "hCD2-_c1",
              "hCD2+_c2", "hCD2-_c2",
              "hCD2+_c3", "hCD2-_c3"
          )

# make_rarefaction_plots(mids_counts, grep("input", mid_names), rarefy=T)
make_diversity(mids_counts, mid_names)

heatmaps(NULL, NULL, cex = 0.5, tsne = T)
heatmaps("(GFP|input)", "GFP+") # GFP+ and inputs
heatmaps("GFP", "GFP-only") # GFP+
heatmaps("hCD2", "hCD2") # Ctrls, hCD2- and hCD2+
heatmaps("hCD2.*_r", "hCD2_nocontrol") # hCD2- and hCD2+, without controls
heatmaps("hCD2-", "hCD2-") # Ctrls and hCD2-
heatmaps("hCD2-_r", "hCD2-_nocontrol") # hCD2-, without controls
heatmaps("hCD2\\+", "hCD2+") # Ctrls and hCD2+
heatmaps("hCD2\\+_r", "hCD2+_nocontrol") # hCD2+, without controls

## radarcharts

col = brewer.pal(6, "Set1")

for (tp in c(1, 0.85)) {
    radarcharts("(GFP\\+_r|input_r)1", "GFP+_r1_input", save = T, colours = col)
    radarcharts("(GFP\\+_r|input_r)2", "GFP+_r2_input", save = T, colours = col)
    radarcharts("(GFP\\+_r|input_r)3", "GFP+_r3_input", save = T, colours = col)

    radarcharts("GFP\\+_r1", "GFP+_r1", save = T, colours = col)
    radarcharts("GFP\\+_r2", "GFP+_r2", save = T, colours = col)
    radarcharts("GFP\\+_r3", "GFP+_r3", save = T, colours = col)

    radarcharts("GFP.*_m1", "GFP+_m1", save = T, colours = col)
    radarcharts("GFP.*_m2", "GFP+_m2", save = T, colours = col)
    radarcharts("GFP.*_m3", "GFP+_m3", save = T, colours = col)

    radarcharts("_c1", "ctrl_1", save = T, colours = col)
    radarcharts("_c2", "ctrl_2", save = T, colours = col)
    radarcharts("_c3", "ctrl_3", save = T, colours = col)

    radarcharts("hCD2.*_r1", "hCD2_r1", save = T, colours = col)
    radarcharts("hCD2.*_r2", "hCD2_r2", save = T, colours = col)
    radarcharts("hCD2.*_r3", "hCD2_r3", save = T, colours = col)

    radarcharts("hCD2.*(_r1_m1|_c2)", "rand_hCD2", save = T, colours = col)
}
