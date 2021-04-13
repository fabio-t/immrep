
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c(

  "1_FoxP3-_hCD2-_Input_1",
  "2_FoxP3-_hCD2-_Input_2",
  "3_FoxP3+_hCD2+_Input_1",
  "4_FoxP3+_hCD2+_Input_2",

  "1_FoxP3+_GFP1_m1_TCRa",
  "2_FoxP3+_GFP2_m1_TCRa",
  "3_FoxP3+_GFP3_m1_TCRa",
  "4_FoxP3-_Teff1_m1_TCRa",
  "5_FoxP3+_Teff2_m1_TCRa",
  "6_FoxP3+_Teff3_m1_TCRa",
  "7_FoxP3+_iTreg_m1_TCRa",

  "1_FoxP3+_GFP1_m2_TCRa",
  "2_FoxP3+_GFP2_m2_TCRa",
  "3_FoxP3+_GFP3_m2_TCRa",
  "4_FoxP3-_Teff1_m2_TCRa",
  "5_FoxP3+_Teff2_m2_TCRa",
  "6_FoxP3+_Teff3_m2_TCRa",
  "7_FoxP3+_iTreg_m2_TCRa",

  "1_FoxP3+_GFP1_m1_TCRb",
  "2_FoxP3+_GFP2_m1_TCRb",
  "3_FoxP3+_GFP3_m1_TCRb",
  "4_FoxP3-_Teff1_m1_TCRb",
  "5_FoxP3+_Teff2_m1_TCRb",
  "6_FoxP3+_Teff3_m1_TCRb",
  "7_FoxP3+_iTreg_m1_TCRb",

  "1_FoxP3+_GFP1_m2_TCRb",
  "2_FoxP3+_GFP2_m2_TCRb",
  "3_FoxP3+_GFP3_m2_TCRb",
  "4_FoxP3-_Teff1_m2_TCRb",
  "5_FoxP3+_Teff2_m2_TCRb",
  "6_FoxP3+_Teff3_m2_TCRb",
  "7_FoxP3+_iTreg_m2_TCRb"
)

# make_rarefaction_plots(mids_counts, rarefy=T)
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("TCRa$", "TCRa")
heatmaps("(TCRa|Input_.)$", "TCRa-inputs")
heatmaps("TCRb$", "TCRb")
heatmaps("m1_TCRa$", "M1_TCRa")
heatmaps("(m1_TCRa|Input_.)$", "M1_TCRa-inputs")
heatmaps("m1_TCRb$", "M1_TCRb")
heatmaps("m2_TCRa$", "M2_TCRa")
heatmaps("(m2_TCRa|Input_.)$", "M2_TCRa-inputs")
heatmaps("m2_TCRb$", "M2_TCRb")

## radarcharts

for (tp in c(1, 0.85)) {
    radarcharts("m1_TCRa$", "m1_TCRa", save = T)
    radarcharts("(m1_TCRa|Input_.)$", "m1_TCRa-inputs", save = T)
    radarcharts("m1_TCRb$", "m1_TCRb", save = T)

    radarcharts("m2_TCRa$", "m2_TCRa", save = T)
    radarcharts("(m2_TCRa|Input_.)$", "m2_TCRa-inputs", save = T)
    radarcharts("m2_TCRb$", "m2_TCRb", save = T)
}
