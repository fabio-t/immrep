
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

# MID1    MID2    MID3    MID4    MID5    MID6    MID7    MID8    MID9    MID10
# MID11   MID13   MID14   MID15   MID16   MID17   MID19   MID20   MID21   MID22
# MID23   MID24   MID25   MID26   MID27   MID28   MID29   MID30
# MID111  MID121  MID211  MID221  MID311  MID321  MID411  MID421

mid_names = c(
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
  "7_FoxP3+_iTreg_m2_TCRb",

  # now the four single cells
  "sc_GFP_m1_TCRa",
  "sc_GFP_m1_TCRb",
  "sc_CD45_m1_TCRa",
  "sc_CD45_m1_TCRb",
  "sc_GFP_m2_TCRa",
  "sc_GFP_m2_TCRb",
  "sc_CD45_m2_TCRa",
  "sc_CD45_m2_TCRb"
)

# make_rarefaction_plots(mids_counts, rarefy=T)
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("TCRa$", "TCRa")
heatmaps("TCRb$", "TCRb")
heatmaps("m1_TCRa$", "M1_TCRa")
heatmaps("m1_TCRb$", "M1_TCRb")
heatmaps("m2_TCRa$", "M2_TCRa")
heatmaps("m2_TCRb$", "M2_TCRb")

## radarcharts

for (tp in c(1, 0.85)) {
    radarcharts("m1_TCRa$", "m1_TCRa", save = T)
    radarcharts("m1_TCRb$", "m1_TCRb", save = T)

    radarcharts("m2_TCRa$", "m2_TCRa", save = T)
    radarcharts("m2_TCRb$", "m2_TCRb", save = T)
}
