
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
# make_diversity(mids_counts, mid_names)

## heatmaps
#
# heatmaps(NULL, NULL, cex = 0.5)
# heatmaps("TCRa$", "TCRa")
# heatmaps("TCRb$", "TCRb")
# heatmaps("m1_TCRa$", "M1_TCRa")
# heatmaps("m1_TCRb$", "M1_TCRb")
# heatmaps("m2_TCRa$", "M2_TCRa")
# heatmaps("m2_TCRb$", "M2_TCRb")

## radarcharts

# for (tp in c(1, 0.85)) {
#     radarcharts("m1_TCRa$", "m1_TCRa", save = T)
#     radarcharts("m1_TCRb$", "m1_TCRb", save = T)
#
#     radarcharts("m2_TCRa$", "m2_TCRa", save = T)
#     radarcharts("m2_TCRb$", "m2_TCRb", save = T)
# }

# overall frequencies
data_percs <- sweep(mids_counts, 2, colSums(mids_counts), "/")
colnames(data_percs) = mid_names
data_percs_top <- get_top(data_percs, topN=20, min_clones=1)
data_percs_top[data_percs_top == 0] = NA
sc_cols <- grep("^sc_", mid_names, value=T)
data_percs_top[sc_cols] = data_percs[rownames(data_percs_top), sc_cols] # restore single cell values
data_percs_long <- melt(data_percs_top, id=sc_cols)
data_percs_long <- separate(data=data_percs_long, col = variable, into = c("id", "foxp3", "type", "mouse", "chain"), sep="_")

library(ggpubr)

# GFP m1
df <- subset(data_percs_long, grepl("GFP", type) & mouse == "m1", select=c(sc_GFP_m1_TCRa, sc_GFP_m1_TCRb, id:value))
df <- melt(df, id=c("id", "foxp3", "type", "mouse", "chain", "value"), value.name="y", variable.name="singlecell")
df <- subset(df, (chain == "TCRa" & singlecell == "sc_GFP_m1_TCRa") | (chain == "TCRb" & singlecell == "sc_GFP_m1_TCRb"))
p <- ggscatter(df, x="value", y="y", add = "reg.line", conf.int = T,
               color = "singlecell", shape = "singlecell", palette = c("#CB181D", "#2171B5")) +
          stat_cor(aes(color = singlecell))
ggsave(p, file="gfp_m1.svg")
ggsave(p + xscale("log10", .format = TRUE) + yscale("log10", .format = TRUE), file="gfp_m1.log.svg")

# GFP m2
df <- subset(data_percs_long, grepl("GFP", type) & mouse == "m2", select=c(sc_GFP_m2_TCRa, sc_GFP_m2_TCRb, id:value))
df <- melt(df, id=c("id", "foxp3", "type", "mouse", "chain", "value"), value.name="y", variable.name="singlecell")
df <- subset(df, (chain == "TCRa" & singlecell == "sc_GFP_m2_TCRa") | (chain == "TCRb" & singlecell == "sc_GFP_m2_TCRb"))
p <- ggscatter(df, x="value", y="y", add = "reg.line", conf.int = T,
               color = "singlecell", shape = "singlecell", palette = c("#CB181D", "#2171B5")) +
          stat_cor(aes(color = singlecell))
ggsave(p, file="gfp_m2.svg")
ggsave(p + xscale("log10", .format = TRUE) + yscale("log10", .format = TRUE), file="gfp_m2.log.svg")

# CD45 m1
df <- subset(data_percs_long, grepl("(Teff)", type) & mouse == "m1", select=c(sc_CD45_m1_TCRa, sc_CD45_m1_TCRb, id:value))
df <- melt(df, id=c("id", "foxp3", "type", "mouse", "chain", "value"), value.name="y", variable.name="singlecell")
df <- subset(df, (chain == "TCRa" & singlecell == "sc_CD45_m1_TCRa") | (chain == "TCRb" & singlecell == "sc_CD45_m1_TCRb"))
p <- ggscatter(df, x="value", y="y", add = "reg.line", conf.int = T,
               color = "singlecell", shape = "singlecell", palette = c("#CB181D", "#2171B5")) +
          stat_cor(aes(color = singlecell))
ggsave(p, file="cd45_m1.svg")
ggsave(p + xscale("log10", .format = TRUE) + yscale("log10", .format = TRUE), file="cd45_m1.log.svg")

# CD45 m2
df <- subset(data_percs_long, grepl("(Teff)", type) & mouse == "m2", select=c(sc_CD45_m2_TCRa, sc_CD45_m2_TCRb, id:value))
df <- melt(df, id=c("id", "foxp3", "type", "mouse", "chain", "value"), value.name="y", variable.name="singlecell")
df <- subset(df, (chain == "TCRa" & singlecell == "sc_CD45_m2_TCRa") | (chain == "TCRb" & singlecell == "sc_CD45_m2_TCRb"))
p <- ggscatter(df, x="value", y="y", add = "reg.line", conf.int = T,
               color = "singlecell", shape = "singlecell", palette = c("#CB181D", "#2171B5")) +
          stat_cor(aes(color = singlecell))
ggsave(p, file="cd45_m2.svg")
ggsave(p + xscale("log10", .format = TRUE) + yscale("log10", .format = TRUE), file="cd45_m2.log.svg")
