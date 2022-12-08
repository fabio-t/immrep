
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names = 1)
print(mid_labels)

idata <- immload(which = "not full")

mids_counts <- read.csv("mids_counts.csv", sep = "\t", row.names = 1)
mid_names <- mid_labels[colnames(mids_counts), "label", drop = T]

overview(immdata = idata)
idata2 <- idata
idata2$data <- idata$data[c("MID20", "MID21", "MID22", "MID28", "MID29", "MID30", "MID36", "MID37", "MID38", "MID39")]
overview(immdata = idata2, dirname = "overview2")
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
#  20,21,22,28,29,30,36,37,38 and 39
heatmaps("(SI2|SI4|SI3|SI9|SI17|SI20|SI21)", "comparison")

track_clones(NULL, NULL, immdata = idata)
track_clones("(MM9|MM9_SPF)$", "MM9_vs_MM9SPF", immdata = idata)
track_clones("(E.coli|E.coli_SPF)$", "Ecoli_vs_EcoliSPF", immdata = idata)
track_clones("MM9$", "MM9", immdata = idata)
track_clones("MM9_SPF$", "MM9_SPF", immdata = idata)
track_clones("E.coli$", "E.coli", immdata = idata)
track_clones("E.coli_SPF$", "E.coli_SPF", immdata = idata)
track_clones("(SI2|SI4|SI3|SI9|SI17|SI20|SI21)", "comparison", immdata = idata)
