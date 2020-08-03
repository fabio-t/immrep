
source("../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names <- c("culture_1810b", "culture_1611a", "culture_1611b",
               "input1", "culture_1810a", "input2", "input3", "input4", "input5",
               "SPL_1", "SPL_2", "SPL_3", "SPL_4", "SPL_5",
               "LI_1", "LI_2", "LI_3", "LI_4", "LI_5"
           )

# rarefaction
# make_rarefaction_plots(mids_counts, rarefy=T)

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, prefix = "full", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# TODO: need Valpha8 together with full repertoire

# diversity
make_diversity(mids_counts, mid_names)
