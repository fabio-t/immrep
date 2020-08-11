
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts <- read.csv("mids_counts.csv", sep = "\t", header = F)
rownames(mids_counts) <- mids_counts$V1
mids_counts$V1 <- NULL
colnames(mids_counts) <- c(
    "mid22", "mid23", "mid24", "mid25",
    "mid26", "mid27", "mid28", "mid29", "mid30",
    "mid31", "mid32", "mid33", "mid34", "mid35",
    "mid36", "mid37", "mid38", "mid39", "mid40"
)
f <- colwise(get_freqs)
mids_freqs <- f(mids_counts)

mid_names <- c("culture_1810a", "culture_1810b", "culture_1611a", "culture_1611b",
               "input1", "input2", "input3", "input4", "input5",
               "SPL_1", "SPL_2", "SPL_3", "SPL_4", "SPL_5",
               "LI_1", "LI_2", "LI_3", "LI_4", "LI_5"
           )

# rarefaction
# make_rarefaction_plots(mids_counts, c(1, 2, 3, 4))

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, prefix = "full", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# TODO: need Valpha8 together with full repertoire

# diversity
divers_df <- as.data.frame(mid_names)
rownames(divers_df) <- colnames(mids_counts)
colnames(divers_df) <- "Sample"
divers_df$Shannon <- diversity(t(mids_counts), index = "shannon")
divers_df$InvSimpson <- diversity(t(mids_counts), index = "invsimpson")
write.csv(divers_df, file = "diversity.csv")
