
source("../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

colnames(mids_counts) <- c(
    "mid1", "mid2", "mid3", "mid4", "mid5",
    "mid6", "mid7", "mid8", "mid9",
    "mid10", "mid11", "mid13", "mid14", "mid15",
    "mid16", "mid17", "mid19", "mid20", "mid21",
    "mid22", "mid23", "mid24", "mid25",
    "mid26", "mid27", "mid28", "mid29", "mid30",
    "mid31", "mid32", "mid33", "mid34", "mid35",
    "mid36", "mid37", "mid38", "mid39", "mid40"
)

mid_names <- c("1_culture_1810b", "1_culture_1611a", "1_culture_1611b",
               "1_input1", "1_culture_1810a", "1_input2", "1_input3", "1_input4", "1_input5",
               "1_SPL_1", "1_SPL_2", "1_SPL_3", "1_SPL_4", "1_SPL_5",
               "1_LI_1", "1_LI_2", "1_LI_3", "1_LI_4", "1_LI_5",
               "2_culture_1810a", "2_culture_1810b", "2_culture_1611a", "2_culture_1611b",
               "2_input1", "2_input2", "2_input3", "2_input4", "2_input5",
               "2_SPL_1", "2_SPL_2", "2_SPL_3", "2_SPL_4", "2_SPL_5",
               "2_LI_1", "2_LI_2", "2_LI_3", "2_LI_4", "2_LI_5"
           )

### removing MIDs 4, 6, 7, 8 because they are failed samples
mids_counts <- mids_counts[-c(4,6,7,8)]
mid_names <- mid_names[-c(4,6,7,8)]

# rarefaction
# make_rarefaction_plots(mids_counts, c(1, 2, 3, 5))

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, prefix = "full", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# only cultures and input 5
mid_ids <- grep("(culture|input5)", mid_names)
make_full_heatmap(mids_counts[mid_ids], "horn", mid_names[mid_ids], prefix = "cultures", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts[mid_ids], "intersect", mid_names[mid_ids], prefix = "cultures", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts[mid_ids], "intersect", mid_names[mid_ids], prefix = "cultures_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# no cultures, no inputs
mid_ids <- grep("(culture|input)", mid_names, invert = T)
make_full_heatmap(mids_counts[mid_ids], "horn", mid_names[mid_ids], prefix = "organs", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts[mid_ids], "intersect", mid_names[mid_ids], prefix = "organs", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts[mid_ids], "intersect", mid_names[mid_ids], prefix = "organs_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# TODO: need Valpha8 together with full repertoire

# diversity
make_diversity(mids_counts, mid_names)

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours = brewer.pal(6, "Set1")

run <- function(pattern, dirname) {
    print(dirname)

    column_indices = grep(pattern, mid_names)

    column_names   = colnames(mids_counts)[column_indices]
    column_names_2 = mid_names[column_indices]

    print(column_names_2)

    data = mids_counts[column_names]

    save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)
}

for (tp in c(1, 0.85)) {
    run("1810a",  "1810a")
    run("1810b",  "1810b")
    run("1611a",  "1611a")
    run("1611b",  "1611b")
    run("input5", "input5")
}
