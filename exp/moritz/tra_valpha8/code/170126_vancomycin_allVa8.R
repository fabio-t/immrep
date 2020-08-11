
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = c("mid1",  "mid2",  "mid3",

                          "mid4",  "mid5",  "mid6",
                          "mid7",  "mid8",  "mid9",
                          "mid10", "mid11", "mid13",

                          "mid14", "mid15", "mid16",
                          "mid17", "mid19", "mid20",
                          "mid21", "mid22", "mid23",

                          "mid24", "mid25", "mid26",
                          "mid27", "mid28", "mid29",
                          "mid32",

                          "mid33", "mid34", "mid35",
                          "mid36", "mid37", "mid38",
                          "mid39", "mid40", "mid41",

                          "mid42", "mid43", "mid44",

                          "mid45", "mid46", "mid47",
                          "mid48", "mid49", "mid50",

                          "mid51", "mid52", "mid53",
                          "mid54", "mid55", "mid56"
                      )

f = colwise(get_freqs)
mids_freqs = f(mids_counts)

mid_names_va8 <- c("input1", "input2", "input3",

                  "Ctrl1_MLN_FoxP3-", "Ctrl2_MLN_FoxP3-", "Ctrl3_MLN_FoxP3-",
                  "Ctrl1_SPL_FoxP3-", "Ctrl2_SPL_FoxP3-", "Ctrl3_SPL_FoxP3-",
                  "Ctrl1_LI_FoxP3-",  "Ctrl2_LI_FoxP3-",  "Ctrl3_LI_FoxP3-",

                  "Ctrl1_MLN_FoxP3+", "Ctrl2_MLN_FoxP3+", "Ctrl3_MLN_FoxP3+",
                  "Ctrl1_SPL_FoxP3+", "Ctrl2_SPL_FoxP3+", "Ctrl3_SPL_FoxP3+",
                  "Ctrl1_LI_FoxP3+",  "Ctrl2_LI_FoxP3+",  "Ctrl3_LI_FoxP3+",

                  "Van1_MLN_FoxP3-",  "Van2_MLN_FoxP3-",  "Van3_MLN_FoxP3-",
                  "Van1_SPL_FoxP3-",  "Van2_SPL_FoxP3-",  "Van3_SPL_FoxP3-",
                  "Van3_LI_FoxP3-",

                  "Van1_MLN_FoxP3+",  "Van2_MLN_FoxP3+",  "Van3_MLN_FoxP3+",
                  "Van1_SPL_FoxP3+",  "Van2_SPL_FoxP3+",  "Van3_SPL_FoxP3+",
                  "Van1_LI_FoxP3+",   "Van2_LI_FoxP3+",   "Van3_LI_FoxP3+"
                )

mid_names_full <- c("input1", "input2", "input3",

                    "Ctrl1_LI_FoxP3-",  "Ctrl2_LI_FoxP3-",  "Ctrl3_LI_FoxP3-",
                    "Ctrl1_LI_FoxP3+",  "Ctrl2_LI_FoxP3+",  "Ctrl3_LI_FoxP3+",

                    "Van1_LI_FoxP3-",   "Van2_LI_FoxP3-",   "Van3_LI_FoxP3-",
                    "Van1_LI_FoxP3+",   "Van2_LI_FoxP3+",   "Van3_LI_FoxP3+"
                )

mid_names <- c(paste0("va8_", mid_names_va8), paste0("full_", mid_names_full))

# rarefaction
# make_rarefaction_plots(mids_counts, c(1, 2, 3, 5))

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, prefix = "full", binary = F, types = "square", cexRow = 0.4, cexCol = 0.4, tsne = T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = 0.4, cexCol = 0.4)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.4, cexCol = 0.4)

# only common
mid_ids <- grep("(input|_LI_)", mid_names)
make_full_heatmap(mids_counts[mid_ids], "horn", mid_names[mid_ids], prefix = "inputLI", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts[mid_ids], "intersect", mid_names[mid_ids], prefix = "inputLI", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts[mid_ids], "intersect", mid_names[mid_ids], prefix = "inputLI_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# only common, no inputs
mid_ids <- grep("_LI_", mid_names)
make_full_heatmap(mids_counts[mid_ids], "horn", mid_names[mid_ids], prefix = "li", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts[mid_ids], "intersect", mid_names[mid_ids], prefix = "li", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts[mid_ids], "intersect", mid_names[mid_ids], prefix = "li_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)


quit(save="no")

# diversity
divers_df <- as.data.frame(mid_names)
rownames(divers_df) <- colnames(mids_counts)
colnames(divers_df) <- "Sample"
divers_df$Shannon <- diversity(t(mids_counts), index = "shannon")
divers_df$InvSimpson <- diversity(t(mids_counts), index = "invsimpson")
write.csv(divers_df, file = "diversity.csv")

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
    run("Ctrl.*LI_FoxP3-", "ctrl_li_foxp3-")
    run("Ctrl.*LI_FoxP3+", "ctrl_li_foxp3+")
    run("Van.*LI_FoxP3-",  "van_li_foxp3-")
    run("Van.*LI_FoxP3+",  "van_li_foxp3+")
}
