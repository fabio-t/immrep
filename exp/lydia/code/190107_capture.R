
source("../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names <- c(paste0("m78_", c("IL17+", "IFNy-_IL17-", "IFNy+", "Foxp3+")),
               paste0("m80_", c("IL17+", "IFNy-_IL17-", "IFNy+", "Foxp3+")),
               paste0("m81_", c("IL17+", "IFNy-_IL17-", "IFNy+", "Foxp3+")),
               paste0("m88_", c("IL17+", "IFNy-_IL17-", "IFNy+", "Foxp3+"))
           )

# rarefaction
make_rarefaction_plots(mids_counts, rarefy = T)

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, prefix = "full", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# full heatmap without IFNy- IL17-
columns_i = grep("IFNy-_IL17-", mid_names, invert = T)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn", columns_n, prefix = "noIFNy-_IL17-", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(data, "intersect", columns_n, prefix = "noIFNy-_IL17-", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(data, "intersect", columns_n, prefix = "noIFNy-_IL17-_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# full heatmap without IFNy- IL17-
columns_i = grep("(m81_|IFNy-_IL17-)", mid_names, invert = T)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn", columns_n, prefix = "no_m81_IFNy-_IL17-", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(data, "intersect", columns_n, prefix = "no_m81_IFNy-_IL17-", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(data, "intersect", columns_n, prefix = "no_m81_IFNy-_IL17-_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

quit(save="no")

# TODO: need Valpha8 together with full repertoire

# diversity
make_diversity(mids_counts, mid_names)

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours = brewer.pal(9, "Set1")

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
    for (m in c(78, 80, 81, 88)) {
        run(paste0("m", m, "_(IFNy|IL17|Foxp3)\\+$"), paste0("m", m, "_noIFNy-_IL17-"))
        run(paste0("m", m, "_"),  paste0("m", m))
    }

    run("IFNy\\+", "IFNy+")
    run("IL17\\+", "IL17+")
    run("Foxp3", "Foxp3+")
}
