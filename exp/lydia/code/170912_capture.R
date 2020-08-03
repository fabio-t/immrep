
source("../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names <- c(paste0("m1_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")),
               paste0("m2_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")),
               paste0("m4_", c("IFNy-_IL17-", "IFNy+_IL17+",    "Foxp3+")),
               paste0("m5_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")),
               paste0("m6_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")),
               paste0("m7_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+"))
           )

# rarefaction
# make_rarefaction_plots(mids_counts, rarefy = T)

# full heatmap
make_full_heatmap(mids_counts, "horn",      mid_names, prefix = "full", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# full heatmap without IFNy- IL17-
columns_i = grep("IFNy-_IL17-", mid_names, invert = T)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix = "noIFNy-_IL17-", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(data, "intersect", columns_n, prefix = "noIFNy-_IL17-", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(data, "intersect", columns_n, prefix = "noIFNy-_IL17-_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# requested on 6/02/2020: heatmap with m1, m2, m5, m6, m7 only IFN+, IL-17+ and Foxp3+
# so, essentially, exclude anything with [+-]_IL17 pattern and everything from m4
columns_i = grep("([\\+-]_IL17|m4_)", mid_names, invert = T)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix = "special", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(data, "intersect", columns_n, prefix = "special", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(data, "intersect", columns_n, prefix = "special_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# requested on 14/05/2020
columns_i = grep("^(m1|m2|m5|)_(IFNy\\+|IL17\\+|Foxp3\\+)", mid_names)
columns_n = mid_names[columns_i]
print(columns_n)
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix = "m1-2-5_only-plus", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(data, "intersect", columns_n, prefix = "m1-2-5_only-plus", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(data, "intersect", columns_n, prefix = "m1-2-5_only-plus_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

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
    for (m in c(1,2,4,5,6,7)) {
        run(paste0("m", m, "_(IFNy|IL17|Foxp3)\\+$"), paste0("m", m, "_noIFNy-_IL17-"))
        run(paste0("m", m, "_"),  paste0("m", m))
    }

    run("IFNy\\+", "IFNy+")
    run("IL17\\+", "IL17+")
    run("Foxp3", "Foxp3+")
}
