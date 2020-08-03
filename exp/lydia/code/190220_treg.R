
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

library(matrixStats)

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

celltypes = c("hCD2+", "hCD2-", "GFP+")

mid_names = c(paste0(celltypes, "_m25_LI_tcra"),
              paste0(celltypes, "_m27_LI_tcra"),
              paste0(celltypes, "_m21_LI_tcra"),
              paste0(celltypes, "_m22_LI_tcra"),
              paste0(celltypes, "_m23_LI_tcra"),
              paste0(celltypes, "_m24_LI_tcra"),

              paste0(celltypes, "_m25_LI_tcrb"),
              paste0(celltypes, "_m27_LI_tcrb"),
              paste0(celltypes, "_m21_LI_tcrb"),
              paste0(celltypes, "_m22_LI_tcrb"),
              paste0(celltypes, "_m23_LI_tcrb"),
              paste0(celltypes, "_m24_LI_tcrb")
          )

# make_rarefaction_plots(mids_counts, rarefy=T)

## All
make_full_heatmap(mids_counts, "horn",      mid_names, prefix="full", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## All TRA
columns_i = grep("_tcra", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="tcra", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="tcra", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="tcra_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## All TRB
columns_i = grep("_tcrb", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="tcrb", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="tcrb", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="tcrb_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## GFP+ TRA
columns_i = grep("GFP.*_tcra", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="GFP+_tcra", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="GFP+_tcra", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="GFP+_tcra_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## GFP+ TRB
columns_i = grep("GFP.*_tcrb", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="GFP+_tcrb", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="GFP+_tcrb", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="GFP+_tcrb_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## hCD2- and hCD2+ TRA
columns_i = grep("hCD2.*_tcra", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="hCD2_tcra", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2_tcra", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2_tcra_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## hCD2- and hCD2+ TRB
columns_i = grep("hCD2.*_tcrb", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="hCD2_tcrb", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2_tcrb", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2_tcrb_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# diversity
make_diversity(mids_counts, mid_names)

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours = brewer.pal(6, "Set1")

run <- function(regexp, dirname, save = F) {
    print(dirname)

    column_indices = grep(regexp, mid_names)

    column_names   = colnames(mids_counts)[column_indices]
    column_names_2 = mid_names[column_indices]

    print(column_names_2)

    data = mids_counts[column_names]

    l <- save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)

    if (save && (length(l$clones) > 0)) {
        d1 <- sweep(data, 2, colSums(data), "/")[l$clones, ]
        d2 <- data[l$clones, ]
        d3 <- as.data.frame(colRanks(as.matrix(d2), preserveShape=T, ties.method="average"))
        colnames(d1) <- column_names_2
        colnames(d2) <- column_names_2
        colnames(d3) <- column_names_2; rownames(d3) <- l$clones
        write.csv(d1, file = paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname, "/data_perc.csv"))
        write.csv(d2, file = paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname, "/data_counts.csv"))
        write.csv(d3, file = paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname, "/data_ranks.csv"))
    }
}

for (tp in c(1, 0.85)) {
    run("_m25_LI_tcra", "m25_tcra", save = T)
    run("_m27_LI_tcra", "m27_tcra", save = T)
    run("_m21_LI_tcra", "m21_tcra", save = T)
    run("_m22_LI_tcra", "m22_tcra", save = T)
    run("_m23_LI_tcra", "m23_tcra", save = T)
    run("_m24_LI_tcra", "m24_tcra", save = T)

    run("_m25_LI_tcrb", "m25_tcrb", save = T)
    run("_m27_LI_tcrb", "m27_tcrb", save = T)
    run("_m21_LI_tcrb", "m21_tcrb", save = T)
    run("_m22_LI_tcrb", "m22_tcrb", save = T)
    run("_m23_LI_tcrb", "m23_tcrb", save = T)
    run("_m24_LI_tcrb", "m24_tcrb", save = T)
}
