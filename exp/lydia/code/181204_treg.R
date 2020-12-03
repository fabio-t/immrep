
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

library(matrixStats)

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

celltypes = c("hCD2+", "hCD2-", "GFP+")

mid_names = c("input_r1",
              paste0(paste0(celltypes, "_r1"), "_m1"),
              paste0(paste0(celltypes, "_r1"), "_m2"),
              paste0(paste0(celltypes, "_r1"), "_m3"),

              "input_r2_a", "input_r2_b",
              paste0(paste0(celltypes, "_r2"), "_m1"),
              paste0(paste0(celltypes, "_r2"), "_m2"),
              paste0(paste0(celltypes, "_r2"), "_m3"),

              paste0(paste0(celltypes, "_r3"), "_m1"),
              paste0(paste0(celltypes, "_r3"), "_m2"),
              paste0(paste0(celltypes, "_r3"), "_m3"),

              "hCD2-_c1",
              "hCD2+_c2", "hCD2-_c2",
              "hCD2+_c3", "hCD2-_c3"
          )

make_rarefaction_plots(mids_counts, grep("input", mid_names), rarefy=T)

## All
make_full_heatmap(mids_counts, "horn",      mid_names, prefix="full", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## GFP+ and inputs
columns_i = grep("(GFP|input)", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="GFP+", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="GFP+", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="GFP+_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## GFP+
columns_i = grep("GFP", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="GFP-only", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="GFP-only", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="GFP-only_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## Ctrls, hCD2- and hCD2+
columns_i = grep("hCD2", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="hCD2", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## Ctrls and hCD2-
columns_i = grep("hCD2-", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="hCD2-", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2-", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2-_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## hCD2-, without controls
columns_i = grep("hCD2-_r", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="hCD2-_nocontrol", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2-_nocontrol", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2-_nocontrol_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## Ctrls and hCD2+
columns_i = grep("hCD2\\+", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="hCD2+", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2+", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2+_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## hCD2+, without controls
columns_i = grep("hCD2\\+_r", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="hCD2+_nocontrol", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2+_nocontrol", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="hCD2+_nocontrol_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

quit(save="no")

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
    run("(GFP\\+_r|input_r)1", "GFP+_r1_input")
    run("(GFP\\+_r|input_r)2", "GFP+_r2_input")
    run("(GFP\\+_r|input_r)3", "GFP+_r3_input")

    run("GFP\\+_r1", "GFP+_r1", save = T)
    run("GFP\\+_r2", "GFP+_r2", save = T)
    run("GFP\\+_r3", "GFP+_r3", save = T)

    run("GFP.*_m1", "GFP+_m1", save = T)
    run("GFP.*_m2", "GFP+_m2", save = T)
    run("GFP.*_m3", "GFP+_m3", save = T)

    run("_c1", "ctrl_1")
    run("_c2", "ctrl_2")
    run("_c3", "ctrl_3")

    run("hCD2.*_r1", "hCD2_r1")
    run("hCD2.*_r2", "hCD2_r2")
    run("hCD2.*_r3", "hCD2_r3")

    run("hCD2.*(_r1_m1|_c2)", "rand_hCD2")
}
