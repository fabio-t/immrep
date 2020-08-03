
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

library(matrixStats)

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c("LI_r1_m1", "LI_r1_m2", "LI_r1_m3",
              "LI_r2_m4", "LI_r2_m5", "LI_r2_m6",
              "LI_r3_m7", "LI_r3_m8", "LI_r3_m9",

              "MLN_r1_m1", "MLN_r1_m2", "MLN_r1_m3",
              "MLN_r2_m4", "MLN_r2_m5", "MLN_r2_m6",
              "MLN_r3_m7", "MLN_r3_m8", "MLN_r3_m9"
          )

# make_rarefaction_plots(mids_counts, rarefy=T)

## All
make_full_heatmap(mids_counts, "horn",      mid_names, prefix="full", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

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
    run("r1", "r1_LI-MLN", save = T)
    run("r2", "r2_LI-MLN", save = T)
    run("r3", "r3_LI-MLN", save = T)

    run("m1", "m1_LI-MLN", save = T)
    run("m2", "m2_LI-MLN", save = T)
    run("m3", "m3_LI-MLN", save = T)
    run("m4", "m4_LI-MLN", save = T)
    run("m5", "m5_LI-MLN", save = T)
    run("m6", "m6_LI-MLN", save = T)
    run("m7", "m7_LI-MLN", save = T)
    run("m8", "m8_LI-MLN", save = T)
    run("m9", "m9_LI-MLN", save = T)
}
