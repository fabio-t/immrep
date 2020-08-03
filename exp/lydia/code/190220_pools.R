
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

library(matrixStats)

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c("iTreg_pool1", "NegTGFb_pool1", "Neg_pool1",
              "iTreg_pool2", "NegTGFb_pool2", "Neg_pool2",
              "iTreg_pool3", "NegTGFb_pool3", "Neg_pool3"
          )

# rarefaction
# make_rarefaction_plots(mids_counts, rarefy=T)

# full heatmap
make_full_heatmap(mids_counts, "horn",      mid_names, prefix="full", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# pool 1
make_full_heatmap(mids_counts[1:3], "horn",      mid_names[1:3], prefix="pool1", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[1:3], "intersect", mid_names[1:3], prefix="pool1", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[1:3], "intersect", mid_names[1:3], prefix="pool1_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# pool 2
make_full_heatmap(mids_counts[4:6], "horn",      mid_names[4:6], prefix="pool2", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[4:6], "intersect", mid_names[4:6], prefix="pool2", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[4:6], "intersect", mid_names[4:6], prefix="pool2_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# pool 3
make_full_heatmap(mids_counts[7:9], "horn",      mid_names[7:9], prefix="pool3", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[7:9], "intersect", mid_names[7:9], prefix="pool3", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[7:9], "intersect", mid_names[7:9], prefix="pool3_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# diversity
make_diversity(mids_counts, mid_names)

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours = c(rgb(0.2,0.5,0.5), rgb(0.8,0.2,0), rgb(0.7,0.5,0.1))

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
        colnames(d3) <- column_names_2
        rownames(d3) <- l$clones

        write.csv(d1, file = paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname, "/data_perc.csv"))
        write.csv(d2, file = paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname, "/data_counts.csv"))
        write.csv(d3, file = paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname, "/data_ranks.csv"))
    }
}

for (tp in c(1, 0.85)) {
    run("pool1", "pool1", save = T)
    run("pool2", "pool2", save = T)
    run("pool3", "pool3", save = T)

    run("iTreg", "iTreg", save = T)

    run("(iTreg|NegTGFb)_pool1", "pool1_12", save = T)
    run("(iTreg|Neg)_pool1",     "pool1_13", save = T)
    run("(NegTGFb|Neg)_pool1",   "pool1_23", save = T)

    run("(iTreg|NegTGFb)_pool2", "pool2_12", save = T)
    run("(iTreg|Neg)_pool2",     "pool2_13", save = T)
    run("(NegTGFb|Neg)_pool2",   "pool2_23", save = T)

    run("(iTreg|NegTGFb)_pool3", "pool3_12", save = T)
    run("(iTreg|Neg)_pool3",     "pool3_13", save = T)
    run("(NegTGFb|Neg)_pool3",   "pool3_23", save = T)
}
