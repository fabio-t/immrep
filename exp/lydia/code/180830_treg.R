
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

library(matrixStats)

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c("input_r1",

              # LI
              "GFP+_r1_m1_LI",
              "GFP+_r1_m2_LI",
              "GFP+_r1_m3_LI",

              "input_r2_a", "input_r2_b",
              "GFP+_r2_m4_LI",
              "GFP+_r2_m5_LI",
              "GFP+_r2_m6_LI",

              "GFP+_r3_m7_LI",
              "GFP+_r3_m8_LI",
              "GFP+_r3_m9_LI",

              # MLN
              "GFP+_r1_m1_MLN",
              "GFP+_r1_m2_MLN",
              "GFP+_r1_m3_MLN",

              "GFP+_r2_m4_MLN",
              "GFP+_r2_m5_MLN",
              "GFP+_r2_m6_MLN",

              "GFP+_r3_m7_MLN",
              "GFP+_r3_m8_MLN",
              "GFP+_r3_m9_MLN",

              # MLN tra
              "GFP+_r1_m1_MLN_tra",
              "GFP+_r1_m2_MLN_tra",
              "GFP+_r1_m3_MLN_tra",

              "GFP+_r2_m4_MLN_tra",
              "GFP+_r2_m5_MLN_tra",
              "GFP+_r2_m6_MLN_tra",

              "GFP+_r3_m7_MLN_tra",
              "GFP+_r3_m8_MLN_tra",
              "GFP+_r3_m9_MLN_tra"
          )

make_rarefaction_plots(mids_counts, grep("input", mid_names), rarefy=T)

## All
make_full_heatmap(mids_counts, "horn",      mid_names, prefix="full", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## LI plus MLN
columns_i = grep("(MLN$|_LI)", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="LI-MLN-noinput", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="LI-MLN-noinput", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="LI-MLN-noinput_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## no inputs
columns_i = grep("GFP", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="full-noinput", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="full-noinput", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="full-noinput_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## LIs, no inputs
columns_i = grep("(LI|input)", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="LI", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="LI", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="LI_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## LIs plus inputs
columns_i = grep("LI", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="LI-noinput", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="LI-noinput", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="LI-noinput_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## MLN, no inputs
columns_i = grep("(MLN$|input)", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="MLN", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## MLN plus inputs
columns_i = grep("MLN$", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="MLN-noinput", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN-noinput", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN-noinput_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## MLN tra, no inputs
columns_i = grep("(MLN_tra|input)", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="MLN_tra", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN_tra", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN_tra_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## MLN tra plus inputs
columns_i = grep("MLN_tra", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="MLN_tra-noinput", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN_tra-noinput", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN_tra-noinput_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

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
    run("(r1_.*_LI|input_r1)", "r1_LI_input")
    run("(r2_.*_LI|input_r2)", "r2_LI_input")

    run("r1_.*_LI", "r1_LI", save = T)
    run("r2_.*_LI", "r2_LI", save = T)
    run("r3_.*_LI", "r3_LI", save = T)

    run("(r1_.*_MLN$|input_r1)", "r1_MLN_input")
    run("(r2_.*_MLN$|input_r2)", "r2_MLN_input")

    run("r1_.*_MLN$", "r1_MLN", save = T)
    run("r2_.*_MLN$", "r2_MLN", save = T)
    run("r3_.*_MLN$", "r3_MLN", save = T)

    run("(r1_.*_MLN_tra|input_r1)", "r1_MLN_tra_input")
    run("(r2_.*_MLN_tra|input_r2)", "r2_MLN_tra_input")

    run("r1_.*_MLN_tra", "r1_MLN_tra", save = T)
    run("r2_.*_MLN_tra", "r2_MLN_tra", save = T)
    run("r3_.*_MLN_tra", "r3_MLN_tra", save = T)

    run("r1_.*_(LI|MLN$)", "r1_LI-MLN", save = T)
    run("r2_.*_(LI|MLN$)", "r2_LI-MLN", save = T)
    run("r3_.*_(LI|MLN$)", "r3_LI-MLN", save = T)

    run("m1_(LI|MLN$)", "m1_LI-MLN", save = T)
    run("m2_(LI|MLN$)", "m2_LI-MLN", save = T)
    run("m3_(LI|MLN$)", "m3_LI-MLN", save = T)
    run("m4_(LI|MLN$)", "m4_LI-MLN", save = T)
    run("m5_(LI|MLN$)", "m5_LI-MLN", save = T)
    run("m6_(LI|MLN$)", "m6_LI-MLN", save = T)
    run("m7_(LI|MLN$)", "m7_LI-MLN", save = T)
    run("m8_(LI|MLN$)", "m8_LI-MLN", save = T)
    run("m9_(LI|MLN$)", "m9_LI-MLN", save = T)

    run("r1_m[12]_LI", "r1_m1-m2_LI", save = T)
    run("r1_m[13]_LI", "r1_m1-m3_LI", save = T)
    run("r1_m[23]_LI", "r1_m2-m3_LI", save = T)

    run("r1_m[12]_MLN$", "r1_m1-m2_MLN", save = T)
    run("r1_m[13]_MLN$", "r1_m1-m3_MLN", save = T)
    run("r1_m[23]_MLN$", "r1_m2-m3_MLN", save = T)

    run("r1_m1_(LI|MLN$)", "r1_m1_LI-MLN", save = T)
    run("r1_m2_(LI|MLN$)", "r1_m2_LI-MLN", save = T)
    run("r1_m3_(LI|MLN$)", "r1_m3_LI-MLN", save = T)

    run("r2_m[45]_LI", "r2_m4-m5_LI", save = T)
    run("r2_m[46]_LI", "r2_m4-m6_LI", save = T)
    run("r2_m[56]_LI", "r2_m5-m6_LI", save = T)

    run("r2_m[45]_MLN$", "r2_m4-m5_MLN", save = T)
    run("r2_m[46]_MLN$", "r2_m4-m6_MLN", save = T)
    run("r2_m[56]_MLN$", "r2_m5-m6_MLN", save = T)

    run("r2_m4_(LI|MLN$)", "r2_m4_LI-MLN", save = T)
    run("r2_m5_(LI|MLN$)", "r2_m5_LI-MLN", save = T)
    run("r2_m6_(LI|MLN$)", "r2_m6_LI-MLN", save = T)

    run("r3_m[78]_LI", "r3_m7-m8_LI", save = T)
    run("r3_m[79]_LI", "r3_m7-m9_LI", save = T)
    run("r3_m[89]_LI", "r3_m8-m9_LI", save = T)

    run("r3_m[78]_MLN$", "r3_m7-m8_MLN", save = T)
    run("r3_m[79]_MLN$", "r3_m7-m9_MLN", save = T)
    run("r3_m[89]_MLN$", "r3_m8-m9_MLN", save = T)

    run("r3_m7_(LI|MLN$)", "r3_m7_LI-MLN", save = T)
    run("r3_m8_(LI|MLN$)", "r3_m8_LI-MLN", save = T)
    run("r3_m9_(LI|MLN$)", "r3_m9_LI-MLN", save = T)
}
