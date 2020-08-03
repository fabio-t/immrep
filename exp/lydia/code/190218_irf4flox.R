
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c("input1", "input2", "input3",
              "CRE+_m857_LI", "CRE+_m857_MLN", "CRE+_m857_SPL",
              "CRE+_m861_LI", "CRE+_m861_MLN", "CRE+_m861_SPL",
              "CRE+_m862_LI", "CRE+_m862_MLN", "CRE+_m862_SPL",
              "CRE-_m865_LI", "CRE-_m865_MLN", "CRE-_m865_SPL",
              "CRE-_m866_LI", "CRE-_m866_MLN", "CRE-_m866_SPL",
              "CRE-_m868_LI", "CRE-_m868_MLN", "CRE-_m868_SPL"
          )

# make_rarefaction_plots(mids_counts, 1:3, rarefy=T)

## All
make_full_heatmap(mids_counts, "horn",      mid_names, prefix="full", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## LI with inputs
columns_i = grep("(input|_LI)", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="LI-inputs", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="LI-inputs", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="LI-inputs_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## LI without inputs
columns_i = grep("_LI", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="LI", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="LI", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="LI_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## MLN with inputs
columns_i = grep("(input|_MLN)", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="MLN-inputs", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN-inputs", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN-inputs_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## MLN without inputs
columns_i = grep("_MLN", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="MLN", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="MLN_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## SPL with inputs
columns_i = grep("(input|_SPL)", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="SPL-inputs", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="SPL-inputs", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="SPL-inputs_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

## SPL without inputs
columns_i = grep("_SPL", mid_names)
columns_n = mid_names[columns_i]
data = mids_counts[columns_i]
make_full_heatmap(data, "horn",      columns_n, prefix="SPL", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(data, "intersect", columns_n, prefix="SPL", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(data, "intersect", columns_n, prefix="SPL_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# diversity
make_diversity(mids_counts, mid_names)

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours = brewer.pal(6, "Set1")

run <- function(regexp, dirname) {
    print(dirname)

    column_indices = grep(regexp, mid_names)

    column_names   = colnames(mids_counts)[column_indices]
    column_names_2 = mid_names[column_indices]

    print(column_names_2)

    data = mids_counts[column_names]

    save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)
}

for (tp in c(1, 0.85)) {
    run("input",  "inputs")
    run("_LI",  "LI")
    run("_MLN", "MLN")
    run("_SPL", "SPL")
}
