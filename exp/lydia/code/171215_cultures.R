
source("../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c("Pool_1.1_a", "Pool_1.2_a", "Pool_1.3_a",
              "Pool_2.1_a", "Pool_2.2_a", "Pool_2.3_a",
              "Pool_3.1_a", "Pool_3.2_a", "Pool_3.3_a",
              "Pool_4.1_a", "Pool_4.2_a", "Pool_4.3_a",
              "Pool_5.1_a", "Pool_5.2_a", "Pool_5.3_a",

              "Pool_1.1_b", "Pool_1.2_b", "Pool_1.3_b",
              "Pool_2.1_b", "Pool_2.2_b", "Pool_2.3_b",
              "Pool_3.1_b", "Pool_3.2_b", "Pool_3.3_b",
              "Pool_4.1_b", "Pool_4.2_b", "Pool_4.3_b",
              "Pool_5.1_b", "Pool_5.2_b", "Pool_5.3_b")

# rarefaction
make_rarefaction_plots(mids_counts, rarefy=T)

# full heatmap
make_full_heatmap(mids_counts, "horn",      mid_names, prefix="full", binary=F, types="square", cexRow=0.9, cexCol=0.9, tsne=T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# full heatmap - A
make_full_heatmap(mids_counts[1:15], "horn",      mid_names[1:15], prefix="full-A", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[1:15], "intersect", mid_names[1:15], prefix="full-A", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[1:15], "intersect", mid_names[1:15], prefix="full-A_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# full heatmap - B
make_full_heatmap(mids_counts[16:30], "horn",      mid_names[16:30], prefix="full-B", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[16:30], "intersect", mid_names[16:30], prefix="full-B", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts[16:30], "intersect", mid_names[16:30], prefix="full-B_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# diversity
make_diversity(mids_counts, mid_names)

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

## By Pool

colours = c(rgb(0.2,0.5,0.5), rgb(0.8,0.2,0), rgb(0.7,0.5,0.1))

for (tp in c(1, 0.85))
{
  for (c in c("a", "b"))
  {
    for (i in 1:5)
    {
      dirname = paste0("pool_", i, "_", c)

      print(dirname)

      column_indices = grep(paste0("Pool_", i, ".[123]_", c), mid_names)

      column_names   = colnames(mids_counts)[column_indices]
      column_names_2 = mid_names[column_indices]

      data = mids_counts[column_names]

      make_full_heatmap(data, "horn",      column_names_2, prefix=dirname, binary=F, types="square", cexRow=0.9, cexCol=0.9)
      make_full_heatmap(data, "intersect", column_names_2, prefix=dirname, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
      make_full_heatmap(data, "intersect", column_names_2, prefix=paste0(dirname, "_top", tp), topPerc=tp, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

      save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)
    }
  }
}

## Cross-pools

# 3 shades of blue for one group, 3 shades of red for the other
colours = c(rgb(0, 0, as.integer(seq(150, 210, length.out=3)), maxColorValue=255), rgb(as.integer(seq(150, 210, length.out=3)), 0, 0, maxColorValue=255))

for (tp in c(1, 0.85))
{
  for (i in 1:5)
  {
    dirname = paste0("pool_", i)

    print(dirname)

    column_indices = grep(paste0("Pool_", i, ".[123]_[ab]"), mid_names)

    column_names   = colnames(mids_counts)[column_indices]
    column_names_2 = mid_names[column_indices]

    data = mids_counts[column_names]

    make_full_heatmap(data, "horn",      column_names_2, prefix=dirname, binary=F, types="square", cexRow=0.9, cexCol=0.9)
    make_full_heatmap(data, "intersect", column_names_2, prefix=dirname, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
    make_full_heatmap(data, "intersect", column_names_2, prefix=paste0(dirname, "_top", tp), topPerc=tp, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

    save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)
  }
}
