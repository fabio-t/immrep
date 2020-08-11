
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = paste0("mid", 40:60)

f = colwise(get_freqs)
mids_freqs = f(mids_counts)

mid_names = c("p1_inflamed_s1", "p1_not-inflamed_s2", "p5_not-inflamed_s1", "p5_inflamed_s2",
              "p10_not-inflamed_s1", "p10_not-inflamed_s2", "p24_inflamed_s1", "p24_not-inflamed_s2", "p24_not-inflamed_s3",
              "p25_not-inflamed_s1", "p25_not-inflamed_s2", "p29_not-inflamed_s1", "p29_inflamed_s2",
              "p44_inflamed_s1", "p44_not-inflamed_s2", "p53_female_s1", "p53_female_s2",
              "p58_not-inflamed_s1", "p58_inflamed_s2", "p59_inflamed_s1", "p59_not-inflamed_s2"
          )

# make_rarefaction_plots(mids_counts, rarefy=T)

# full heatmap
make_full_heatmap(mids_counts, "horn",      mid_names, prefix="full", binary=F, types="square", cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", topPerc=0.85, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

# diversity
divers_df = as.data.frame(mid_names)
rownames(divers_df) = colnames(mids_counts)
colnames(divers_df) = "Sample"
divers_df$Shannon = diversity(t(mids_counts),index="shannon")
divers_df$InvSimpson = diversity(t(mids_counts),index="invsimpson")
write.csv(divers_df, file="diversity.csv")

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours = brewer.pal(5, "Set1")

for (tp in c(1, 0.85))
{
    ## input vs corresponding mouse

    for (p_i in c(1, 5, 10, 24, 25, 29, 44, 53, 58, 59))
    {
        dirname = paste0("p", p_i)

        print(dirname)

        column_indices = grep(paste0("^p", p_i, "_"), mid_names)

        column_names   = colnames(mids_counts)[column_indices]
        column_names_2 = mid_names[column_indices]

        data = mids_counts[column_names]

        make_full_heatmap(data, "horn",      column_names_2, prefix=dirname, binary=F, types="square", cexRow=0.9, cexCol=0.9)
        make_full_heatmap(data, "intersect", column_names_2, prefix=dirname, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)
        make_full_heatmap(data, "intersect", column_names_2, prefix=paste0(dirname, "_top", tp), topPerc=tp, types="square", vlim=c(0, NA), cexRow=0.9, cexCol=0.9)

        save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)
    }
}
