
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)

idata <- immload(which="not full")

mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)
mid_names <- mid_labels[colnames(mids_counts), "label", drop=T]

overview(immdata=idata)

# make_rarefaction_plots(mids_counts)
make_diversity(mids_counts, mid_names)

## heatmaps
heatmaps(NULL, NULL, cex = 0.5, tsne = T)

## radarcharts
# for (tp in c(0.85)) {
for (tp in c(1, 0.85)) {
  radarcharts("^(Healthy_p5_iga|IBD_p5_iga|Healthy_p10_1_iga)$", "mid2-3-4",  save=T)
  radarcharts("IBD.*iga", "iga_IBD",  save=T)
  radarcharts("SI.*iga", "iga_SI",  save=T)
}

track_clones(NULL,  NULL,  immdata=idata)
track_clones("iga", "iga", immdata=idata)
track_clones("^(Healthy_p5_iga|IBD_p5_iga|Healthy_p10_1_iga)$", "mid2-3-4", immdata=idata)
track_clones("IBD.*iga", "iga_IBD", immdata=idata)
track_clones("SI.*iga", "iga_SI",   immdata=idata)

# library(treemap)
# dirname = make_path("treemaps")
# for (column in colnames(mids_counts)) {
#   mid_df <- cbind(clones=rownames(mids_counts), mids_counts[column])
#   mid_df$freq <- mid_df[,column] / sum(mid_df[,column])
#
#   print(head(mid_df))
#
#   mid_df2 <- get_top(mid_df, topN=10, exclude=c("clones", column))
#   print(head(mid_df2))
#   svg(paste0(dirname, column, "_top10.svg"))
#   treemap(mid_df2, index="clones", vSize="freq")
#   dev.off()
#
#   mid_df2 <- get_top(mid_df, topPerc=0.15, exclude=c("clones", column))
#   print(head(mid_df2))
#   svg(paste0(dirname, column, "_top15perc.svg"))
#   treemap(mid_df2, index="clones", vSize="freq")
#   dev.off()
# }
