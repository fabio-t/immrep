
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)
mids_counts = mids_counts[grep("MID34", colnames(mids_counts), invert=T)] # asked to exclude MID34

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)
mid_names <- mid_labels[colnames(mids_counts), "label", drop=T]

# make_rarefaction_plots(mids_counts, rarefy=T)

## diversity

make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)

## radarcharts

for (tp in c(1, 0.85)) {
    radarcharts("IgA_ID", "IgA_ID",  save=T)
    radarcharts("IgA_SI", "IgA_SI",  save=T)
    radarcharts("IgG_ID", "IgG_ID",  save=T)
    radarcharts("IgG_SI", "IgG_SI",  save=T)
}

library(treemap)
dirname = make_path("treemaps")
for (column in colnames(mids_counts)) {
    mid_df <- cbind(clones=rownames(mids_counts), mids_counts[column])
    mid_df$freq <- mid_df[,column] / sum(mid_df[,column])

    print(head(mid_df))

    mid_df2 <- get_top_clones(mid_df, topN=10, exclude=c("clones", column))
    print(head(mid_df2))
    svg(paste0(dirname, column, "_top10.svg"))
    treemap(mid_df2, index="clones", vSize="freq")
    dev.off()

    mid_df2 <- get_top_clones(mid_df, topPerc=0.15, exclude=c("clones", column))
    print(head(mid_df2))
    svg(paste0(dirname, column, "_top15perc.svg"))
    treemap(mid_df2, index="clones", vSize="freq")
    dev.off()
}
