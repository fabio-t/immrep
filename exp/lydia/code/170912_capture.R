
source("../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names <- c(paste0("m1_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")),
               paste0("m2_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")),
               paste0("m4_", c("IFNy-_IL17-", "IFNy+_IL17+",    "Foxp3+")),
               paste0("m5_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")),
               paste0("m6_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")),
               paste0("m7_", c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+"))
           )

# rarefaction
# make_rarefaction_plots(mids_counts, rarefy = T)

heatmaps(NULL, NULL)
# full heatmap without IFNy- IL17-
heatmaps("IFNy-_IL17-", "noIFNy-_IL17-", invert = T)
# requested on 6/02/2020: heatmap with m1, m2, m5, m6, m7 only IFN+, IL-17+ and Foxp3+
# so, essentially, exclude anything with [+-]_IL17 pattern and everything from m4
heatmaps("([\\+-]_IL17|m4_)", "special", invert = T)
# requested on 14/05/2020
heatmaps("^(m1|m2|m5|)_(IFNy\\+|IL17\\+|Foxp3\\+)", "m1-2-5_only-plus")

# diversity
make_diversity(mids_counts, mid_names)

## one-vs-all radarcharts
colours = brewer.pal(9, "Set1")

for (tp in c(1, 0.85)) {
    for (m in c(1,2,4,5,6,7)) {
        radarcharts(paste0("m", m, "_IFNy-_IL17-$"), paste0("m", m, "_noIFNy-_IL17-"), invert = T, save = T)
        radarcharts(paste0("m", m, "_"),  paste0("m", m), save = T)
    }

    radarcharts("IFNy\\+", "IFNy+", save = T)
    radarcharts("IL17\\+", "IL17+", save = T)
    radarcharts("Foxp3", "Foxp3+",  save = T)
}
