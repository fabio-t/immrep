
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts_base   = read.csv("../180830_sort/mids_counts.csv", sep="\t", row.names = 1)
colnames(mids_counts_base) <- paste0("base_", colnames(mids_counts_base))
mids_counts_base$rows <- rownames(mids_counts_base)
mids_counts_noUMIs = read.csv("../180830_sort_noUMIs/mids_counts.csv", sep="\t", row.names = 1)
colnames(mids_counts_noUMIs) <- paste0("noUMIs_", colnames(mids_counts_noUMIs))
mids_counts_noUMIs$rows <- rownames(mids_counts_noUMIs)
mids_counts_strict = read.csv("../180830_sort_strict/mids_counts.csv", sep="\t", row.names = 1)
colnames(mids_counts_strict) <- paste0("strict_", colnames(mids_counts_strict))
mids_counts_strict$rows <- rownames(mids_counts_strict)

mid_names = c("10k_sorted_a",
              "10k_sorted_b",
              "10k_sorted_c",
              "10k_sorted_d",
              "10k_sorted_e",

              "10k_seeding_a",
              "10k_seeding_b",
              "10k_seeding_c",
              "10k_seeding_d",
              "10k_seeding_e"
          )

mid_names <- c(paste0("base_", mid_names), paste0("noUMIs_", mid_names), paste0("strict_", mid_names))

# joining all dataframes
mids_counts <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "rows", all = TRUE), list(mids_counts_base, mids_counts_noUMIs, mids_counts_strict))
rownames(mids_counts) <- mids_counts$rows
mids_counts$rows = NULL

mids_counts[is.na(mids_counts)] <- 0

write.csv(mids_counts, file="mids_counts.csv")

## heatmaps
heatmaps(NULL, NULL, cex = 0.5)

# diversity
make_diversity(mids_counts, mid_names)

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours = brewer.pal(6, "Set1")

for (tp in c(1, 0.85)) {
    radarcharts("sorted_a", "sorted_a", save = T)
    radarcharts("sorted_b", "sorted_b", save = T)
    radarcharts("sorted_c", "sorted_c", save = T)
    radarcharts("sorted_d", "sorted_d", save = T)
    radarcharts("sorted_e", "sorted_e", save = T)
    
    radarcharts("sorted_[ab]", "sorted_ab", save = T)
    radarcharts("sorted_[cd]", "sorted_cd", save = T)

    radarcharts("seeding_a", "seeding_a", save = T)
    radarcharts("seeding_b", "seeding_b", save = T)
    radarcharts("seeding_c", "seeding_c", save = T)
    radarcharts("seeding_d", "seeding_d", save = T)
    radarcharts("seeding_e", "seeding_e", save = T)
    
    radarcharts("seeding_[ab]", "seeding_ab", save = T)
    radarcharts("seeding_[cd]", "seeding_cd", save = T)
}
