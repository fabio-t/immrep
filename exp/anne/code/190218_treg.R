
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts <- read.csv("mids_counts.csv", sep = "\t", header = F)
rownames(mids_counts) <- mids_counts$V1
mids_counts$V1 <- NULL
print(dim(mids_counts))
colnames(mids_counts) <- c(25:39, 5:9, 45:49, 40, 21, 55, 56, 59)

print(colnames(mids_counts))
length(colnames(mids_counts))

f <- colwise(get_freqs)
mids_freqs <- f(mids_counts)

mid_names <- paste0(1:5, rep(c("A", "B", "C", "D", "O", "E"), each = 5))

print(colnames(mid_names))
length(colnames(mid_names))

make_rarefaction_plots(mids_counts, rarefy = T)

## All
make_full_heatmap(mids_counts, "horn", mid_names, prefix = "full", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

# All without E and O (requested by Anne on 19/03/2019)
columns_i <- grep("[ABCD]", mid_names)
columns_n <- mid_names[columns_i]
data <- mids_counts[columns_i]
make_full_heatmap(data, "horn", columns_n, prefix = "full_noEO", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(data, "intersect", columns_n, prefix = "full_noEO", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(data, "intersect", columns_n, prefix = "full_noEO_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

## by mouse
for (m_i in 1:5) {
    columns_i <- grep(m_i, mid_names)
    columns_n <- mid_names[columns_i]
    data <- mids_counts[columns_i]
    make_full_heatmap(data, "horn", columns_n, prefix = paste0("M", m_i), binary = F, types = "square", cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i), types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i, "_top0.85"), topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

    # Mouse without E and O (requested by Anne on 19/03/2019)
    columns_i <- grep(paste0(m_i, "[ABCD]"), mid_names)
    columns_n <- mid_names[columns_i]
    data <- mids_counts[columns_i]
    make_full_heatmap(data, "horn", columns_n, prefix = paste0("M", m_i, "_noEO"), binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i, "_noEO"), types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i, "_noEO_top0.85"), topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
}

# diversity
divers_df <- as.data.frame(mid_names)
rownames(divers_df) <- colnames(mids_counts)
colnames(divers_df) <- "Sample"
divers_df$Shannon <- diversity(t(mids_counts), index = "shannon")
divers_df$InvSimpson <- diversity(t(mids_counts), index = "invsimpson")
write.csv(divers_df, file = "diversity.csv")

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours <- brewer.pal(6, "Set1")

run <- function(regexp, dirname) {
    print(dirname)

    column_indices <- grep(regexp, mid_names)

    column_names <- colnames(mids_counts)[column_indices]
    column_names_2 <- mid_names[column_indices]

    print(column_names_2)

    data <- mids_counts[column_names]

    save_multiple_radarcharts(data, column_names, labels = column_names_2, topPerc = tp, percentage = T, use_log = F, dirname = paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours = colours, show_axis = F, show_legend = T)
}

for (tp in c(1, 0.85)) {
    run("1[ABEO]", "M1_ABEO")
    run("2[ABEO]", "M2_ABEO")
    run("3[ABEO]", "M3_ABEO")
    run("4[ABEO]", "M4_ABEO")
    run("5[ABEO]", "M5_ABEO")

    run("1[AB]", "M1_AB")
    run("2[AB]", "M2_AB")
    run("3[AB]", "M3_AB")
    run("4[AB]", "M4_AB")
    run("5[AB]", "M5_AB")

    run("1[AC]", "M1_AC")
    run("2[AC]", "M2_AC")
    run("3[AC]", "M3_AC")
    run("4[AC]", "M4_AC")
    run("5[AC]", "M5_AC")

    run("1[CD]", "M1_CD")
    run("2[CD]", "M2_CD")
    run("3[CD]", "M3_CD")
    run("4[CD]", "M4_CD")
    run("5[CD]", "M5_CD")

    run("A", "M12345_A")
    run("B", "M12345_B")
    run("C", "M12345_C")
    run("D", "M12345_D")
    run("E", "M12345_E")
    run("O", "M12345_O")
}
