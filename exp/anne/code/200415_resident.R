
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)

mid_names = c(
              "Foxp3+HD2+_ILN-L_m1", "Foxp3+HD2-_ILN-L_m1", "Foxp3-HD2+_ILN-L_m1", "Foxp3-HD2-_ILN-L_m1",
			  "Foxp3+HD2+_ILN-R_m1", "Foxp3+HD2-_ILN-R_m1", "Foxp3-HD2+_ILN-R_m1", "Foxp3-HD2-_ILN-R_m1",
              "Foxp3+HD2+_CoMLN_m1", "Foxp3+HD2-_CoMLN_m1", "Foxp3-HD2+_CoMLN_m1", "Foxp3-HD2-_CoMLN_m1",

			  "Foxp3+HD2+_ILN-L_m2", "Foxp3+HD2-_ILN-L_m2", "Foxp3-HD2+_ILN-L_m2", "Foxp3-HD2-_ILN-L_m2",
              "Foxp3+HD2+_ILN-R_m2", "Foxp3+HD2-_ILN-R_m2", "Foxp3-HD2+_ILN-R_m2", "Foxp3-HD2-_ILN-R_m2",
              "Foxp3+HD2+_CoMLN_m2", "Foxp3+HD2-_CoMLN_m2", "Foxp3-HD2+_CoMLN_m2", "Foxp3-HD2-_CoMLN_m2",

			  "Foxp3+HD2+_ILN-L_m3", "Foxp3+HD2-_ILN-L_m3", "Foxp3-HD2+_ILN-L_m3", "Foxp3-HD2-_ILN-L_m3",
              "Foxp3+HD2+_ILN-R_m3", "Foxp3+HD2-_ILN-R_m3", "Foxp3-HD2+_ILN-R_m3", "Foxp3-HD2-_ILN-R_m3",
              "Foxp3+HD2+_CoMLN_m3", "Foxp3+HD2-_CoMLN_m3", "Foxp3-HD2+_CoMLN_m3", "Foxp3-HD2-_CoMLN_m3"
	      )

# make_rarefaction_plots(mids_counts, 1:3, rarefy=T)

print(colnames(mid_names))
length(colnames(mid_names))

make_rarefaction_plots(mids_counts, rarefy = T)

## All
make_full_heatmap(mids_counts, "horn",      mid_names, prefix = "full", binary = F, types = "square", cexRow = 0.9, cexCol = 0.9, tsne = T)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

## by mouse
for (m_i in 1:3) {
    columns_i <- grep(paste0("_m", m_i), mid_names)
    columns_n <- mid_names[columns_i]
    data <- mids_counts[columns_i]
    make_full_heatmap(data, "horn",      columns_n, prefix = paste0("M", m_i), binary = F, types = "square", cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i), types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i, "_top0.85"), topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

    # Foxp3+ HD2+
    columns_i <- grep(paste0("Foxp3\\+HD2\\+_.*_m", m_i), mid_names)
    columns_n <- mid_names[columns_i]
    data <- mids_counts[columns_i]
    make_full_heatmap(data, "horn",      columns_n, prefix = paste0("M", m_i, "_plus-plus"), binary = F, types = "square", cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i, "_plus-plus"), types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i, "_plus-plus_top0.85"), topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)

    # Foxp3- HD2+
    columns_i <- grep(paste0("Foxp3-HD2\\+_.*_m", m_i), mid_names)
    columns_n <- mid_names[columns_i]
    data <- mids_counts[columns_i]
    make_full_heatmap(data, "horn",      columns_n, prefix = paste0("M", m_i, "_minus-plus"), binary = F, types = "square", cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i, "_minus-plus"), types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0("M", m_i, "_minus-plus_top0.85"), topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = 0.9, cexCol = 0.9)
}

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
    run("_ILN-L_m1", "M1_ILN-L")
    run("_ILN-R_m1", "M1_ILN-R")
    run("_CoMLN_m1", "M1_CoMLN")

    run("_ILN-L_m2", "M2_ILN-L")
    run("_ILN-R_m2", "M2_ILN-R")
    run("_CoMLN_m2", "M2_CoMLN")

    run("_ILN-L_m3", "M3_ILN-L")
    run("_ILN-R_m3", "M3_ILN-R")
    run("_CoMLN_m3", "M3_CoMLN")

    run("Foxp3\\+HD2\\+_.*_m1", "M1_plus-plus")
    run("Foxp3\\+HD2-_.*_m1",   "M1_plus-minus")
    run("Foxp3-HD2\\+_.*_m1",   "M1_minus-plus")
    run("Foxp3-HD2-_.*_m1",     "M1_minus-minus")

    run("Foxp3\\+HD2\\+_.*_m2", "M2_plus-plus")
    run("Foxp3\\+HD2-_.*_m2",   "M2_plus-minus")
    run("Foxp3-HD2\\+_.*_m2",   "M2_minus-plus")
    run("Foxp3-HD2-_.*_m2",     "M2_minus-minus")

    run("Foxp3\\+HD2\\+_.*_m3", "M3_plus-plus")
    run("Foxp3\\+HD2-_.*_m3",   "M3_plus-minus")
    run("Foxp3-HD2\\+_.*_m3",   "M3_minus-plus")
    run("Foxp3-HD2-_.*_m3",     "M3_minus-minus")
}
