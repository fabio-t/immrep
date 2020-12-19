
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = c("mid31",
                          "mid34", "mid35", "mid36",

                          "mid37",
                          "mid40", "mid41")
f = colwise(get_freqs)
mids_freqs = f(mids_counts)

mid_names = c("TCR1_Ctrl1",
              "TCR1_Col1", "TCR1_Col2", "TCR1_Col3",
              
              "TCR2_Ctrl1",
              "TCR2_Col1", "TCR2_Col2")

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, binary=F)

# TCR1
tcr1_i = grep("TCR1_", mid_names)
make_full_heatmap(mids_counts[tcr1_i], "horn",      mid_names[tcr1_i], prefix="TCR1", binary=F)
make_full_heatmap(mids_counts[tcr1_i], "intersect", mid_names[tcr1_i], prefix="TCR1", range=c(0, NA))
make_full_heatmap(mids_counts[tcr1_i], "intersect", mid_names[tcr1_i], prefix="TCR1_top0.85", range=c(0, NA), topPerc=0.85)

# TCR2
tcr2_i = grep("TCR2_", mid_names)
make_full_heatmap(mids_counts[tcr2_i], "horn",      mid_names[tcr2_i], prefix="TCR2", binary=F)
make_full_heatmap(mids_counts[tcr2_i], "intersect", mid_names[tcr2_i], prefix="TCR2", range=c(0, NA))
make_full_heatmap(mids_counts[tcr2_i], "intersect", mid_names[tcr2_i], prefix="TCR2_top0.85", range=c(0, NA), topPerc=0.85)

# Ctrl
ctrl_i = grep("_Ctrl", mid_names)
make_full_heatmap(mids_counts[ctrl_i], "horn",      mid_names[ctrl_i], prefix="Ctrl", binary=F)
make_full_heatmap(mids_counts[ctrl_i], "intersect", mid_names[ctrl_i], prefix="Ctrl", range=c(0, NA))
make_full_heatmap(mids_counts[ctrl_i], "intersect", mid_names[ctrl_i], prefix="Ctrl_top0.85", range=c(0, NA), topPerc=0.85)

# Colonised
col_i = grep("_Col", mid_names)
make_full_heatmap(mids_counts[col_i], "horn",      mid_names[col_i], prefix="Col", binary=F)
make_full_heatmap(mids_counts[col_i], "intersect", mid_names[col_i], prefix="Col", range=c(0, NA))
make_full_heatmap(mids_counts[col_i], "intersect", mid_names[col_i], prefix="Col_top0.85", range=c(0, NA), topPerc=0.85)

# diversity
divers_df = as.data.frame(mid_names)
rownames(divers_df) = colnames(mids_counts)
colnames(divers_df) = "Sample"
divers_df$Shannon = diversity(t(mids_counts), index="shannon")
divers_df$InvSimpson = diversity(t(mids_counts), index="invsimpson")
write.csv(divers_df, file="diversity.csv")

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

topPerc = 0.85

# 1 shade of blue for one group, 3 shades of red for the other
colours = c(rgb(0, 0, 150, maxColorValue=255), rgb(as.integer(seq(150, 210, length.out=3)), 0, 0, maxColorValue=255))

# TCR1
dirname="TCR1"

column_indices = tcr1_i

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# TCR1
dirname="TCR2"

column_indices = tcr2_i

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

