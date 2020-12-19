
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = c("mid1",  "mid2",  "mid3",

                          "mid4",  "mid5",  "mid6",  "mid7",
                          "mid8",  "mid9",  "mid10", "mid11",
                          "mid13", "mid14", "mid15", "mid16")
f = colwise(get_freqs)
mids_freqs = f(mids_counts)

mid_names = c("input1", "input2", "input3",

              "M1_IFN-_IL17-", "M1_IFN+",  "M1_IL17+", "M1_FoxP3+",
              "M2_IFN-_IL17-", "M2_IFN+",  "M2_IL17+", "M2_FoxP3+",
              "M3_IFN-_IL17-", "M3_IFN+",  "M3_IL17+", "M3_FoxP3+"
              )

input_columns = c("mid1", "mid2", "mid3")

input1 = rownames(mids_counts)[mids_counts$mid1>0]
input2 = rownames(mids_counts)[mids_counts$mid2>0]
input3 = rownames(mids_counts)[mids_counts$mid3>0]

print("input1")
print(length(input1))
print("input2")
print(length(input2))
print("input3")
print(length(input3))

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, binary=F)

# full heatmap (no input)
make_full_heatmap(mids_counts[-c(1,2,3)], "horn", mid_names[-c(1,2,3)], prefix="full-noinput", binary=F)

# inputs
make_full_heatmap(mids_counts[1:3], "intersect", mid_names[1:3], prefix="inputs", range=c(0, NA))
make_full_heatmap(mids_counts[1:3], "intersect", mid_names[1:3], prefix="inputs_top0.85", range=c(0, NA), topPerc=0.85)

# IFN- IL17-
make_full_heatmap(mids_counts[c(4,8,12)], "horn",      mid_names[c(4,8,12)], prefix="ifn-_il17-", binary=F)
make_full_heatmap(mids_counts[c(4,8,12)], "intersect", mid_names[c(4,8,12)], prefix="ifn-_il17-", range=c(0, NA))
make_full_heatmap(mids_counts[c(4,8,12)], "intersect", mid_names[c(4,8,12)], prefix="ifn-_il17-_top0.85", range=c(0, NA), topPerc=0.85)

# IFN+
make_full_heatmap(mids_counts[c(5,9,13)], "horn",      mid_names[c(5,9,13)], prefix="ifn+", binary=F)
make_full_heatmap(mids_counts[c(5,9,13)], "intersect", mid_names[c(5,9,13)], prefix="ifn+", range=c(0, NA))
make_full_heatmap(mids_counts[c(5,9,13)], "intersect", mid_names[c(5,9,13)], prefix="ifn+_top0.85", range=c(0, NA), topPerc=0.85)

# IL17+
make_full_heatmap(mids_counts[c(6,10,14)], "horn",      mid_names[c(6,10,14)], prefix="il17+", binary=F)
make_full_heatmap(mids_counts[c(6,10,14)], "intersect", mid_names[c(6,10,14)], prefix="il17+", range=c(0, NA))
make_full_heatmap(mids_counts[c(6,10,14)], "intersect", mid_names[c(6,10,14)], prefix="il17+_top0.85", range=c(0, NA), topPerc=0.85)

# FoxP3+
make_full_heatmap(mids_counts[c(7,11,15)], "horn",      mid_names[c(7,11,15)], prefix="foxp3+", binary=F)
make_full_heatmap(mids_counts[c(7,11,15)], "intersect", mid_names[c(7,11,15)], prefix="foxp3+", range=c(0, NA))
make_full_heatmap(mids_counts[c(7,11,15)], "intersect", mid_names[c(7,11,15)], prefix="foxp3+_top0.85", range=c(0, NA), topPerc=0.85)

# M1
make_full_heatmap(mids_counts[4:7], "horn",      mid_names[4:7], prefix="m1", binary=F)
make_full_heatmap(mids_counts[4:7], "intersect", mid_names[4:7], prefix="m1", range=c(0, NA))
make_full_heatmap(mids_counts[4:7], "intersect", mid_names[4:7], prefix="m1_top0.85", range=c(0, NA), topPerc=0.85)

# M2
make_full_heatmap(mids_counts[8:11], "horn",      mid_names[8:11], prefix="m2", binary=F)
make_full_heatmap(mids_counts[8:11], "intersect", mid_names[8:11], prefix="m2", range=c(0, NA))
make_full_heatmap(mids_counts[8:11], "intersect", mid_names[8:11], prefix="m2_top0.85", range=c(0, NA), topPerc=0.85)

# M3
make_full_heatmap(mids_counts[12:15], "horn",      mid_names[12:15], prefix="m3", binary=F)
make_full_heatmap(mids_counts[12:15], "intersect", mid_names[12:15], prefix="m3", range=c(0, NA))
make_full_heatmap(mids_counts[12:15], "intersect", mid_names[12:15], prefix="m3_top0.85", range=c(0, NA), topPerc=0.85)

# rarefaction
# make_rarefaction_plots(mids_counts, input_columns)

# diversity
divers_df = as.data.frame(mid_names)
rownames(divers_df) = colnames(mids_counts)
colnames(divers_df) = "Sample"
divers_df$Shannon = diversity(t(mids_counts),index="shannon")
divers_df$InvSimpson = diversity(t(mids_counts),index="invsimpson")
write.csv(divers_df, file="diversity.csv")

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

# INPUTS

topPerc = 0.85

colours = c(rgb(0.2,0.5,0.5), rgb(0.8,0.2,0), rgb(0.7,0.5,0.1))

dirname="inputs"

column_indices = 2:3

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)

# ORGANS

# 3 shades of blue for one group, 3 shades of red for the other
colours = c(rgb(0, 0, as.integer(seq(150, 210, length.out=3)), maxColorValue=255), rgb(as.integer(seq(150, 210, length.out=3)), 0, 0, maxColorValue=255))

# IFN- IL17-
dirname="ifn-_il17-"

column_indices = c(4,8,12)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# IFN+
dirname="ifn+"

column_indices = c(5,9,13)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# IL17+
dirname="il17+"

column_indices = c(6,10,14)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# FoxP3+
dirname="foxp3+"

column_indices = c(7,11,15)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

## by mouse

# M1
dirname="m1"

column_indices = 4:7

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# M2
dirname="m2"

column_indices = 8:11

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# M3
dirname="m3"

column_indices = 12:15

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)
