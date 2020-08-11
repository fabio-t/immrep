
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = c("mid17", "mid19", "mid20", "mid21", "mid22", "mid23", "mid24", "mid25", "mid26",
                          "mid27", "mid28", "mid29", "mid30", "mid31", "mid32", "mid33", "mid34", "mid35",
                          "mid36", "mid37", "mid38", "mid39", "mid40", "mid41", "mid42", "mid43", "mid44",
                          "mid45", "mid46",          "mid47", "mid48",          "mid50", "mid51", "mid52",
                          "mid53", "mid54", "mid55")
f = colwise(get_freqs)
mids_freqs = f(mids_counts)

mid_names = c("NO_S1-", "NO_S2-", "NO_S3-", "NO_M1-", "NO_M2-", "NO_M3-", "NO_C1-", "NO_C2-", "NO_C3-",
              "NO_S1+", "NO_S2+", "NO_S3+", "NO_M1+", "NO_M2+", "NO_M3+", "NO_C1+", "NO_C2+", "NO_C3+",
              "PF_S1-", "PF_S2-", "PF_S3-", "PF_M1-", "PF_M2-", "PF_M3-", "PF_C1-", "PF_C2-", "PF_C3-",
                        "PF_S2+", "PF_S3+", "PF_M1+", "PF_M2+",           "PF_C1+", "PF_C2+", "PF_C3+",
              "input1", "input2", "input3")

input_columns = c("mid53", "mid54", "mid55")

input1 = rownames(mids_counts)[mids_counts$mid53>0]
input2 = rownames(mids_counts)[mids_counts$mid54>0]
input3 = rownames(mids_counts)[mids_counts$mid55>0]

input_clones = union(union(input1, input2), input3)
common_inputs  = intersect(intersect(input1, input2), input3)
uncommon_inputs = setdiff(input_clones, common_inputs)
make_input_histograms(mids_counts, common_inputs, uncommon_inputs, input_columns)

set.seed(42)
colours = as.character(distinctColorPalette(length(common_inputs)))

# check that the colour "black" is not included in the palette
repeat
{
  if("#000000" %in% colours)
  {
    print("black colour generated in palette, trying again..")
    colours = as.character(distinctColorPalette(length(common_inputs)))
  }
  else
  {
    break
  }
}

# colours = rainbow(length(common_inputs))
clone_colours = data.frame(common_inputs, colours)
rownames(clone_colours) = common_inputs
clone_colours$common_inputs = NULL

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, binary=F)

# NO
make_full_heatmap(mids_counts[1:18], "horn", mid_names[1:18], prefix="NO", binary=F)

# NO FoxP3-
make_full_heatmap(mids_counts[1:9], "horn", mid_names[1:9], prefix="NO_foxp3-", binary=F)
make_full_heatmap(mids_counts[1:9], "intersect", mid_names[1:9], prefix="NO_foxp3-", range=c(0, NA))
make_full_heatmap(mids_counts[1:9], "intersect", mid_names[1:9], prefix="NO_foxp3-_top0.85", range=c(0, NA), topPerc=0.85)
make_abundance_matrix(mids_counts[1:9], topPerc=0.85, prefix="NO_foxp3-_top0.85")

# NO SPLs
make_full_heatmap(mids_counts[c(1:3, 10:12)], "horn", mid_names[c(1:3, 10:12)], prefix="NO_spl", binary=F)

# NO MLNs
make_full_heatmap(mids_counts[c(4:6, 13:15)], "horn", mid_names[c(4:6, 13:15)], prefix="NO_mln", binary=F)

# NO LIs
make_full_heatmap(mids_counts[c(7:9, 16:18)], "horn", mid_names[c(7:9, 16:18)], prefix="NO_li", binary=F)
make_full_heatmap(mids_counts[c(7:9, 16:18)], "intersect", mid_names[c(7:9, 16:18)], prefix="NO_li", range=c(0, NA))
make_full_heatmap(mids_counts[c(7:9, 16:18)], "intersect", mid_names[c(7:9, 16:18)], prefix="NO_li_top0.85", range=c(0, NA), topPerc=0.85)
make_abundance_matrix(mids_counts[c(7:9, 16:18)], topPerc=0.85, prefix="NO_li_top0.85")

# PF
make_full_heatmap(mids_counts[19:34], "horn", mid_names[19:34], prefix="PF", binary=F)

# All FoxP3-
make_full_heatmap(mids_counts[c(1:9, 19:27)], "horn", mid_names[c(1:9, 19:27)], prefix="foxp3-", binary=F)

# All FoxP3+
make_full_heatmap(mids_counts[c(10:18, 28:34)], "horn", mid_names[c(10:18, 28:34)], prefix="foxp3+", binary=F)

# SPLs
make_full_heatmap(mids_counts[c(1:3, 10:12, 19:21, 28, 29)], "horn", mid_names[c(1:3, 10:12, 19:21, 28, 29)], prefix="spl", binary=F)

# SPLs FoxP3-
make_full_heatmap(mids_counts[c(1:3, 19:21)], "horn", mid_names[c(1:3, 19:21)], prefix="spl_foxp3-", binary=F)

# MLNs
make_full_heatmap(mids_counts[c(4:6, 13:15, 22:24, 30, 31)], "horn", mid_names[c(4:6, 13:15, 22:24, 30, 31)], prefix="mln", binary=F)

# LIs
make_full_heatmap(mids_counts[c(7:9, 16:18, 25:27, 32:34)], "horn", mid_names[c(7:9, 16:18, 25:27, 32:34)], prefix="li", binary=F)

# LIs FoxP3-
make_full_heatmap(mids_counts[c(7:9, 25:27)], "horn", mid_names[c(7:9, 25:27)], prefix="li_foxp3-", binary=F)
make_full_heatmap(mids_counts[c(7:9, 25:27)], "intersect", mid_names[c(7:9, 25:27)], prefix="li_foxp3-", range=c(0, NA))
make_full_heatmap(mids_counts[c(7:9, 25:27)], "intersect", mid_names[c(7:9, 25:27)], prefix="li_foxp3-_top0.85", range=c(0, NA), topPerc=0.85)
make_abundance_matrix(mids_counts[c(7:9, 25:27)], topPerc=0.85, prefix="li_foxp3-_top0.85")

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

column_indices = 35:37

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)

# ORGANS

# 3 shades of blue for one group, 3 shades of red for the other
colours = c(rgb(0, 0, as.integer(seq(150, 210, length.out=3)), maxColorValue=255), rgb(as.integer(seq(150, 210, length.out=3)), 0, 0, maxColorValue=255))

# SPL FoxP3-
dirname="spl_foxp3-"

column_indices = c(1:3, 19:21)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# SPL FoxP3+
dirname="spl_foxp3+"

column_indices = c(10:12, 28, 29)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# NO SPL
dirname="NO_spl"

column_indices = c(1:3, 10:12)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# MLN FoxP3-
dirname="mln_foxp3-"

column_indices = c(4:6, 22:24)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# MLN FoxP3+
dirname="mln_foxp3+"

column_indices = c(13:15, 30, 31)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# NO MLN
dirname="NO_mln"

column_indices = c(4:6, 13:15)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# LI FoxP3-
dirname="li_foxp3-"

column_indices = c(7:9, 25:27)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# LI FoxP3+
dirname="li_foxp3+"

column_indices = c(16:18, 32:34)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# NO LI
dirname="NO_li"

column_indices = c(7:9, 16:18)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

## by-mouse comparisons, three colour shades

colours = c(rgb(0, 0, as.integer(seq(150, 210, length.out=2)), maxColorValue=255), rgb(0, as.integer(seq(150, 210, length.out=2)), 0, maxColorValue=255), rgb(as.integer(seq(150, 210, length.out=2)), 0, 0, maxColorValue=255))

# NO 1
dirname="NO_1"

column_indices = seq(1, 18, 3)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# NO 2
dirname="NO_2"

column_indices = seq(2, 18, 3)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# NO 3
dirname="NO_3"

column_indices = seq(3, 18, 3)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)
