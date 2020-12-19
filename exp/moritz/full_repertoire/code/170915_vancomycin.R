
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = c("mid17", "mid19", "mid20",

                          "mid21", "mid22", "mid23",
                          "mid24", "mid25",
                          "mid26", "mid27", "mid28",
                          "mid29", "mid30")
f = colwise(get_freqs)
mids_freqs = f(mids_counts)

mid_names = c("input1", "input2", "input3",

              "LI_V1_FoxP3-",    "LI_V2_FoxP3-", "LI_V3_FoxP3-",
              "LI_Ctrl1_FoxP3-", "LI_Ctrl2_FoxP3-",
              "LI_V1_FoxP3+",    "LI_V2_FoxP3+", "LI_V3_FoxP3+",
              "LI_Ctrl1_FoxP3+", "LI_Ctrl2_FoxP3+")

input_columns = c("mid17", "mid19", "mid20")

input1 = rownames(mids_counts)[mids_counts$mid17>0]
input2 = rownames(mids_counts)[mids_counts$mid19>0]
input3 = rownames(mids_counts)[mids_counts$mid20>0]

print("input1")
print(length(input1))
print("input2")
print(length(input2))
print("input3")
print(length(input3))

input_clones = union(union(input1, input2), input3)
common_inputs  = intersect(intersect(input1, input2), input3)
uncommon_inputs = setdiff(input_clones, common_inputs)
make_input_histograms(mids_counts, common_inputs, uncommon_inputs, input_columns)

print(paste("input_clones:", length(input_clones)))
print(paste("common_inputs:", length(common_inputs)))
print(paste("uncommon_inputs:", length(uncommon_inputs)))

# full heatmap
make_full_heatmap(mids_counts, "horn", mid_names, binary=F)

# CLTR
ctrl_i = grep("_Ctrl", mid_names)
make_full_heatmap(mids_counts[ctrl_i], "horn", mid_names[ctrl_i], prefix="ctrl", binary=F)

# Vancomycin
van_i = grep("_V", mid_names)
make_full_heatmap(mids_counts[van_i], "horn", mid_names[van_i], prefix="van", binary=F)

# FoxP3-
foxp3m_i = grep("_FoxP3-", mid_names)
make_full_heatmap(mids_counts[foxp3m_i], "horn", mid_names[foxp3m_i], prefix="foxp3-", binary=F)

# FoxP3+
foxp3p_i = grep("_FoxP3\\+", mid_names)
make_full_heatmap(mids_counts[foxp3p_i], "horn", mid_names[foxp3p_i], prefix="foxp3+", binary=F)

# LI
li_i = grep("LI_", mid_names)
make_full_heatmap(mids_counts[li_i], "horn", mid_names[li_i], prefix="li", binary=F)

# LI FoxP3-
lifoxp3m_i = intersect(li_i, foxp3m_i)
make_full_heatmap(mids_counts[lifoxp3m_i], "horn",      mid_names[lifoxp3m_i], prefix="li_foxp3-", binary=F)
make_full_heatmap(mids_counts[lifoxp3m_i], "intersect", mid_names[lifoxp3m_i], prefix="li_foxp3-", range=c(0, NA))
make_full_heatmap(mids_counts[lifoxp3m_i], "intersect", mid_names[lifoxp3m_i], prefix="li_foxp3-_top0.85", range=c(0, NA), topPerc=0.85)
make_abundance_matrix(mids_counts[lifoxp3m_i], topPerc=0.85, prefix="li_foxp3-_top0.85")

# rarefaction
# make_rarefaction_plots(mids_counts, input_columns)

# diversity
divers_df = as.data.frame(mid_names)
rownames(divers_df) = colnames(mids_counts)
colnames(divers_df) = "Sample"
divers_df$Shannon = diversity(t(mids_counts), index="shannon")
divers_df$InvSimpson = diversity(t(mids_counts), index="invsimpson")
write.csv(divers_df, file="diversity.csv")

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

# INPUTS

topPerc = 0.85

colours = c(rgb(0.2,0.5,0.5), rgb(0.8,0.2,0), rgb(0.7,0.5,0.1))

dirname="inputs"

column_indices = grep("input", mid_names)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)

# ORGANS

# 3 shades of blue for one group, 2 shades of red for the other
colours = c(rgb(0, 0, as.integer(seq(150, 210, length.out=3)), maxColorValue=255), rgb(as.integer(seq(150, 210, length.out=2)), 0, 0, maxColorValue=255))

# LI FoxP3-
dirname="li_foxp3-"

column_indices = lifoxp3m_i

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# LI FoxP3+
dirname="li_foxp3+"

column_indices = intersect(li_i, foxp3p_i)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# ctrl LI
dirname="ctrl_li"

column_indices = intersect(ctrl_i, li_i)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

## by-mouse comparisons, 2 shades of blue for one and 2 shades of red for the other

colours = c(rgb(0, 0, as.integer(seq(150, 210, length.out=2)), maxColorValue=255), rgb(as.integer(seq(150, 210, length.out=2)), 0, 0, maxColorValue=255))

# ctrl
dirname="ctrl"

column_indices = ctrl_i

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)
