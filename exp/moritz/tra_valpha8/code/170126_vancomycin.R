
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = c("mid1",  "mid2",  "mid3",

                          "mid4",  "mid5",  "mid6",
                          "mid7",  "mid8",  "mid9",
                          "mid10", "mid11", "mid13",

                          "mid14", "mid15", "mid16",
                          "mid17", "mid19", "mid20",
                          "mid21", "mid22", "mid23",

                          "mid24", "mid25", "mid26",
                          "mid27", "mid28", "mid29",
                          "mid30", "mid31", "mid32",

                          "mid33", "mid34", "mid35",
                          "mid36", "mid37", "mid38",
                          "mid39", "mid40", "mid41")
f = colwise(get_freqs)
mids_freqs = f(mids_counts)

mid_names = c("input1", "input2", "input3",

              "Ctrl1_MLN_FoxP3-", "Ctrl2_MLN_FoxP3-", "Ctrl3_MLN_FoxP3-",
              "Ctrl1_SPL_FoxP3-", "Ctrl2_SPL_FoxP3-", "Ctrl3_SPL_FoxP3-",
              "Ctrl1_LI_FoxP3-",  "Ctrl2_LI_FoxP3-",  "Ctrl3_LI_FoxP3-",

              "Ctrl1_MLN_FoxP3+", "Ctrl2_MLN_FoxP3+", "Ctrl3_MLN_FoxP3+",
              "Ctrl1_SPL_FoxP3+", "Ctrl2_SPL_FoxP3+", "Ctrl3_SPL_FoxP3+",
              "Ctrl1_LI_FoxP3+",  "Ctrl2_LI_FoxP3+",  "Ctrl3_LI_FoxP3+",

              "Van1_MLN_FoxP3-", "Van2_MLN_FoxP3-", "Van3_MLN_FoxP3-",
              "Van1_SPL_FoxP3-", "Van2_SPL_FoxP3-", "Van3_SPL_FoxP3-",
              "Van1_LI_FoxP3-",  "Van2_LI_FoxP3-",  "Van3_LI_FoxP3-",

              "Van1_MLN_FoxP3+", "Van2_MLN_FoxP3+", "Van3_MLN_FoxP3+",
              "Van1_SPL_FoxP3+", "Van2_SPL_FoxP3+", "Van3_SPL_FoxP3+",
              "Van1_LI_FoxP3+",  "Van2_LI_FoxP3+",  "Van3_LI_FoxP3+")

input_columns = c("mid1", "mid2", "mid3")

input1 = rownames(mids_counts)[mids_counts$mid1>0]
input2 = rownames(mids_counts)[mids_counts$mid2>0]
input3 = rownames(mids_counts)[mids_counts$mid3>0]

input_clones = union(union(input1, input2), input3)
common_inputs  = intersect(intersect(input1, input2), input3)
uncommon_inputs = setdiff(input_clones, common_inputs)
make_input_histograms(mids_counts, common_inputs, uncommon_inputs, input_columns)

print(paste("input_clones:", length(input_clones)))
print(paste("common_inputs:", length(common_inputs)))
print(paste("uncommon_inputs:", length(uncommon_inputs)))

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

# Ctrl
ctrl_i = grep("Ctrl", mid_names)
make_full_heatmap(mids_counts[ctrl_i], "horn", mid_names[ctrl_i], prefix="ctrl", binary=F)

# Vancomycin
van_i = grep("Van", mid_names)
make_full_heatmap(mids_counts[ctrl_i], "horn", mid_names[ctrl_i], prefix="van", binary=F)

# FoxP3-
foxp3m_i = grep("FoxP3-", mid_names)
make_full_heatmap(mids_counts[foxp3m_i], "horn", mid_names[foxp3m_i], prefix="foxp3-", binary=F)

# FoxP3+
foxp3p_i = grep("FoxP3\\+", mid_names)
make_full_heatmap(mids_counts[foxp3p_i], "horn", mid_names[foxp3p_i], prefix="foxp3+", binary=F)

# MLN
mln_i = grep("MLN", mid_names)
make_full_heatmap(mids_counts[mln_i], "horn", mid_names[mln_i], prefix="mln", binary=F)

# SPL
spl_i = grep("SPL", mid_names)
make_full_heatmap(mids_counts[spl_i], "horn", mid_names[spl_i], prefix="spl", binary=F)

# SPL FoxP3-
splfoxp3m_i = intersect(spl_i, foxp3m_i)
make_full_heatmap(mids_counts[splfoxp3m_i], "horn", mid_names[splfoxp3m_i], prefix="spl_foxp3-", binary=F)

# LI
li_i = grep("LI", mid_names)
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
divers_df$Shannon = diversity(t(mids_counts),index="shannon")
divers_df$InvSimpson = diversity(t(mids_counts),index="invsimpson")
write.csv(divers_df, file="diversity.csv")

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

# INPUTS

topPerc = 0.85

colours = c(rgb(0.2,0.5,0.5), rgb(0.8,0.2,0), rgb(0.7,0.5,0.1))

dirname="inputs"

column_indices = 1:3

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)

# ORGANS

# 3 shades of blue for one group, 3 shades of red for the other
colours = c(rgb(0, 0, as.integer(seq(150, 210, length.out=3)), maxColorValue=255), rgb(as.integer(seq(150, 210, length.out=3)), 0, 0, maxColorValue=255))

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

column_indices = c(13:15, 31:33)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# ctrl MLN
dirname="ctrl_mln"

column_indices = c(4:6, 13:15)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# SPL FoxP3-
dirname="spl_foxp3-"

column_indices = c(7:9, 25:27)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# SPL FoxP3+
dirname="spl_foxp3+"

column_indices = c(16:18, 34:36)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# ctrl SPL
dirname="ctrl_spl"

column_indices = c(7:9, 16:18)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# LI FoxP3-
dirname="li_foxp3-"

column_indices = c(10:12, 28:30)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# LI FoxP3+
dirname="li_foxp3+"

column_indices = c(19:21, 37:39)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# ctrl LI
dirname="ctrl_li"

column_indices = c(10:12, 19:21)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

## by-mouse comparisons

colours = c(rgb(0.2,0.5,0.5), rgb(0.8,0.2,0), rgb(0.7,0.5,0.1))

# ctrl1
dirname="ctrl_1"

column_indices = grep("Ctrl1_(MLN|SPL|LI)_FoxP3-", mid_names)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# ctrl2
dirname="ctrl_2"

column_indices = grep("Ctrl2_(MLN|SPL|LI)_FoxP3-", mid_names)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# ctrl3
dirname="ctrl_3"

column_indices = grep("Ctrl3_(MLN|SPL|LI)_FoxP3-", mid_names)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# Vanco1
dirname="van_1"

column_indices = grep("Van1_(MLN|SPL|LI)_FoxP3-", mid_names)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# Vanco2
dirname="van_2"

column_indices = grep("Van2_(MLN|SPL|LI)_FoxP3-", mid_names)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)

# Vanco3
dirname="van_3"

column_indices = grep("Van3_(MLN|SPL|LI)_FoxP3-", mid_names)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)
