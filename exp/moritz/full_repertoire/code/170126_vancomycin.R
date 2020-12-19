
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = c("mid42", "mid43", "mid44",

                          "mid45", "mid46", "mid47",
                          "mid48", "mid49", "mid50",

                          "mid51", "mid52", "mid53",
                          "mid54", "mid55", "mid56")
f = colwise(get_freqs)
mids_freqs = f(mids_counts)

mid_names = c("input1", "input2", "input3",

              "Ctrl1_LI_FoxP3-",  "Ctrl2_LI_FoxP3-",  "Ctrl3_LI_FoxP3-",
              "Ctrl1_LI_FoxP3+",  "Ctrl2_LI_FoxP3+",  "Ctrl3_LI_FoxP3+",

              "Van1_LI_FoxP3-",  "Van2_LI_FoxP3-",  "Van3_LI_FoxP3-",
              "Van1_LI_FoxP3+",  "Van2_LI_FoxP3+",  "Van3_LI_FoxP3+")

input_columns = c("mid42", "mid43", "mid44")

input1 = rownames(mids_counts)[mids_counts$mid42>0]
input2 = rownames(mids_counts)[mids_counts$mid43>0]
input3 = rownames(mids_counts)[mids_counts$mid44>0]

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

# CLTR
make_full_heatmap(mids_counts[4:9], "horn", mid_names[4:9], prefix="ctrl", binary=F)

# Vancomycin
make_full_heatmap(mids_counts[10:15], "horn", mid_names[10:15], prefix="van", binary=F)

# LIs
make_full_heatmap(mids_counts[4:15], "horn", mid_names[4:15], prefix="li", binary=F)

# LIs FoxP3-
make_full_heatmap(mids_counts[c(4:6, 10:12)], "horn", mid_names[c(4:6, 10:12)], prefix="li_foxp3-", binary=F)

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

# LI FoxP3-
dirname="li_foxp3-"

column_indices = c(4:6, 10:12)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)
# make_chord_diag(data, dirname=paste0("./chords/topPerc", topPerc), prefix=dirname, topPerc=0.85)

# LI FoxP3+
dirname="li_foxp3+"

column_indices = c(7:9, 13:15)

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)
# make_chord_diag(data, dirname=paste0("./chords/topPerc", topPerc), prefix=dirname, topPerc=0.85)

# ctrl LI
dirname="ctrl_li"

column_indices = 4:9

column_names   = colnames(mids_counts)[column_indices]
column_names_2 = mid_names[column_indices]

data = mids_counts[column_names]

save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, show_legend=F)
save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=topPerc, percentage=T, use_log=T, dirname=paste0("./radarcharts/one-vs-all/topPerc", topPerc, "/", dirname), use_colours=colours, join_plots=T)
# make_chord_diag(data, dirname=paste0("./chords/topPerc", topPerc), prefix=dirname, topPerc=0.85)
