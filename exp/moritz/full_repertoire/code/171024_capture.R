
source("../../../../code/util.R")

# calculate frequencies from counts
mids_counts = read.csv("mids_counts.csv", sep="\t", header=F)
rownames(mids_counts) = mids_counts$V1
mids_counts$V1 = NULL
colnames(mids_counts) = c("mid40", "mid41", "mid42", "mid43", "mid44",
                          "mid45", "mid46", "mid47", "mid48",
                          "mid49", "mid50", "mid51", "mid52",
                          "mid53", "mid54", "mid55", "mid56",
                          "mid57", "mid58", "mid59", "mid60",
                          "mid61", "mid62", "mid63", "mid64"
                      )

f = colwise(get_freqs)
mids_freqs = f(mids_counts)

organs = c("IFNy-_IL17-", "IFNy+", "IL17+", "Foxp3+")

mid_names = c("M1_input", "M2_input", "M3_input", "M4_input", "M5_input",
              paste0("M1_", organs),
              paste0("M2_", organs),
              paste0("M3_", organs),
              paste0("M4_", organs),
              paste0("M5_", organs)
          )

# rarefaction
# make_rarefaction_plots(mids_counts, 1:5)

# full heatmap
make_full_heatmap(mids_counts, "horn",      mid_names, binary=F)
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full", range=c(0, NA))
make_full_heatmap(mids_counts, "intersect", mid_names, prefix="full_top0.85", range=c(0, NA), topPerc=0.85)

# diversity
divers_df = as.data.frame(mid_names)
rownames(divers_df) = colnames(mids_counts)
colnames(divers_df) = "Sample"
divers_df$Shannon = diversity(t(mids_counts),index="shannon")
divers_df$InvSimpson = diversity(t(mids_counts),index="invsimpson")
write.csv(divers_df, file="diversity.csv")

## one-vs-all radarcharts

print("one-vs-all radarcharts..")

colours = brewer.pal(5, "Set1")

for (tp in c(1, 0.85))
{
    ## input vs corresponding mouse

    for (m_i in 1:5)
    {
        dirname = paste0("M", m_i)

        print(dirname)

        column_indices = grep(paste0("M", m_i, "_"), mid_names)

        column_names   = colnames(mids_counts)[column_indices]
        column_names_2 = mid_names[column_indices]

        data = mids_counts[column_names]

        make_full_heatmap(data, "horn",      column_names_2, prefix=dirname, binary=F)
        make_full_heatmap(data, "intersect", column_names_2, prefix=dirname, range=c(0, NA))
        make_full_heatmap(data, "intersect", column_names_2, prefix=paste0(dirname, "_top", tp), range=c(0, NA), topPerc=tp)

        save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)

        ## only the "plus" organs together
        dirname2 = paste0(dirname, "_only-plus")

        print(dirname2)

        column_indices = grep(paste0("^", dirname, "_.*\\+"), mid_names)

        column_names   = colnames(mids_counts)[column_indices]
        column_names_2 = mid_names[column_indices]

        data = mids_counts[column_names]

        make_full_heatmap(data, "horn",      column_names_2, prefix=dirname2, binary=F)
        make_full_heatmap(data, "intersect", column_names_2, prefix=dirname2, range=c(0, NA))
        make_full_heatmap(data, "intersect", column_names_2, prefix=paste0(dirname2, "_top", tp), range=c(0, NA), topPerc=tp)

        save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname2), use_colours=colours, show_axis=F, show_legend=T)
    }

    ## organs by themselves

    for (organ in organs)
    {
        dirname = paste0(organ)

        print(dirname)

        column_indices = grep(paste0("M\\d+_", gsub("+", "\\+", organ, fixed=T), "$"), mid_names)

        column_names   = colnames(mids_counts)[column_indices]
        column_names_2 = mid_names[column_indices]

        data = mids_counts[column_names]

        make_full_heatmap(data, "horn",      column_names_2, prefix=dirname, binary=F)
        make_full_heatmap(data, "intersect", column_names_2, prefix=dirname, range=c(0, NA))
        make_full_heatmap(data, "intersect", column_names_2, prefix=paste0(dirname, "_top", tp), range=c(0, NA), topPerc=tp)

        save_multiple_radarcharts(data, column_names, labels=column_names_2, topPerc=tp, percentage=T, use_log=F, dirname=paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname), use_colours=colours, show_axis=F, show_legend=T)
    }

    ## only plus organs, all mice

    dirname = "only-plus"

    print(dirname)

    column_indices = grep(paste0("^M\\d+_.*\\+"), mid_names)

    column_names   = colnames(mids_counts)[column_indices]
    column_names_2 = mid_names[column_indices]

    data = mids_counts[column_names]

    make_full_heatmap(data, "horn",      column_names_2, prefix=dirname, binary=F)
    make_full_heatmap(data, "intersect", column_names_2, prefix=dirname, range=c(0, NA))
    make_full_heatmap(data, "intersect", column_names_2, prefix=paste0(dirname, "_top", tp), range=c(0, NA), topPerc=tp)
}
