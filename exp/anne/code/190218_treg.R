
library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)

idata <- immload(which="full")

mids_counts = read.csv("mids_counts.csv", sep="\t", row.names = 1)
mid_names <- mid_labels[colnames(mids_counts), "label", drop=T]

overview(immdata=idata)
# make_rarefaction_plots(mids_counts, rarefy=T)
make_diversity(mids_counts, mid_names)

## heatmaps

heatmaps(NULL, NULL, cex = 0.5)
heatmaps("[ABCD]", "ABCD") # All without E and O (requested by Anne on 19/03/2019)
heatmaps("[AO]",   "AO") # All only A and O (requested by Milas on 24/06/2021)

## by mouse
for (m_i in 1:5) {
  heatmaps(m_i, paste0("M", m_i))
  # Mouse without E and O (requested by Anne on 19/03/2019)
  heatmaps(paste0(m_i, "[ABCD]"), paste0("M", m_i, "_ABCD"))
}

colours <- brewer.pal(6, "Set1")

for (tp in c(1, 0.85)) {
    radarcharts("1[ABEO]", "M1_ABEO")
    radarcharts("2[ABEO]", "M2_ABEO")
    radarcharts("3[ABEO]", "M3_ABEO")
    radarcharts("4[ABEO]", "M4_ABEO")
    radarcharts("5[ABEO]", "M5_ABEO")

    radarcharts("1[AB]", "M1_AB")
    radarcharts("2[AB]", "M2_AB")
    radarcharts("3[AB]", "M3_AB")
    radarcharts("4[AB]", "M4_AB")
    radarcharts("5[AB]", "M5_AB")

    radarcharts("1[AC]", "M1_AC")
    radarcharts("2[AC]", "M2_AC")
    radarcharts("3[AC]", "M3_AC")
    radarcharts("4[AC]", "M4_AC")
    radarcharts("5[AC]", "M5_AC")

    radarcharts("1[CD]", "M1_CD")
    radarcharts("2[CD]", "M2_CD")
    radarcharts("3[CD]", "M3_CD")
    radarcharts("4[CD]", "M4_CD")
    radarcharts("5[CD]", "M5_CD")

    radarcharts("A", "M12345_A")
    radarcharts("B", "M12345_B")
    radarcharts("C", "M12345_C")
    radarcharts("D", "M12345_D")
    radarcharts("E", "M12345_E")
    radarcharts("O", "M12345_O")
}

track_clones(NULL, NULL, immdata=idata)
