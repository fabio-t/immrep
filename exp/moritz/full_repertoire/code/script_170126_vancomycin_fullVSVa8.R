
source("../../../../code/util.R")

# Va8 extracted from full repertoire: vancomycin
mids_counts_full = read.csv("../170126_ql15_onlyVa8/mids_counts.csv", sep="\t", header=F)
rownames(mids_counts_full) = mids_counts_full$V1
mids_counts_full$V1 = NULL
colnames(mids_counts_full) = c("input1", "input2", "input3",

                               "Ctrl1_LI_FoxP3-",  "Ctrl2_LI_FoxP3-",  "Ctrl3_LI_FoxP3-",
                               "Ctrl1_LI_FoxP3+",  "Ctrl2_LI_FoxP3+",  "Ctrl3_LI_FoxP3+",

                               "Van1_LI_FoxP3-",   "Van2_LI_FoxP3-",   "Van3_LI_FoxP3-",
                               "Van1_LI_FoxP3+",   "Van2_LI_FoxP3+",   "Van3_LI_FoxP3+")

# Va8-only experiment: vancomycin
mids_counts_va8 = read.csv("../../tra_valpha8/170126_ql15/mids_counts.csv", sep="\t", header=F)
rownames(mids_counts_va8) = mids_counts_va8$V1
mids_counts_va8$V1 = NULL
colnames(mids_counts_va8) = c("input1", "input2", "input3",

                              "Ctrl1_MLN_FoxP3-", "Ctrl2_MLN_FoxP3-", "Ctrl3_MLN_FoxP3-",
                              "Ctrl1_SPL_FoxP3-", "Ctrl2_SPL_FoxP3-", "Ctrl3_SPL_FoxP3-",
                              "Ctrl1_LI_FoxP3-",  "Ctrl2_LI_FoxP3-",  "Ctrl3_LI_FoxP3-",
                
                              "Ctrl1_MLN_FoxP3+", "Ctrl2_MLN_FoxP3+", "Ctrl3_MLN_FoxP3+",
                              "Ctrl1_SPL_FoxP3+", "Ctrl2_SPL_FoxP3+", "Ctrl3_SPL_FoxP3+",
                              "Ctrl1_LI_FoxP3+",  "Ctrl2_LI_FoxP3+",  "Ctrl3_LI_FoxP3+",
                
                              "Van1_MLN_FoxP3-",  "Van2_MLN_FoxP3-",  "Van3_MLN_FoxP3-",
                              "Van1_SPL_FoxP3-",  "Van2_SPL_FoxP3-",  "Van3_SPL_FoxP3-",
                              "Van1_LI_FoxP3-",   "Van2_LI_FoxP3-",   "Van3_LI_FoxP3-",
                
                              "Van1_MLN_FoxP3+",  "Van2_MLN_FoxP3+",  "Van3_MLN_FoxP3+",
                              "Van1_SPL_FoxP3+",  "Van2_SPL_FoxP3+",  "Van3_SPL_FoxP3+",
                              "Van1_LI_FoxP3+",   "Van2_LI_FoxP3+",   "Van3_LI_FoxP3+")

# Va8-only experiment: ampicillin (for "control")
mids_counts_amp = read.csv("../../tra_valpha8/160203_ql15/mids_counts.csv", sep="\t", header=F)
rownames(mids_counts_amp) = mids_counts_amp$V1
mids_counts_amp$V1 = NULL
colnames(mids_counts_amp) = c("Ctrl1_SPL_FoxP3-", "Ctrl1_SPL_FoxP3+", "Ctrl2_SPL_FoxP3-",
                              "Ctrl2_SPL_FoxP3+", "Ctrl3_SPL_FoxP3-", "Ctrl3_SPL_FoxP3+",
                              "Ctrl1_MLN_FoxP3-", "Ctrl1_MLN_FoxP3+", "Ctrl2_MLN_FoxP3-",
                              "Ctrl2_MLN_FoxP3+", "Ctrl3_MLN_FoxP3-", "Ctrl3_MLN_FoxP3+",
                              "Ctrl1_LI_FoxP3-",  "Ctrl1_LI_FoxP3+",  "Ctrl2_LI_FoxP3-",
                              "Ctrl2_LI_FoxP3+",  "Ctrl3_LI_FoxP3-",  "Ctrl3_LI_FoxP3+",
                              "Amp1_SPL_FoxP3-",  "Amp1_SPL_FoxP3+",  "Amp2_SPL_FoxP3-",
                              "Amp2_SPL_FoxP3+",  "Amp3_SPL_FoxP3-",  "Amp3_SPL_FoxP3+",
                              "Amp1_MLN_FoxP3-",  "Amp1_MLN_FoxP3+",  "Amp2_MLN_FoxP3-",
                              "Amp2_MLN_FoxP3+",  "Amp3_MLN_FoxP3-",  "Amp3_MLN_FoxP3+",
                              "Amp1_LI_FoxP3-",   "Amp1_LI_FoxP3+",   "Amp2_LI_FoxP3-",
                              "Amp2_LI_FoxP3+",   "Amp3_LI_FoxP3-",   "Amp3_LI_FoxP3+",
                              "input1",           "input2",           "input3")


all_clones = union(union(rownames(mids_counts_full), rownames(mids_counts_va8)), rownames(mids_counts_amp))

print("all clones:")
print(length(all_clones))

mids_counts = data.frame(full_input1=rep(0, length(all_clones)), full_input2=rep(0, length(all_clones)), full_input3=rep(0, length(all_clones)),
                         va8_input1=rep(0, length(all_clones)),  va8_input2=rep(0, length(all_clones)),  va8_input3=rep(0, length(all_clones)),
                         amp_input1=rep(0, length(all_clones)),  amp_input2=rep(0, length(all_clones)),  amp_input3=rep(0, length(all_clones)),
                         full_Ctrl1_LI_FoxP3_minus=rep(0, length(all_clones)), full_Ctrl2_LI_FoxP3_minus=rep(0, length(all_clones)), full_Ctrl3_LI_FoxP3_minus=rep(0, length(all_clones)),
                         va8_Ctrl1_LI_FoxP3_minus=rep(0, length(all_clones)), va8_Ctrl2_LI_FoxP3_minus=rep(0, length(all_clones)), va8_Ctrl3_LI_FoxP3_minus=rep(0, length(all_clones)),
                         amp_Ctrl1_LI_FoxP3_minus=rep(0, length(all_clones)), amp_Ctrl2_LI_FoxP3_minus=rep(0, length(all_clones)), amp_Ctrl3_LI_FoxP3_minus=rep(0, length(all_clones)),
                         row.names=all_clones)

mids_counts[rownames(mids_counts_full), c("full_input1", "full_input2", "full_input3")] = mids_counts_full[c("input1", "input2", "input3")]
mids_counts[rownames(mids_counts_va8),  c("va8_input1",  "va8_input2",  "va8_input3")]  = mids_counts_va8[c("input1", "input2", "input3")]
mids_counts[rownames(mids_counts_amp),  c("amp_input1",  "amp_input2",  "amp_input3")]  = mids_counts_amp[c("input1", "input2", "input3")]
mids_counts[rownames(mids_counts_full), c("full_Ctrl1_LI_FoxP3_minus", "full_Ctrl2_LI_FoxP3_minus", "full_Ctrl3_LI_FoxP3_minus")] = mids_counts_full[c("Ctrl1_LI_FoxP3-", "Ctrl2_LI_FoxP3-", "Ctrl3_LI_FoxP3-")]
mids_counts[rownames(mids_counts_va8),  c("va8_Ctrl1_LI_FoxP3_minus", "va8_Ctrl2_LI_FoxP3_minus", "va8_Ctrl3_LI_FoxP3_minus")]    = mids_counts_va8[c("Ctrl1_LI_FoxP3-", "Ctrl2_LI_FoxP3-", "Ctrl3_LI_FoxP3-")]
mids_counts[rownames(mids_counts_amp),  c("amp_Ctrl1_LI_FoxP3_minus", "amp_Ctrl2_LI_FoxP3_minus", "amp_Ctrl3_LI_FoxP3_minus")]    = mids_counts_amp[c("Ctrl1_LI_FoxP3-", "Ctrl2_LI_FoxP3-", "Ctrl3_LI_FoxP3-")]

selected_rows = rownames(mids_counts[which(rowSums(mids_counts) != 0),])

mids_counts = mids_counts[selected_rows, ]

make_full_heatmap(mids_counts[1:9],   method="horn",      dirname="./", prefix="inputs",         binary=F)
make_full_heatmap(mids_counts[1:9],   method="intersect", dirname="./", prefix="inputs",         range=c(0, NA))

make_full_heatmap(mids_counts[10:18], method="horn",      dirname="./", prefix="ctrl_li_foxp3-", binary=F)
make_full_heatmap(mids_counts[10:18], method="intersect", dirname="./", prefix="ctrl_li_foxp3-", range=c(0, NA))

# top perc 85%

make_full_heatmap(mids_counts[1:9],   method="horn",      dirname="./", prefix="inputs_top85",         topPerc=0.85, binary=F)
make_full_heatmap(mids_counts[1:9],   method="intersect", dirname="./", prefix="inputs_top85",         topPerc=0.85, range=c(0, NA))
make_abundance_matrix(mids_counts[1:9], topPerc=0.85, prefix="inputs_top85")

make_full_heatmap(mids_counts[10:18], method="horn",      dirname="./", prefix="ctrl_li_foxp3-_top85", topPerc=0.85, binary=F)
make_full_heatmap(mids_counts[10:18], method="intersect", dirname="./", prefix="ctrl_li_foxp3-_top85", topPerc=0.85, range=c(0, NA))
make_abundance_matrix(mids_counts[10:18], topPerc=0.85, prefix="inputs_top85")
