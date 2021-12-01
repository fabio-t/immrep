library(rprojroot)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)

clones2groups(overwrite=T, savefasta=T, dirname="fasta", downsample=1000)
clones2groups(overwrite=T, savefasta=T, dirname="fasta_groups", join_by="mouse", collapse=T, downsample=1000)
