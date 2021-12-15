library(rprojroot)

library("optparse")

option_list = list(
  make_option(c("-d", "--downsample"), type="character", default=NULL,
              help="whether to downsample to a fixed number of cells. Can either be a number or just True"),
  make_option(c("-g", "--group-by"), type="character", default=NULL,
              help="if set, in addition to the normal midX_clones.csv it will also pool together samples according to the specified variable (eg, mouse)")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)

root <- is_vcs_root$make_fix_file()

source(root("code/util.R"))

mid_labels <- read.csv("mid_labels.csv", row.names=1)
print(mid_labels)

if (is.null(opt$downsample)) {
  downsample = F
} else if (opt$downsample == "True" || opt$downsample == "T") {
  downsample = T
} else {
  downsample = as.integer(opt$downsample)
}

clones2groups(overwrite=T, savefasta=T, dirname="fasta", downsample=downsample)

if (!is.null(opt$'group-by')) {
  clones2groups(overwrite=T, savefasta=T, dirname="fasta_groups", join_by=opt$'group-by', collapse=T, downsample=downsample)
}
