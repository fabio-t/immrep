
library(DECIPHER)

args <- commandArgs(trailingOnly=T)

if (length(args) != 2) {
  stop("Input and output FASTA files must be provided as arguments.n", call.=FALSE)
}

seqs <- readDNAStringSet(args[1])
aligned <- AlignSeqs(seqs)
writeXStringSet(aligned, args[2])
