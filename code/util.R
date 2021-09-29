
set.seed(42)

library(methods)
library(divo)
library(vegan)
library(gplots)
library(venneuler)
library(plyr)
library(dplyr) # always after plyr
library(randomcoloR)
library(fmsb)
library(circlize)
library(RColorBrewer)
library(Rtsne)
library(matrixStats)
library(reshape2)
library(tidyr)
library(ggplot2)

library(devtools)
install_github("fabio-t/mixed.utils")
library(mixed.utils)

heatmaps <- function(regexp, prefix = NULL, invert = F, cex = 0.9, tsne = F) {
  if (is.null(regexp)) {
    ## All
    make_full_heatmap(mids_counts, "horn", mid_names, prefix = "full", binary = F, types = "square", cexRow = cex, cexCol = cex, tsne = tsne)
    make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full", types = "square", vlim = c(0, NA), cexRow = cex, cexCol = cex)
    make_full_heatmap(mids_counts, "intersect", mid_names, prefix = "full_top0.85", topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = cex, cexCol = cex)
  } else {
    # prefix is expected to be set, now
    columns_i <- grep(regexp, mid_names, perl = T, invert = invert)
    columns_n <- mid_names[columns_i]
    data <- mids_counts[columns_i]
    make_full_heatmap(data, "horn", columns_n, prefix = prefix, binary = F, types = "square", cexRow = cex, cexCol = cex, tsne = tsne)
    make_full_heatmap(data, "intersect", columns_n, prefix = prefix, types = "square", vlim = c(0, NA), cexRow = cex, cexCol = cex)
    make_full_heatmap(data, "intersect", columns_n, prefix = paste0(prefix, "_top0.85"), topPerc = 0.85, types = "square", vlim = c(0, NA), cexRow = cex, cexCol = cex)
  }
}

radarcharts <- function(regexp, dirname, save = F, invert = F, colours = NULL) {
  print(dirname)

  column_indices <- grep(regexp, mid_names, perl = T, invert = invert)

  column_names <- colnames(mids_counts)[column_indices]
  column_names_2 <- mid_names[column_indices]

  print(column_names_2)

  data <- mids_counts[column_names]

  dirname <- paste0("./radarcharts/one-vs-all/topPerc", tp, "/", dirname, "/")

  l <- save_multiple_radarcharts(data, column_names, labels = column_names_2,
                                 topPerc = tp, percentage = T, use_log = F,
                                 dirname = dirname, use_colours = colours,
                                 show_axis = F, show_legend = T)

  if (save && (length(l$common_clones) > 0)) {
    # common clones
    d1 <- sweep(data, 2, colSums(data), "/")[l$common_clones, , drop=F]
    d2 <- data[l$common_clones, , drop=F]
    d3 <- as.data.frame(colRanks(as.matrix(d2), preserveShape = T, ties.method = "average"))
    colnames(d1) <- column_names_2
    colnames(d2) <- column_names_2
    colnames(d3) <- column_names_2
    rownames(d3) <- l$common_clones
    write.csv(d1, file = paste0(dirname, "common_clones_perc.csv"))
    write.csv(d2, file = paste0(dirname, "common_clones_counts.csv"))
    write.csv(d3, file = paste0(dirname, "common_clones_ranks.csv"))
  }

  if (save && (length(l$total_clones) > 0)) {
    # common clones
    d1 <- sweep(data, 2, colSums(data), "/")[l$total_clones, , drop=F]
    d2 <- data[l$total_clones, , drop=F]
    d3 <- as.data.frame(colRanks(as.matrix(d2), preserveShape = T, ties.method = "average"))
    colnames(d1) <- column_names_2
    colnames(d2) <- column_names_2
    colnames(d3) <- column_names_2
    rownames(d3) <- l$total_clones
    write.csv(d1, file = paste0(dirname, "total_clones_perc.csv"))
    write.csv(d2, file = paste0(dirname, "total_clones_counts.csv"))
    write.csv(d3, file = paste0(dirname, "total_clones_ranks.csv"))
  }
}

# from https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
write.fasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

consensus_nt <- function(x){consensusString(DNAStringSet(x))}
consensus_aa <- function(x){consensusString(AAStringSet(x))}
transl_nt2aa <- function(x){as.character(translate(DNAString(x), if.fuzzy.codon="solve"))}

# by default, two clones are split if they have < 85% similarity
cloneclust <- function(x, h=0.15) {
  if(length(x) == 1){return(1)}

  aa_seqs <- as.matrix(as.AAbin(AAStringSet(x)))
  seq_dist <- dist.aa(aa_seqs, scaled=T)
  cl <- hclust(seq_dist, method="average")
  cl$height = round(cl$height, digits=8)
  # plot(cl)
  # rect.hclust(cl, h=h, border="red")
  cutree(cl, h=h)
}

immload <- function(which="full") {
  library(immunarch)

  if (which == "full") {
    immdata <- repLoad(paste0("full/", rownames(mid_labels), ".csv"), .coding=F)
  } else {
    immdata <- repLoad(paste0(tolower(rownames(mid_labels)), "_clones.csv"), .coding=F)
    names(immdata$data) <- toupper(gsub("_clones", "", names(immdata$data)))
  }

  meta <- cbind(Sample=rownames(mid_labels), mid_labels)
  immdata$meta <- as_tibble(meta)

  return(immdata)
}

clones2groups <- function(immdata = NULL, overwrite = F, savefasta = F) {
  if (is.null(immdata)) {
    immdata <- immload()
  }

  library(dplyr)
  library(stringr)
  library(Biostrings)
  library(ape)

  for (name in names(immdata$data)) {
    d <- immdata$data[[name]]

    d2 <- d %>%
          mutate(len=str_length(CDR3.nt)) %>%
          group_by(V.name, J.name, len) %>%
          mutate(id=str_c(cur_group_id(), cloneclust(CDR3.aa, 0.15), sep=".")) %>%
          group_by(id, .add=T)
    print(d2)

    if (savefasta) {
      ## (careful with cumsum, because of how clones are grouped)

      ## keep only the top groups
      d3 <- d2 %>% mutate(n=n()) %>% ungroup() %>% mutate(gn = dense_rank(dplyr::desc(n))) %>% filter(gn %in% 1:20)
      ## keeps only subclones whose abundance is at least 5% of the whole
      ## lineage (it keeps all lineages but prunes them)
      # d3 <- d2 %>% filter(Clones/sum(Clones) >= 0.05)
      ## keeps only subclones that sumup to 50% total abundance
      ## (drops most singleton subclones and many lineages, but it may be too harsh,
      ## and still let many singleton LINEAGES in the final set)
      # d3 <- d2 %>% ungroup() %>% arrange(desc(Clones)) %>% filter(cumsum(Proportion) <= 0.75)
      ## keeps only lineages that sumup to 50% total abundance
      ## (will potentially keep enormous trees full of singletons..
      ## but we can farther filter out later)
      ### TBD

      print(d3)
      for (clone_id in unique(d3$id)) {
        print(clone_id)
        d4 <- d3 %>%
              ungroup() %>%
              filter(id == clone_id) %>%
              arrange(desc(Clones)) %>%
              slice_head(n=50) # keep from having too-large trees
        if (nrow(d4) < 2) next # drop empty clones, but also singletons
        print(d4)
        seqid <- paste(d4$CDR3.nt, d4$CDR3.aa, d4$Clones, sep="|")
        seq <- substr(d4$Sequence, 1, ceiling(d4$V.end/3)*3) # need a multiple of 3
        # seq <- substr(d4$Sequence, 1, ceiling(d4$J.start/3)*3-3) # need a multiple of 3
        fasta_df <- data.frame(name=seqid, seq=seq)
        dirname = make_path(paste0("fasta/", name))
        print(dirname)
        filename = paste0(dirname, d4$V.name[1], "_", d4$J.name[1], "_", d4$len[1], "_", d4$id[1], ".fasta")
        print(filename)
        write.fasta(fasta_df, filename)
      }
    }

    # FIXME D.name, as well as V.end/J.start/D.start/D.end, will not necessarily
    # match when V, J, CDR3length and CDR3@85% end up grouped together
    d3 <- d2 %>%
          summarise(cloneId=str_c(unique(id)), n=n(),
                    Clones=sum(Clones), Proportion=sum(Proportion),
                    CDR3.nt=consensus_nt(CDR3.nt),
                    # CDR3.aa=consensus_aa(CDR3.aa),
                    CDR3.aa=transl_nt2aa(CDR3.nt),

                    # some extra stuff to caracterise the lineages
                    # expsha=exp(diversity(Clones, index="shannon")),
                    diversity=diversity(Clones, index="invsim"),
                    # shan_even=diversity(Clones, index="shannon")/log(n()),
                    # shan_even2=exp(diversity(Clones, index="shannon"))/n(),
                    evenness=diversity(Clones, index="invsim")/n(),

                    # stuff for immunarch
                    D.name=str_c(unique(D.name), collapse=","),
                    V.end=str_c(unique(V.end), collapse=","),
                    D.start=str_c(unique(D.start), collapse=","),
                    D.end=str_c(unique(D.end), collapse=","),
                    J.start=str_c(unique(J.start), collapse=","),
                    VJ.ins=str_c(unique(VJ.ins), collapse=","),
                    VD.ins=str_c(unique(VD.ins), collapse=","),
                    DJ.ins=str_c(unique(DJ.ins), collapse=","),
                    Sequence=str_c(unique(Sequence), collapse=",")
                    ) %>%
          # slice_head() %>%
          # ungroup() %>%
          arrange(desc(Clones))

    immdata$data[[name]] <- d3
  }

  if (overwrite) {
    for (name in names(immdata$data)) {
      d <- immdata$data[[name]]
      d <- d %>%
           # ungroup() %>%
           select(cloneCount = Clones,
                  nSeqCDR3 = CDR3.nt, aaSeqCDR3 = CDR3.aa,
                  bestVHit = V.name, bestJHit = J.name,
                  bestDHit = D.name,
                  CDR3Length = len,
                  cloneSize = n,
                  cloneDiversity = diversity,
                  cloneEvenness = evenness,
                  cloneId,
                  cloneFraction = Proportion) %>%
           as.data.frame()

      d$bestDHit = NA # removing D because of multiple hits

      write.table(d, file=paste0(tolower(name), "_clones.csv"), quote=F, row.names=F, sep="\t")
      # write.table(d, file=paste0("full/", name, ".csv"), quote=F, row.names=F, sep="\t")
    }
  }

  return(immdata)
}

overview <- function(dirname = "overview", immdata = NULL, exclude = NULL) {
  if (is.null(immdata)) {
    immdata <- immload()
  }

  # downsamples to smallest clonality (NOT clonotypes)
  immdata2 <- immdata
  immdata2$data <- repSample(immdata$data, .method="downsample")
  n_clntps <- min(sapply(immdata2$data, nrow))
  n_clones <- min(sapply(immdata2$data, function(x){sum(x$Clones)}))

  print(n_clntps)
  print(n_clones)

  if (n_clones < 10) {
    raref_step <- 1
  } else {
    raref_step <- NA
  }

  dirname <- make_path(dirname)
  print(dirname)

 # + theme_minimal(base_size=12) + theme(axis.text.x = element_text(angle = 90))

  t1 <-  repExplore(immdata$data,  .method = "volume")
  t1b <- repExplore(immdata2$data, .method = "volume")
  t2 <-  repExplore(immdata$data,  .method = "count")
  t2b <- repExplore(immdata2$data, .method = "count")
  p1 <-  vis(t1)
  p1b <- vis(t1b)
  p2 <-  vis(t2)
  p2b <- vis(t2b)
  p1 + p2 + p1b + p2b
  ggsave(filename = paste0(dirname, "basic.pdf"), width = 16, height = (9/16) * 16)

  # hill curves
  hill_1 <- repDiversity(immdata$data,  "hill", .col="nt+v+j")
  hill_2 <- repDiversity(immdata2$data, "hill", .col="nt+v+j")
  p1 <- vis(hill_1, .meta = immdata$meta)
  p2 <- vis(hill_2, .meta = immdata2$meta)
  p1 + p2
  ggsave(filename = paste0(dirname, "hill.pdf"), width = 16, height = (9/16) * 16)

  # # true diversity
  # div_1 <- repDiversity(immdata$data,  "div", .col="nt+v+j")
  # div_2 <- repDiversity(immdata2$data, "div", .col="nt+v+j")
  # p1 <- vis(div_1, .meta = immdata$meta)
  # p2 <- vis(div_2, .meta = immdata2$meta)
  # p1 + p2
  # ggsave(filename = paste0(dirname, "truediv.pdf"), width = 16, height = (9/16) * 16)

  # # gini-simpson index
  # ginisimp_1 <- repDiversity(immdata$data,  "gini.simp", .col="nt+v+j")
  # ginisimp_2 <- repDiversity(immdata2$data, "gini.simp", .col="nt+v+j")
  # p1 <- vis(ginisimp_1, .meta = immdata$meta)
  # p2 <- vis(ginisimp_2, .meta = immdata2$meta)
  # p1 + p2
  # ggsave(filename = paste0(dirname, "ginisimpson.pdf"), width = 16, height = (9/16) * 16)

  # inverse-simpson index
  invsimp_1 <- repDiversity(immdata$data,  "inv.simp", .col="nt+v+j")
  invsimp_2 <- repDiversity(immdata2$data, "inv.simp", .col="nt+v+j")
  p1 <- vis(invsimp_1, .meta = immdata$meta)
  p2 <- vis(invsimp_2, .meta = immdata2$meta)
  p1 + p2
  ggsave(filename = paste0(dirname, "invsimpson.pdf"), width = 16, height = (9/16) * 16)

  # number of clonotypes summing up to 50% abundance
  d50_1 <- repDiversity(immdata$data,  "d50", .col="nt+v+j")
  d50_2 <- repDiversity(immdata2$data, "d50", .col="nt+v+j")
  p1 <- vis(d50_1, .meta = immdata$meta)
  p2 <- vis(d50_2, .meta = immdata2$meta)
  p1 + p2
  ggsave(filename = paste0(dirname, "d50.pdf"), width = 16, height = (9/16) * 16)

  # rarefaction
  raref_1 <- repDiversity(immdata$data,  "raref", .verbose = F, .norm=F, .col = "nt+v+j", .step = raref_step)
  raref_2 <- repDiversity(immdata2$data, "raref", .verbose = F, .norm=F, .col = "nt+v+j", .step = raref_step)
  p1 <-  vis(raref_1, .meta = immdata$meta)
  p1b <- vis(raref_1, .meta = immdata$meta, .log=T)
  p2 <-  vis(raref_2, .meta = immdata2$meta)
  p2b <- vis(raref_2, .meta = immdata2$meta, .log=T)
  p1 + p2
  ggsave(filename = paste0(dirname, "raref.pdf"), width = 16, height = (9/16) * 16)
  p1b + p2b
  ggsave(filename = paste0(dirname, "raref-log.pdf"), width = 16, height = (9/16) * 16)

  # rarefaction (normalised)
  raref_1 <- repDiversity(immdata$data,  "raref", .verbose = F, .norm=T, .col = "nt+v+j", .step = raref_step)
  raref_2 <- repDiversity(immdata2$data, "raref", .verbose = F, .norm=T, .col = "nt+v+j", .step = raref_step)
  p1 <-  vis(raref_1, .meta = immdata$meta)
  p1b <- vis(raref_1, .meta = immdata$meta, .log=T)
  p2 <-  vis(raref_2, .meta = immdata2$meta)
  p2b <- vis(raref_2, .meta = immdata2$meta, .log=T)
  p1 + p2
  ggsave(filename = paste0(dirname, "raref-norm.pdf"), width = 16, height = (9/16) * 16)
  p1b + p2b
  ggsave(filename = paste0(dirname, "raref-log-norm.pdf"), width = 16, height = (9/16) * 16)
}

track_clones <- function(regexp, dirname, invert = F, n = 15, immdata = NULL) {
  if (is.null(immdata)) {
    immdata <- immload()
  }

  if (is.null(regexp)) {
    # one-vs-all
    dirname <- paste0("./tracking/all/")
    dirname <- make_path(dirname)
    print(dirname)

    mid_list <- rownames(mid_labels)
  } else {
    dirname <- paste0("./tracking/", dirname, "/")
    dirname <- make_path(dirname)
    print(dirname)

    column_indices <- grep(regexp, mid_names, perl = T, invert = invert)

    mid_list <- colnames(mids_counts)[column_indices]
    label_list <- mid_names[column_indices]

    print(label_list)

    immdata$data <- immdata$data[column_indices]
  }

  for (mid in mid_list) {
    tc1 <- trackClonotypes(immdata$data, list(mid, n), .col = "v+nt+j")
    vis(tc1) + theme_minimal(base_size=12) + theme(axis.text.x = element_text(angle = 90))
    ggsave(filename = paste0(dirname, mid, ".pdf"), width = 16, height = (9/16) * 16)
  }
}

make_diversity <- function(mids_counts, mid_names) {
  divers_df <- as.data.frame(mid_names)
  rownames(divers_df) <- colnames(mids_counts)
  colnames(divers_df) <- "Sample"
  t_df <- t(mids_counts)
  divers_df$NumberOfCells <- rowSums(t_df)
  divers_df$NumberOfClones <- specnumber(t_df)
  divers_df$Shannon <- diversity(t_df, index = "shannon")
  divers_df$ShannonEvenness <- divers_df$Shannon / log(divers_df$NumberOfClones)
  divers_df$InvSimpson <- diversity(t_df, index = "invsimpson")
  divers_df$InvSimpsonEvenness <- divers_df$InvSimpson / divers_df$NumberOfClones
  write.csv(divers_df, file = "diversity.csv")
}

make_full_heatmap <- function(data, method, labels = NULL, binary = F, dirname = "heatmaps",
                              prefix = "full", clustering = F, topN = NULL, topPerc = NULL, cexRow = 0.5,
                              cexCol = 0.5, types = c("square", "nonsquare"), vlim = NULL, scale = "none", tsne = F) {
  # the returned "data" must overwrite the parameter "data",
  # because we modified some of the cells (see the function body)
  data <- get_top(data, topN = topN, topPerc = topPerc)

  if (!is.null(labels)) {
    colnames(data) <- labels
  }

  print(head(data))
  print(dim(data))

  if ("square" %in% types) {
    k <- ifelse(clustering, 4, 0)
    l <- make_heatmap(data,
      dirname = dirname, prefix = paste0(prefix, "_", method, "_heatmap-square"),
      threshold = 0, k = k, dist_method = method, abs = F, method = "ward.D2", scale = scale,
      square = T, output = cairo_pdf, cexRow = cexRow, cexCol = cexRow, vlim = vlim
    )

    if (tsne) {
      set.seed(42)
      cairo_pdf(paste0(dirname, "/", prefix, "_", method, "_tsne-square.pdf"))
      dist <- vegdist(t(data), method = "horn", diag = T, upper = T, binary = F)
      res <- Rtsne(dist, perplexity = ceiling(ncol(data) / 6), theta = 0.0, max_iter = 5000, eta = 1, is_distance = T)
      plot(res$Y, cex = 0.1, xlim = c(min(res$Y[, 1]) - 0.1, max(res$Y[, 1]) + 15), main = "Tsne (with MHI distance)")
      text(res$Y, colnames(data), cex = 0.5, pos = 4)
      dev.off()
    }
  }

  if ("nonsquare" %in% types) {
    k <- ifelse(clustering, 4, 0)
    k2 <- ifelse(clustering, 4, 0)
    l <- make_heatmap(data,
      dirname = dirname, prefix = paste0(prefix, "_", method, "_heatmap"),
      threshold = 0, k = k, k2 = k2, dist_method = method, abs = F, method = "ward.D2", scale = scale,
      square = F, output = cairo_pdf, cexRow = cexRow, cexCol = cexRow, vlim = vlim
    )
  }

  write.csv(l$data, file = paste0(dirname, "/", prefix, "_", method, "_matrix.csv"))
}

make_radarchart <- function(data, column_names, labels = NULL, topN = NULL, topPerc = NULL, minmax = "global", percentage = T, use_clones = NULL, show_axis = T,
                            use_union = T, sort_by_sum = T, legend_x = "topright", legend_inset = c(0, 0), legend_counts = F, add_clones = NULL, noplot = F, use_log = F,
                            colours_border = c(rgb(0.2, 0.5, 0.6, 0.8), rgb(1.0, 0.2, 0, 0.8), rgb(0.9, 0.7, 0.1, 0.8), rgb(0.2, 0.6, 0.3, 0.8)),
                            colours_in = c(rgb(0.2, 0.5, 0.6, 0.4), rgb(1.0, 0.2, 0, 0.4), rgb(0.9, 0.7, 0.1, 0.4), rgb(0.2, 0.6, 0.3, 0.4)), vlabels = NA) {
  if (is.null(labels)) {
    labels <- column_names
  }

  if (is.null(topN)) {
    topN <- Inf
  }

  data <- data[column_names]

  # backup for later
  orig_data <- data

  if (percentage) {
    # takes the sum of all reads for each organ,
    # then divides the single clonal abudance by the total (again by organ).
    # the output presents relative abundance values (hence, each organ sums up to 1)
    data_sums <- colSums(data)
    data <- sweep(data, 2, data_sums, "/")
  }

  common_rows <- c()

  if (is.null(use_clones)) {
    # since clones aren't provided, we use either the union or the intersection
    # of all our organs' clones
    for (column in column_names)
    {
      df <- subset(data[column], data[column] > 0)
      df <- df[order(df[, , drop = F], decreasing = T), , drop = F]

      if (is.null(topPerc)) {
        df_rows <- rownames(df[1:min(nrow(df), topN), , drop = F])
      }
      else {
        print(sum(df))
        abundance <- sum(df) * topPerc
        print(abundance)
        df_rows <- rownames(df)[cumsum(df[df_rows, ]) <= abundance]
        print(df[df_rows, ])
      }

      if (use_union) {
        common_rows <- union(common_rows, df_rows)
      }
      else {
        common_rows <- intersect(common_rows, df_rows)
      }
    }
  }
  else {
    # clones are provided, so we use them
    # note: topN is ignored

    common_rows <- use_clones
  }

  print(paste0("#clones: ", length(common_rows)))

  if (!is.null(add_clones)) {
    # trick: putting add_clones first compels the ordering to be based on add_clones
    common_rows <- union(add_clones, common_rows)
  }

  # order selected clones by total abundance
  if (sort_by_sum) {
    # for each of the selected clones, take the sum across organs
    clone_abund <- rowSums(data[common_rows, ])
    # sort the clones by their row-sum, bigger-to-smaller
    clone_abund <- sort(clone_abund, decreasing = T)

    # trick: intersect fixes the order using the first argument
    common_rows <- intersect(names(clone_abund), common_rows)
  }

  data <- as.data.frame(t(data[common_rows, ]))
  colnames(data) <- 1:ncol(data)

  # if counts are needed in the legend, we must go back
  # to the original data (ie, before getting percentanges)
  orig_data <- as.data.frame(t(orig_data[common_rows, ]))
  colnames(orig_data) <- 1:ncol(orig_data)
  organ_abund_counts <- rowSums(orig_data, na.rm=T)
  organ_abund_perc <- rowSums(data, na.rm=T)

  print("organ abundance:")
  print(organ_abund_counts)
  print(organ_abund_perc)

  if (minmax == "global") {
    # global minimum and maximum: this makes it easier to see
    # print(data)
    minv <- min(sapply(data, FUN=function(x){x=x[is.finite(x)]; min(x)}))
    print(minv)
    maxv <- max(sapply(data, FUN=function(x){x=x[is.finite(x)]; max(x)}))
    print(maxv)
    data <- rbind(maxv, minv, data)
  }
  else if (minmax == "by-column") {
    # by-organ max/min: might be hard to visualise
    data <- rbind(colMax(data, na.rm=T), colMin(data, na.rm=T), data)
  }
  else if (minmax == "by-row") {
    # by-organ max/min: might be hard to visualise
    data <- rbind(rowMax(data, na.rm=T), rowMin(data, na.rm=T), data)
  }
  else if (length(minmax) == 2) {
    # [min, max]

    data <- rbind(rep(minmax[2], ncol(data)), rep(minmax[1], ncol(data)), data)
  }
  else {
    stop("ERROR: minmax argument wrong")
  }

  print(data)

  if (!noplot) {
    if (use_log) {
      if (percentage) {
        data <- data * 100
      }

      # zero becomes one, so when log'ed it is still zero
      data <- log(data + 1)
    }

    print(data)

    if (minmax == "global" || length(minmax == 2)) {
      center_labels <- seq(data[2, 1], data[1, 1], length.out = 5)

      # un-log it
      if (use_log) {
        center_labels <- sprintf("%.5f", exp(center_labels) - 1)
      }
      else {
        center_labels <- sprintf("%.5f", center_labels)
      }

      # add the percent symbol
      if (percentage) {
        center_labels <- paste0(center_labels, "%")
      }
    }
    else {
      center_labels <- NULL
    }

    print(center_labels)

    title <- paste("#clones:", length(common_rows))

    if (show_axis) {
      axistype <- 1
    }
    else {
      axistype <- 0

      title <- paste(title, "| axis:", paste(center_labels, collapse = " "))
    }

    # radarchart(data, axistype=axistype, pcol=colours_border , pfcol=colours_in , plwd=4 , plty=1, cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, vlcex=0.8)
    radarchart(data,
      axistype = axistype, pcol = colours_border, pfcol = colours_in, plwd = 2, plty = 1, caxislabels = center_labels,
      cglcol = "grey", cglty = 1, axislabcol = "darkslategrey", cglwd = 0.8, vlcex = 0.8, maxmin = T, vlabels = vlabels, title = title
    )

    if (!is.null(legend_x)) {
      if (legend_counts || !percentage) {
        legend_values <- paste0(labels, " (", round(organ_abund_counts, digits = 3), ")")
      }
      else {
        legend_values <- paste0(labels, " (", round(organ_abund_perc * 100, digits = 3), "%)")
      }

      legend(legend_x, inset = legend_inset, legend = legend_values, bty = "n", pch = 20, col = colours_in, text.col = "grey", cex = 0.85, pt.cex = 2)
    }
  }

  return(list(
    clones = common_rows,
    data = data,
    organ_abund_counts = organ_abund_counts,
    organ_abund_perc = organ_abund_perc
  ))
}

## Add an alpha value to a colour
add.alpha <- function(col, alpha = 1) {
  if (missing(col)) {
    stop("Please provide a vector of colours.")
  }

  apply(sapply(col, col2rgb) / 255, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

save_multiple_radarcharts <- function(data, column_names, labels = NULL, topN = NULL, topPerc = NULL, dirname = "./radarcharts/one-vs-all", minmax = "global", percentage = T, use_log = F, use_colours = NULL, join_plots = F, show_legend = T, show_axis = T, legend_counts = F) {
  print("save_multiple_radarcharts:")
  print(column_names)

  dirname <- make_path(dirname)

  if (is.null(labels)) {
    labels <- column_names
  }

  if (is.null(topN)) {
    topN <- Inf
  }

  if (is.null(use_colours) || length(use_colours) < length(column_names)) {
    colours <- distinctColorPalette(length(column_names))
  }
  else {
    colours <- use_colours
  }

  colours_border <- add.alpha(colours, 0.8)
  colours_in <- add.alpha(colours, 0.08)

  if (join_plots) {
    svg(paste0(dirname, "all.svg"))

    par(xpd = TRUE, mar = c(0, 1, 0, 2), mfcol = c(ceiling(length(column_names) / 2), 2))
  }

  all_counts <- data.frame()

  ret <- list()

  # start with rownames to keep the ordering, then subset with "intersect"
  common_clones <- rownames(data)

  # we accumulate the selected clones of each radarchart, then at the end
  # we intersect with rownames to keep the ordering
  total_clones <- c()

  for (i in 1:length(column_names)) {
    column <- column_names[i]
    label <- labels[i]

    print(label)

    if (!join_plots) {
      svg(paste0(dirname, label, ".svg"))
    }

    df <- subset(data[column], data[column] > 0)
    df <- df[order(df[, , drop = F], decreasing = T), , drop = F]

    if (is.null(topPerc)) {
      df_rows <- rownames(df[1:min(nrow(df), topN), , drop = F])
    }
    else {
      abundance <- sum(df) * topPerc
      df_rows <- rownames(df)[cumsum(df) <= abundance]
    }

    if (length(df_rows) < 3) {
      print(paste0("ERROR: ", column, " has ", length(df_rows), " selected, but needs at least 3. Skipping."))
      next()
    }

    if (show_legend) {
      # each radarchart consists of the i-th organ topN clones sorted by abudance,
      # with the other organs following through
      l <- make_radarchart(data, column_names,
        labels = labels, use_clones = df_rows, sort_by_sum = F, minmax = minmax, percentage = percentage,
        use_log = use_log, colours_border = colours_border, colours_in = colours_in, show_axis = show_axis, legend_counts = legend_counts
      )
    }
    else {
      # each radarchart consists of the i-th organ topN clones sorted by abudance,
      # with the other organs following through
      l <- make_radarchart(data, column_names,
        labels = labels, use_clones = df_rows, sort_by_sum = F, minmax = minmax, percentage = percentage,
        use_log = use_log, colours_border = colours_border, colours_in = colours_in, legend_x = NULL, show_axis = show_axis
      )
    }

    organ_abund <- data.frame(
      Index = rep(label, length(labels)),
      Sample = labels,
      Counts = l$organ_abund_counts,
      Percentage = l$organ_abund_perc
    )
    write.csv(organ_abund, file = paste0(dirname, "/", label, ".csv"), row.names = F)

    all_counts <- rbind(all_counts, organ_abund)

    if (!join_plots) {
      dev.off()
    }

    common_clones <- intersect(common_clones, l$clones)
    total_clones <- union(total_clones, l$clones)
  }

  write.csv(all_counts, file = paste0(dirname, "all.csv"), row.names = F)

  if (join_plots) {
    dev.off()
  }

  print(paste0("common clones: ", length(common_clones)))
  print(paste0("total clones: ", length(total_clones)))

  ret$common_clones <- common_clones
  ret$total_clones <- intersect(rownames(data), total_clones)

  return(ret)
}
