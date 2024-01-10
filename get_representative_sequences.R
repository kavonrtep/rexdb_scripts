#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
library(optparse)

# get command line arguments
option_list <-  list(
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="fasta file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output files prefix"),
  make_option(c("-m", "--min_cluster_size"), type="integer", default=4,
              help="minimum cluster size"),
  make_option(c("-c", "--cutoff"), type="numeric", default=0.6,
              help="distance cutoff"),
  make_option(c("-n", "--n_threads"), type="integer", default=8,
              help="number of threads, default 8"),
  make_option(c("-t", "--tree"), type="character", default=NULL,
              help="tree file from mafft (optional)"),
  make_option(c("-A", "--auto_threshold"), type="logical", default=FALSE,
              help="use automatic threshold estimation", action = "store_true"
  )
)


opt_parser <-  OptionParser(option_list=option_list)
opt <-  parse_args(opt_parser)



suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(Biostrings)
  library(parallel)
})

mafft_tree <- function(fasta,cpu=8){
  tmpdir <- tempdir()
  dir.create(tmpdir, showWarnings = FALSE)
  # copy fasta to tmpdir
  ftmp <-paste0(tmpdir,"/", 'seq.fasta')
  file.copy(fasta, ftmp)
  cmd <- paste0("mafft --treeout   --thread ",cpu, " " , ftmp, " > /dev/null")
  print(cmd)
  system(cmd)
  tree <- read.tree(paste0(ftmp, ".tree"))
  unlink(tmpdir, recursive = TRUE)
  return(tree)
}

get_subtree <- function(node, phylo_tree) {
  desc_tips <- Descendants(phylo_tree, node, "tips")
  subtree <- keep.tip(phylo_tree, desc_tips[[1]])
  return(subtree)
}

most_central_node <- function(tree){
  d <-  cophenetic.phylo(tree)
  sd <- colSums(d)
  most_central_node <- which.min(sd)
  return(most_central_node)
}


mafft_alingment <- function(aaseq, cpu=8){
  tmpdir <- tempdir()
  dir.create(tmpdir, showWarnings = FALSE)
  fin <- file.path(tmpdir, "seq.fasta")
  fout <- file.path(tmpdir, "seq.aln")
  writeXStringSet(aaseq, file.path(tmpdir, "seq.fasta"))
  cmd <- paste0("mafft --thread ", cpu, " ", file.path(tmpdir, "seq.fasta"), " >" , fout)
  system(cmd)
  aln <- readAAStringSet(fout)
  unlink(tmpdir, recursive = TRUE)
  return(aln)
}


get_representative <- function(ss, ...){
  aln <- mafft_alingment(ss, ...)
  conseq <- consensusString(aln)
  conseq_char <- strsplit(conseq, "")[[1]]
  aln_char <- strsplit(as.character(aln), "")
  score <- sapply(aln_char, function(x) sum(x == conseq_char)/length(x))
  # score also end gaps
  # count beggining gaps
  end1_gap <- nchar(aln) - nchar(gsub("^-+","", aln))
  end2_gap <- nchar(aln) - nchar(gsub("-+$","", aln))
  end1_gap_conseq <- nchar(conseq) - nchar(gsub("^-+","", conseq))
  end2_gap_conseq <- nchar(conseq) - nchar(gsub("-+$","", conseq))
  end1_gap == end1_gap_conseq
  end2_gap == end2_gap_conseq
  return(names(ss[which.max(score)]))
}


fin <- opt$fasta
if (is.null(opt$tree)) {
  message("No tree provided, calculating tree with mafft")
  phylo_tree <- mafft_tree(fin, cpu = opt$n_threads)
} else {
  message("Using provided tree")
  phylo_tree <- read.tree(opt$tree)
}

s <- readAAStringSet(fin)

# Calculate the distance from the root to each node
node_heights <- node.depth.edgelength(phylo_tree)

# try to estimate suitable cutoff
plot.phylo(phylo_tree, show.node.label = FALSE, cex = 1.5, type= "t")
ltt_data <- ltt.plot.coords(phylo_tree)
save.image("tmp/tmp.RData")
# try first with default parameters, if it fails, try with tol=1e-6

spline_fit <- tryCatch({
  smooth.spline(ltt_data[,1], ltt_data[,2], spar = 1)
}, error = function(e) {
  smooth.spline(ltt_data[,1], ltt_data[,2], spar = 1, tol=1e-6)
})


# spline_fit <- smooth.spline(ltt_data[,1], ltt_data[,2], spar = 1)
# spline_fit <- smooth.spline(ltt_data[,1], ltt_data[,2], spar = 1, tol=1e-6)


second_derivative <- predict(spline_fit, deriv = 2)
inflection_points <- which(diff(sign(second_derivative$y)) != 0)
posible_thresholds <- max(node_heights) + spline_fit$x[inflection_points]

autothreshold <- posible_thresholds[which.min(abs(posible_thresholds - max(node_heights) * 0.5))]

if (opt$auto_threshold) {
  threshold <- autothreshold
} else {
  threshold <- max(node_heights) * opt$cutoff
}
#threshold <- autothreshold

#Identify the internal nodes that are closer to the root than the threshold
internal_nodes <- which(node_heights < threshold & seq_along(node_heights) > Ntip(phylo_tree))

subtrees <- mclapply(internal_nodes, get_subtree, phylo_tree = phylo_tree, mc.cores = opt$n_threads)


# sort trees by size
subtrees <- subtrees[order(sapply(subtrees, Ntip))]

all_tip <- phylo_tree$tip.label
final_trees <- list()
j <- 0
for (i in seq_along(subtrees)) {
  if (all(subtrees[[i]]$tip.label %in% all_tip)) {
    j <- j + 1
    all_tip <- setdiff(all_tip, subtrees[[i]]$tip.label)
    final_trees[[j]] <- subtrees[[i]]
  }
}
cluster_size <- sapply(final_trees, Ntip)

min_cluster_size <- opt$min_cluster_size
# get subset of sequences
s_part <- list()
repr <- character()
repr_table_list <- list()
j <- 0
for (i in seq_along(final_trees)) {
  print(i)
  if (cluster_size[i] < min_cluster_size) {
    next
  }
  j <- j + 1
  idx <- as.numeric(gsub("_.+","", final_trees[[i]]$tip.label))
  s_part[[j]] <- s[idx]
  repr[j] <- get_representative(s_part[[j]], cpu = opt$n_threads)
  repr_table_list[[j]] <- data.frame(representative=repr[j], ids=names(s_part[[j]]), length= nchar(s_part[[j]]))
}


s_repr <- s[names(s) %in% repr]

repr_table <- do.call(rbind, repr_table_list)

df_cluster_size <- data.frame(cluster_size = cluster_size[cluster_size>=min_cluster_size], representative = repr)



# export results
write.table(repr_table, paste0(opt$output, "_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_cluster_size, paste0(opt$output, "_cluster_size.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeXStringSet(s_repr, paste0(opt$output, "_representative.fasta"))

tipcolor <- ifelse(as.numeric(gsub("_.+", "", phylo_tree$tip.label)) %in% which(names(s) %in% repr), 2,"#00000000")
pngout <- paste0(opt$output, "_tree.png")
png(pngout, width = 4000, height = 3000)
plot(phylo_tree, show.node.label = TRUE, tip.color = tipcolor, cex = 1.5, type= "t")
abline(v = threshold, col = "red", lwd = 2)
abline(v = autothreshold, col = "blue")
xmax <- max(node_heights)
# recalculate scale , xmax is at 1
at <- seq(0, xmax, length.out = 11)
lab <-  seq(0, 1, length.out = 11)
axis(1, at = at, labels = lab)
dev.off()

# export tree
write.tree(phylo_tree, paste0(opt$output, "_tree.nwk"))

save.image(paste0(opt$output, "_data.RData"))