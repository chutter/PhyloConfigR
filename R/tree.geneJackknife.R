#' @title analysis.geneJackknife
#'
#' @description Runs a gene jackknife on a set of individual locus alignments.
#'   Each replicate draws loci at random without replacement until a fixed amount
#'   of data is reached (a number of base pairs or a number of genes),
#'   concatenates them with concatenateAlignments(), and estimates a maximum
#'   likelihood tree with analysis.concatenationTree(). The resampling is the
#'   replication, so no per-replicate bootstrap is run. Optionally builds a
#'   majority-rule consensus from the replicate trees. Replicates can be run in
#'   any order or split across array jobs and the result is the same, because
#'   each replicate is seeded from its own index.
#'
#' @param alignment.path path to the folder of individual phylip-format
#'   alignments
#'
#' @param output.directory path to the output directory; replicate trees are
#'   written to output.directory/trees/rep_XXXX.tre
#'
#' @param selection.method how to size each replicate: "basepairs" draws loci
#'   until the concatenated length reaches jackknife.size, "genes" draws exactly
#'   jackknife.size loci (default: "basepairs")
#'
#' @param jackknife.size target size of each replicate: base pairs when
#'   selection.method = "basepairs", or number of loci when
#'   selection.method = "genes"
#'
#' @param replicates total number of jackknife replicates; also sets the width
#'   of the rep_XXXX labels (default: 1000)
#'
#' @param replicate.subset optional vector of replicate indices to run in this
#'   call; NULL runs 1:replicates. Use this to split the run across array jobs,
#'   one index per job (default: NULL)
#'
#' @param locus.lengths optional shortcut to avoid re-reading every alignment:
#'   a path to a two- or three-column file (locus,length[,taxa]) or a data.frame
#'   with those columns. When NULL the phylip header of every alignment is read
#'   to get its length and taxon count (default: NULL)
#'
#' @param locus.completeness minimum percent taxon completeness a locus must
#'   have to be eligible, measured against the most complete locus. Needs a
#'   taxon count per locus, so it is only applied when the alignments are scanned
#'   or locus.lengths carries a taxa column; set to 0 to disable (default: 0)
#'
#' @param min.locus.size loci at or below this length in base pairs are dropped
#'   before sampling (default: 100)
#'
#' @param outgroup outgroup taxon (or taxa) used to root the trees before the
#'   consensus; NULL leaves the consensus unrooted (default: NULL)
#'
#' @param partition.scheme partition handling passed to
#'   analysis.concatenationTree(): "none", "merge", or "file". A jackknife
#'   normally uses "none" with a fixed model, since model selection on hundreds
#'   of matrices costs far more than it changes the topology (default: "none")
#'
#' @param model model passed to IQ-TREE when partition.scheme = "none"
#'   (default: "GTR")
#'
#' @param msub.type substitution model category: "nuclear" or "mitochondrial"
#'   (default: "nuclear")
#'
#' @param codon.partition if TRUE adds -st CODON for codon-aware models
#'
#' @param rcluster percent of partitions used by rcluster when
#'   partition.scheme = "merge"
#'
#' @param threads CPU threads per replicate tree
#'
#' @param memory memory in GB per replicate tree
#'
#' @param iqtree.path path to an IQ-TREE executable or its directory; NULL uses
#'   iqtree2 or iqtree on the system PATH
#'
#' @param consensus if TRUE builds a majority-rule consensus from the replicate
#'   trees after the replicates in this call finish (default: TRUE)
#'
#' @param consensus.p proportion for the majority-rule consensus, passed to
#'   ape::consensus (default: 0.5)
#'
#' @param seed base random seed; replicate i is seeded with seed + i so each
#'   replicate is reproducible on its own (default: 0)
#'
#' @param keep.matrices if TRUE keeps each replicate's concatenated matrix;
#'   FALSE deletes it once the tree is built (default: FALSE)
#'
#' @param quiet if TRUE passes -quiet to IQ-TREE and prints less
#'
#' @param resume if TRUE skips a replicate whose tree already exists, so a
#'   partial run can be topped up (default: TRUE)
#'
#' @param overwrite if TRUE removes an existing output.directory before running
#'
#' @return invisibly returns the paths to the replicate treefiles written or
#'   found in this call
#'
#' @examples
#'
#' #Full run: 1000 replicates of 1,000,000 bp each, unpartitioned GTR+G.
#'
#' analysis.geneJackknife(alignment.path = "alignments/trimmed",
#'                        output.directory = "jackknife",
#'                        selection.method = "basepairs",
#'                        jackknife.size = 1000000,
#'                        replicates = 1000,
#'                        model = "GTR+G",
#'                        outgroup = "Kalophrynus_pleurostigma",
#'                        threads = 4)
#'
#' #One replicate only, no consensus, for a cluster array job:
#'
#' analysis.geneJackknife(alignment.path = "alignments/trimmed",
#'                        output.directory = "jackknife",
#'                        jackknife.size = 1000000,
#'                        replicates = 200,
#'                        replicate.subset = 37,
#'                        locus.lengths = "inputs/locus_lengths.csv",
#'                        model = "GTR+G",
#'                        consensus = FALSE,
#'                        threads = 4)
#'
#' @export

analysis.geneJackknife = function(alignment.path = NULL,
                                  output.directory = NULL,
                                  selection.method = c("basepairs", "genes"),
                                  jackknife.size = NULL,
                                  replicates = 1000,
                                  replicate.subset = NULL,
                                  locus.lengths = NULL,
                                  locus.completeness = 0,
                                  min.locus.size = 100,
                                  outgroup = NULL,
                                  partition.scheme = c("none", "merge", "file"),
                                  model = "GTR",
                                  msub.type = c("nuclear", "mitochondrial"),
                                  codon.partition = FALSE,
                                  rcluster = 100,
                                  threads = 1,
                                  memory = 1,
                                  iqtree.path = NULL,
                                  consensus = TRUE,
                                  consensus.p = 0.5,
                                  seed = 0,
                                  keep.matrices = FALSE,
                                  quiet = FALSE,
                                  resume = TRUE,
                                  overwrite = FALSE) {

  #match.arg picks the first choice when the argument is left at its default,
  #the same reason analysis.concatenationTree() uses it.
  selection.method = match.arg(selection.method)
  partition.scheme = match.arg(partition.scheme)
  msub.type = match.arg(msub.type)

  if (is.null(alignment.path) || dir.exists(alignment.path) == FALSE){ stop("A valid alignment.path is needed.") }
  if (is.null(output.directory)){ stop("An output.directory is needed.") }
  if (is.null(jackknife.size) || length(jackknife.size) != 1 || is.numeric(jackknife.size) == FALSE || jackknife.size <= 0){
    stop("jackknife.size must be a single positive number.")
  }
  if (length(replicates) != 1 || is.numeric(replicates) == FALSE || replicates < 1){ stop("replicates must be a positive number.") }
  if (resume == TRUE && overwrite == TRUE){
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  #Which replicates this call is responsible for. NULL means the whole run;
  #a subset lets one array task own a single index without changing the labels.
  if (is.null(replicate.subset)){
    replicate.subset = seq_len(as.integer(replicates))
  } else {
    replicate.subset = as.integer(replicate.subset)
    if (any(is.na(replicate.subset)) || any(replicate.subset < 1) || any(replicate.subset > replicates)){
      stop("replicate.subset must be between 1 and replicates.")
    }
  }

  #########################################
  # 1. Build the eligible locus set
  #########################################

  #Two ways to get each locus length. Reading the phylip header of every file is
  #self-contained but re-reads the folder on every call; a precomputed
  #locus.lengths table skips that, which matters when an array launches hundreds
  #of tasks against the same alignments.
  if (is.null(locus.lengths)){
    locus.names = list.files(alignment.path, full.names = FALSE, recursive = FALSE)
    if (length(locus.names) == 0){ stop("No alignments were found in ", alignment.path) }
    locus.data = data.table(locus = locus.names, length = as.numeric(NA), taxa = as.numeric(NA))
    for (i in 1:length(locus.names)){
      #The phylip header is "<ntax> <nchar>", so the first line gives both the
      #taxon count and the length without reading the alignment body.
      header = strsplit(trimws(readLines(file.path(alignment.path, locus.names[i]), n = 1)), "\\s+")[[1]]
      if (length(header) < 2){ next }
      locus.data$taxa[i] = suppressWarnings(as.numeric(header[1]))
      locus.data$length[i] = suppressWarnings(as.numeric(header[2]))
    }#end i loop
    locus.data = locus.data[is.na(locus.data$length) == FALSE, ]
  } else {
    if (is.character(locus.lengths) && length(locus.lengths) == 1){
      if (file.exists(locus.lengths) == FALSE){ stop("locus.lengths file not found: ", locus.lengths) }
      locus.data = fread(locus.lengths)
    } else { locus.data = as.data.table(locus.lengths) }
    if (ncol(locus.data) < 2){ stop("locus.lengths needs at least a locus column and a length column.") }
    #Take the first two columns as locus and length, and a "taxa" column if present.
    setnames(locus.data, 1:2, c("locus", "length"))
    if ("taxa" %in% colnames(locus.data) == FALSE){ locus.data$taxa = as.numeric(NA) }
    locus.data$length = as.numeric(locus.data$length)
    locus.data = locus.data[is.na(locus.data$length) == FALSE, ]
  }

  if (nrow(locus.data) == 0){ stop("No usable locus lengths were found.") }

  #########################################
  # A. Filter loci by length and completeness
  #########################################
  locus.data = locus.data[locus.data$length > min.locus.size, ]

  if (locus.completeness > 0){
    if (all(is.na(locus.data$taxa))){
      warning("locus.completeness > 0 but no taxon counts are available; completeness filtering skipped.")
    } else {
      max.taxa = max(locus.data$taxa, na.rm = TRUE)
      locus.data = locus.data[(locus.data$taxa / max.taxa) * 100 > locus.completeness, ]
    }
  }

  if (nrow(locus.data) == 0){ stop("No loci remain after filtering; loosen min.locus.size or locus.completeness.") }

  total.bp = sum(as.numeric(locus.data$length))
  n.loci = nrow(locus.data)

  #Sizing sanity checks, up front so a bad target fails before any tree is built.
  if (selection.method == "genes"){
    if (jackknife.size >= n.loci){ stop("jackknife.size (", jackknife.size, ") must be smaller than the ", n.loci, " eligible loci.") }
  } else {
    if (jackknife.size > total.bp){ stop("jackknife.size (", jackknife.size, " bp) exceeds the ", total.bp, " bp across ", n.loci, " eligible loci.") }
  }

  #########################################
  # 2. Output layout
  #########################################
  if (dir.exists(output.directory) == TRUE){
    if (overwrite == TRUE){ unlink(output.directory, recursive = TRUE) }
  }
  trees.dir = file.path(output.directory, "trees")
  picks.dir = file.path(output.directory, "picks")
  dir.create(trees.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(picks.dir, recursive = TRUE, showWarnings = FALSE)

  #Label width follows the total number of replicates so rep_0007 and rep_1000
  #sort together, independent of which subset this call runs.
  label.width = max(4, nchar(as.character(as.integer(replicates))))
  locus.list = locus.data$locus

  #########################################
  # 3. Run the replicates
  #########################################
  written = c()
  for (rep in replicate.subset){

    rep.tag = paste0("rep_", formatC(rep, width = label.width, flag = "0"))
    tree.out = file.path(trees.dir, paste0(rep.tag, ".tre"))

    if (resume == TRUE && file.exists(tree.out) == TRUE){
      if (quiet == FALSE){ print(paste0(rep.tag, " already has a tree, skipping.")) }
      written = c(written, tree.out)
      next
    }

    #Seeded from the replicate index so any replicate reproduces on its own and
    #the run does not depend on the order tasks are scheduled.
    set.seed(seed + rep)

    #Draw loci without replacement to the target size. cumsum on a shuffled
    #order is the same as adding one locus at a time until the target is met.
    draw.order = sample.int(n.loci)
    if (selection.method == "genes"){
      picked.loci = locus.list[draw.order[seq_len(as.integer(jackknife.size))]]
    } else {
      running.bp = cumsum(as.numeric(locus.data$length[draw.order]))
      n.take = which(running.bp >= jackknife.size)[1]
      picked.loci = locus.list[draw.order[seq_len(n.take)]]
    }
    picked.bp = sum(as.numeric(locus.data$length[locus.data$locus %in% picked.loci]))

    if (quiet == FALSE){
      cat(sprintf("%s: %d loci, %.0f bp (target %.0f)\n", rep.tag, length(picked.loci), picked.bp, jackknife.size))
    }
    writeLines(picked.loci, file.path(picks.dir, paste0(rep.tag, "_loci.txt")))

    #concatenateAlignments does not create its output folder, so make it first.
    #Samples absent from a locus are filled with Ns, which is what makes matrices
    #from different draws comparable.
    work.dir = file.path(output.directory, "work", rep.tag)
    dir.create(work.dir, recursive = TRUE, showWarnings = FALSE)

    concatenateAlignments(alignment.path = alignment.path,
                          alignment.names = picked.loci,
                          file.name = rep.tag,
                          output.dir = work.dir,
                          partition.format = "none")

    matrix.file = file.path(work.dir, paste0(rep.tag, ".phy"))
    if (file.exists(matrix.file) == FALSE || file.size(matrix.file) == 0){
      stop("concatenateAlignments produced no matrix at ", matrix.file)
    }

    #A sample present at none of the drawn loci is simply absent from the matrix,
    #so a replicate can carry fewer than the full set of taxa. That is expected
    #at small targets; recorded here so it shows in the run log.
    header = suppressWarnings(as.integer(strsplit(trimws(readLines(matrix.file, n = 1)), "\\s+")[[1]]))
    n.taxa = header[1]
    n.sites = header[2]
    if (quiet == FALSE){ cat(sprintf("%s: matrix %d taxa x %d sites\n", rep.tag, n.taxa, n.sites)) }

    #uf.bootstrap = 0: the jackknife is the replication, a per-replicate bootstrap
    #would only add runtime. seed + rep gives IQ-TREE its own reproducible seed.
    tree.file = analysis.concatenationTree(alignment.file = matrix.file,
                                           output.directory = file.path(work.dir, "tree"),
                                           output.name = rep.tag,
                                           partition.scheme = partition.scheme,
                                           model = model,
                                           msub.type = msub.type,
                                           codon.partition = codon.partition,
                                           rcluster = rcluster,
                                           uf.bootstrap = 0,
                                           threads = threads,
                                           memory = memory,
                                           iqtree.path = iqtree.path,
                                           seed = seed + rep,
                                           quiet = quiet,
                                           resume = resume)

    file.copy(tree.file, tree.out, overwrite = TRUE)

    write.table(data.frame(replicate = rep,
                           n_loci = length(picked.loci),
                           bp = picked.bp,
                           n_taxa = n.taxa,
                           n_sites = n.sites,
                           model = model),
                file = file.path(picks.dir, paste0(rep.tag, "_summary.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    if (keep.matrices == FALSE){ unlink(work.dir, recursive = TRUE) }

    written = c(written, tree.out)
    print(paste0(rep.tag, " done -> ", tree.out))

  }#end rep loop

  #########################################
  # 4. Majority-rule consensus
  #########################################
  if (consensus == TRUE){
    tree.files = list.files(trees.dir, pattern = "\\.tre$", full.names = TRUE)
    if (length(tree.files) < 2){
      warning("Fewer than 2 replicate trees are present; consensus skipped.")
    } else {
      #Loop rather than read.tree the whole folder at once, so each tree can be
      #rooted on the outgroup before the consensus is taken.
      rep.trees = vector("list", length(tree.files))
      for (i in 1:length(tree.files)){
        temp.tree = ape::read.tree(tree.files[i])
        if (is.null(outgroup) == FALSE && all(outgroup %in% temp.tree$tip.label)){
          temp.tree = ape::root(temp.tree, outgroup = outgroup, resolve.root = TRUE)
        }
        rep.trees[[i]] = temp.tree
      }#end i loop
      class(rep.trees) = "multiPhylo"

      cons.tree = ape::consensus(rep.trees, p = consensus.p, check.labels = TRUE)
      consensus.file = file.path(output.directory, "majority_consensus.tre")
      ape::write.tree(cons.tree, file = consensus.file)
      print(paste0("Majority-rule consensus (p = ", consensus.p, ") of ",
                   length(tree.files), " trees written to ", consensus.file))
    }
  }

  return(invisible(written))

}#end function
