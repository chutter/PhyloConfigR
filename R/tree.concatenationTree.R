#' @title analysis.concatenationTree
#'
#' @description Estimates a maximum likelihood tree from a concatenated
#'   alignment with IQ-TREE. Handles partitioned and unpartitioned runs, finds
#'   the IQ-TREE executable regardless of whether version 2 or 3 is installed,
#'   and makes UFBoot optional so the function can be used inside resampling
#'   procedures that do their own replication.
#'
#' @param alignment.file path to the concatenated alignment file in phylip format
#'
#' @param output.directory path to the parent output directory
#'
#' @param output.name name for the run subdirectory and output files
#'
#' @param partition.file path to a partition file (only used when
#'   partition.scheme = "file")
#'
#' @param partition.scheme how to handle partitions: "file" uses a provided
#'   partition file with ModelFinder, "merge" runs ModelFinder with partition
#'   merging, "none" runs unpartitioned (default: "file")
#'
#' @param model model passed to IQ-TREE when partition.scheme = "none"; ignored
#'   for the partitioned schemes, which use ModelFinder (default: "GTR")
#'
#' @param codon.partition if TRUE adds -st CODON flag for codon-aware model
#'
#' @param seq.type force the data type passed to IQ-TREE with -st (e.g. "DNA",
#'   "AA", "CODON"); NULL lets IQ-TREE auto-detect. Forcing the type avoids the
#'   "Unknown sequence type" failure some IQ-TREE builds hit when auto-detecting
#'   matrices with a lot of missing data (default: NULL)
#'
#' @param program reserved for future use; currently only "IQTREE" is supported
#'
#' @param msub.type substitution model category passed to IQ-TREE -msub flag;
#'   "nuclear" or "mitochondrial" (default: "nuclear")
#'
#' @param uf.bootstrap number of UFBoot replicates; use 0 to disable, otherwise
#'   1000 or more, which is IQ-TREE's own minimum (default: 1000)
#'
#' @param rcluster percentage of partitions used in the rcluster algorithm for
#'   partition merging
#'
#' @param threads number of CPU threads passed to IQ-TREE -nt flag
#'
#' @param memory memory in GB passed to IQ-TREE -mem flag
#'
#' @param iqtree.path path to an IQ-TREE executable or the directory containing
#'   it; use NULL if iqtree2 or iqtree is on the system PATH
#'
#' @param seed optional IQ-TREE random seed for reproducible runs
#'
#' @param quiet if TRUE passes -quiet to IQ-TREE
#'
#' @param resume if TRUE returns without running when a treefile already exists
#'
#' @param overwrite if TRUE removes the existing output directory before running
#'
#' @return invisibly returns the path to the treefile
#'
#' @examples
#'
#' analysis.concatenationTree(alignment.file = "concatenated/my_concat.phy",
#'                            output.directory = "concatenation-trees",
#'                            output.name = "my_concat",
#'                            partition.scheme = "merge",
#'                            uf.bootstrap = 1000,
#'                            threads = 8)
#'
#' #Inside a jackknife or other resampling procedure, where the resampling is
#' #the replication and a per-replicate bootstrap would only add runtime:
#'
#' analysis.concatenationTree(alignment.file = "replicates/rep_0001.phy",
#'                            output.directory = "replicate-trees",
#'                            output.name = "rep_0001",
#'                            partition.scheme = "none",
#'                            model = "GTR+G",
#'                            uf.bootstrap = 0,
#'                            threads = 4)
#'
#' @export

analysis.concatenationTree = function(alignment.file = NULL,
                             output.directory = NULL,
                             output.name = NULL,
                             partition.file = NULL,
                             partition.scheme = c("file", "merge", "none"),
                             model = "GTR",
                             codon.partition = FALSE,
                             seq.type = NULL,
                             program = "IQTREE",
                             msub.type = c("nuclear", "mitochondrial"),
                             uf.bootstrap = 1000,
                             rcluster = 100,
                             threads = 1,
                             memory = 1,
                             iqtree.path = NULL,
                             seed = NULL,
                             quiet = FALSE,
                             resume = TRUE,
                             overwrite = FALSE) {

  #match.arg picks the first choice when the argument is left at its default.
  #Without it the default is the whole vector, and comparing it with == is an
  #error in R 4.2 and later, so the function could not be called on defaults.
  partition.scheme = match.arg(partition.scheme)
  msub.type = match.arg(msub.type)

  if (is.null(alignment.file) || file.exists(alignment.file) == FALSE){ stop("A valid alignment.file is needed.") }
  if (is.null(output.directory)){ stop("An output.directory is needed.") }
  if (is.null(output.name)){ stop("An output.name is needed.") }
  if (length(uf.bootstrap) != 1 || is.numeric(uf.bootstrap) == FALSE || uf.bootstrap < 0){ stop("uf.bootstrap must be 0 or greater.") }
  if (uf.bootstrap > 0 && uf.bootstrap < 1000){ stop("IQ-TREE requires at least 1000 UFBoot replicates. Use uf.bootstrap = 0 to disable bootstrapping.") }
  if (length(memory) != 1 || is.numeric(memory) == FALSE || memory <= 0){ stop("memory must be greater than 0 GB.") }
  if (is.null(seed) == FALSE && (length(seed) != 1 || is.numeric(seed) == FALSE || seed < 1)){ stop("seed must be NULL or a positive number.") }
  if (partition.scheme == "file" && (is.null(partition.file) || file.exists(partition.file) == FALSE)){
    stop("partition.scheme = 'file' needs a partition.file that exists.")
  }
  if (resume == TRUE && overwrite == TRUE){
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }
  if (identical(normalizePath(alignment.file, winslash = "/", mustWork = FALSE),
                normalizePath(output.directory, winslash = "/", mustWork = FALSE))){
    stop("You should not overwrite the original alignments.")
  }

  if (is.numeric(threads) && length(threads) == 1 && threads >= 1){
    threads = as.character(as.integer(threads))
  } else if (length(threads) == 1 && toupper(as.character(threads)) == "AUTO"){
    threads = "AUTO"
  } else { stop("threads must be a positive number or 'AUTO'.") }

  #Finds IQ-TREE, whether it is installed as iqtree2 (version 2) or iqtree (version 3)
  iqtree = findIQTREE(iqtree.path = iqtree.path, quiet = quiet)

  if (dir.exists(output.directory) == TRUE){
    if (overwrite == TRUE){
      unlink(output.directory, recursive = TRUE)
      dir.create(output.directory, recursive = TRUE)
    }
  } else { dir.create(output.directory, recursive = TRUE) }

  run.dir = file.path(output.directory, output.name)
  dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)

  run.alignment = file.path(run.dir, "alignment.phy")
  tree.file = file.path(run.dir, paste0(output.name, ".treefile"))

  if (resume == TRUE && file.exists(tree.file) == TRUE){
    if (quiet == FALSE){ print(paste0(output.name, " already has a treefile, skipping.")) }
    return(invisible(tree.file))
  }

  file.copy(alignment.file, run.alignment, overwrite = TRUE)

  #Sets up parameter type selections from above
  iqtree.args = c("-s", run.alignment,
                  "-pre", file.path(run.dir, output.name),
                  "-nt", threads,
                  "-mem", paste0(memory, "G"))

  #-msub restricts ModelFinder's amino-acid model set, so it only means anything
  #when ModelFinder runs. With an explicit model (partition.scheme = "none") it
  #is a no-op on IQ-TREE 2 and some IQ-TREE builds reject it with a sequence-type
  #error, so add it only for the ModelFinder schemes.
  if (partition.scheme == "none"){
    iqtree.args = c(iqtree.args, "-m", model)
  } else {
    part.scheme = "MFP"
    if (partition.scheme == "merge"){
      part.scheme = paste0(part.scheme, "+MERGE")
      iqtree.args = c(iqtree.args, "-rcluster", rcluster)
    }
    if (partition.scheme == "file"){ iqtree.args = c(iqtree.args, "-spp", partition.file) }
    iqtree.args = c(iqtree.args, "-m", part.scheme, "-msub", msub.type)
  }

  #uf.bootstrap = 0 leaves -bb off entirely. IQ-TREE rejects -bb below 1000, so
  #passing a small number here used to make the run fail rather than run faster.
  if (uf.bootstrap > 0){ iqtree.args = c(iqtree.args, "-bb", as.integer(uf.bootstrap)) }

  #Sequence type. seq.type forces the data type with -st instead of leaving it to
  #IQ-TREE's auto-detection, which some builds fail with "Unknown sequence type"
  #on matrices carrying a lot of missing data. codon.partition is kept as a
  #shortcut for -st CODON. If both are given, seq.type wins.
  st.value = seq.type
  if (codon.partition == TRUE && is.null(st.value)){ st.value = "CODON" }
  if (is.null(st.value) == FALSE && nzchar(st.value)){ iqtree.args = c(iqtree.args, "-st", st.value) }

  if (is.null(seed) == FALSE){ iqtree.args = c(iqtree.args, "-seed", as.integer(seed)) }
  if (quiet == TRUE){ iqtree.args = c(iqtree.args, "-quiet") }
  if (resume == FALSE){ iqtree.args = c(iqtree.args, "-redo") }

  #Runs IQTree
  command = paste(shQuote(c(iqtree$path, iqtree.args)), collapse = " ")
  if (quiet == FALSE){ cat(paste0(command, "\n")) }
  run.status = system(command)

  if (run.status != 0){ stop("IQ-TREE exited with status ", run.status, " for ", output.name) }
  if (file.exists(tree.file) == FALSE){ stop("IQ-TREE produced no treefile for ", output.name) }

  print(paste0(output.name, " finished concatenation tree estimation!"))
  return(invisible(tree.file))

}#end function
