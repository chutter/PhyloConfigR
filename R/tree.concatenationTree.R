#' @title analysis.concatenationTree
#'
#' @description Runs IQ-TREE 2 on a concatenated alignment to estimate a
#'   maximum-likelihood species tree. Supports MFP model selection with optional
#'   partition merging, codon partitioning, and UFBoot branch support. The
#'   alignment is copied into a per-run subdirectory inside output.directory.
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
#'   partition file; "merge" uses MFP+MERGE to find optimal merging;
#'   "none" fits a single GTR model
#'
#' @param codon.partition if TRUE adds -st CODON flag for codon-aware model
#'   fitting (requires in-frame codon alignment)
#'
#' @param program reserved for future use; currently only "IQTREE" is supported
#'
#' @param msub.type substitution model category passed to IQ-TREE -msub flag;
#'   "nuclear" or "mitochondrial"
#'
#' @param uf.bootstrap number of ultrafast bootstrap replicates (default: 100)
#'
#' @param rcluster percentage of partitions used in the rcluster algorithm for
#'   partition model selection (default: 100)
#'
#' @param threads number of CPU threads passed to IQ-TREE -nt flag
#'
#' @param memory memory in GB (currently informational)
#'
#' @param iqtree.path path to the directory containing the iqtree2 executable,
#'   or NULL if iqtree2 is on the system PATH
#'
#' @param resume if TRUE allows IQ-TREE to resume an interrupted run
#'
#' @param overwrite if TRUE removes the existing output directory before running
#'
#' @return IQ-TREE output files are written to output.directory/output.name/;
#'   nothing is returned in R
#'
#' @examples
#'
#' analysis.concatenationTree(alignment.file = "concat_alignment.phy",
#'                             output.directory = "concat-trees",
#'                             output.name = "all-markers",
#'                             partition.scheme = "merge",
#'                             uf.bootstrap = 1000,
#'                             threads = 4)
#'
#' @export

analysis.concatenationTree = function(alignment.file = NULL,
                                      output.directory = NULL,
                                      output.name = NULL,
                                      partition.file = NULL,
                                      partition.scheme = c("file", "merge", "none"),
                                      codon.partition = FALSE,
                                      program = "IQTREE",
                                      msub.type = c("mitochondrial", "nuclear"),
                                      uf.bootstrap = 100,
                                      rcluster = 100,
                                      threads = 1,
                                      memory = 1,
                                      iqtree.path = NULL,
                                      resume = TRUE,
                                      overwrite = FALSE) {

  #Debug
  # alignment.file = alignment.files[i]
  # output.directory = out.path
  # output.name = align.name
  # partition.file = NULL
  # partition.scheme = "merge"
  # codon.partition = FALSE
  # program = "IQTREE"
  # msub.type = "nuclear"
  # uf.bootstrap = uf.bootstrap
  # rcluster = rcluster
  # threads = threads
  # memory = memory
  # iqtree.path = iqtree.path
  # resume = resume
  # overwrite = overwrite

  #Checks and formats path
  if (is.null(iqtree.path) == FALSE){
    b.string = unlist(strsplit(iqtree.path, ""))
    if (b.string[length(b.string)] != "/") {
      iqtree.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { iqtree.path = NULL }

  if (alignment.file == output.directory){ stop("You should not overwrite the original alignments.") }

  # if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }

  #So I don't accidentally delete everything while testing resume
  if (resume == TRUE & overwrite == TRUE){
    overwrite = FALSE
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }
#
#   #Gathers alignments
#   iq.files = list.files(alignment.dir)
#
#   if (length(align.files) == 0) { stop("alignment files could not be found.") }
#
#   #Skips files done already if resume = TRUE
#   if (resume == TRUE){
#     done.files = list.files(output.dir)
#     align.files = align.files[!gsub("\\..*", "", align.files) %in% gsub("\\..*", "", done.files)]
#   }

  #Sets up parameter type selections from above
  part.file = ""
  part.scheme = "MFP"
  if (partition.scheme == "merge"){ part.scheme = paste0(part.scheme, "+MERGE") }
  if (partition.scheme == "file"){ part.file = paste0(" -spp ", partition.file) }
  if (partition.scheme == "none"){ part.scheme = "GTR" }
  if (codon.partition == T){ codon.st = " -st CODON" } else { codon.st = "" }

  dir.create(paste0(output.directory, "/", output.name))
  system(paste0("cp ", alignment.file, " ", output.directory, "/", output.name, "/alignment.phy"))

  #Runs IQTree
  system(paste0(iqtree.path, "iqtree2 -s ", output.directory, "/", output.name, "/alignment.phy", part.file,
                " -bb ", uf.bootstrap,
                " -nt ", threads,
                " -m ", part.scheme, codon.st,
                " -rcluster ", rcluster,
                " -msub ", msub.type))

  print(paste0(output.name, " finished concatenation tree estimation!"))

}#end function

