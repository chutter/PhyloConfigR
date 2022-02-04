#' @title analysis.concatenationTree
#'
#' @description Function for concatenating a large number of alignments
#'
#' @param alignment.folder folder that contains aligmnents to be concatenated
#'
#' @param output.name output file name
#'
#' @param partition.file TRUE to save a partition file
#'
#' @param output.format output file format. "all" to save all three types phylip, nexus, and fasta
#'
#' @param partition.format partition file format. "all" to save both raxml and table format
#'
#' @param overwrite TRUE to overwrite file. Default FALSE.
#'
#' @return saves to file concatenated alignments and partition files delimiting the coordinates of each indidividual marker
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
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

